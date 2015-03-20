#! /usr/bin/env python
__author__ = 'bendall'

# Initial author: Solaiappan Manimaran
# Functions to read alignment file (sam/gnu-sam or bl8), run EM algorithm
# and output report file that can be opened in Excel and also
# output updated alignment file (sam/gnu-sam or bl8)

#	Pathoscope - Predicts strains of genomes in Nextgen seq alignment file (sam/bl8)
#	Copyright (C) 2013  Johnson Lab - Boston University
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, sys
from time import time
import pysam
from pathoscope.pathorna.utils import PSRead
from pathoscope.pathorna.utils import AnnotationLookup
from pathoscope.pathorna.utils import iterread
from pathoscope.pathorna.utils import phred # Used by: updated_alignments
from pathoscope.utils import samUtils


class PathoRNAOptions:
  option_fields = ['verbose','score_cutoff','out_matrix','no_updated_sam','emEpsilon','maxIter','piPrior','thetaPrior',
                   'exp_tag','outdir','ali_format','ali_file','gtf_file', 'out_samfile']

  def __init__(self,alnfile,gtffile,**kwargs):
    # Input options
    self.ali_format      = 'sam'
    self.ali_file        = os.path.abspath(alnfile)
    self.gtf_file        = os.path.abspath(gtffile)
    self.no_feature_key  = '__nofeature__'

    # Output options
    self.verbose         = False
    self.outdir          = '.'
    self.out_samfile     = None
    self.exp_tag         = 'pathorna'
    self.report_all      = False
    self.min_final_guess = 0.01
    self.out_matrix      = False
    self.out_abundance   = False
    self.no_updated_sam  = False

    # Model parameters
    self.piPrior         = 0
    self.thetaPrior      = 0
    self.score_cutoff    = 0.01

    # EM parameters
    self.emEpsilon       = 1e-7
    self.maxIter         = 100

    # Set any other options passed in kwargs
    for k,v in kwargs.iteritems():
      if v is None: continue
      assert hasattr(self,k)
      if getattr(self,k) != v: setattr(self,k,v)

    self.outdir = os.path.abspath(self.outdir)

    if self.out_samfile is None:
      self.out_samfile = self.generate_filename('updated.sam')

  def generate_filename(self,suffix):
    basename = '%s-%s' % (self.exp_tag, suffix)
    return os.path.join(self.outdir,basename)

  def __str__(self):
    _ret = "PathoRNAOptions:\n"
    _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
    return _ret

"""
Wrappers for pathoscope functions
"""

def wrap_pathoscope_em(unique,repeat,genomes,opts):
  from pathoscope.pathoid import PathoID
  return PathoID.pathoscope_em(unique, repeat, genomes, opts.maxIter, opts.emEpsilon, opts.verbose, opts.piPrior, opts.thetaPrior)

def wrap_computeBestHit(unique,repeat,genomes,reads):
  from pathoscope.pathoreport import PathoReportA
  tup = PathoReportA.computeBestHit(unique, repeat, genomes, reads)
  labels = ['bestHitReads','bestHit','level1','level2']
  return dict(zip(labels,tup))

"""
Functions for parsing SAM file
"""
def mp_parse_reads(samfile, flookup, opts=None):
  ''' Multiprocessing version of read parser  '''
  import pathos.multiprocessing as mp

  def _procread(rname,segments):
    r = PSRead(rname,segments)
    if r.is_unmapped:
      return None
    r.assign_feats(refnames, flookup)
    if not r.aligns_to_feat():
      return None
    r.assign_best()
    return r

  def _collect(result):
    if result is not None:
      mapped[result.readname] = result

  _verbose = opts.verbose if opts is not None else True

  refnames = dict(enumerate(samfile.references))
  # counts = {'unmapped':0, 'nofeat':0, 'mapped':0}
  mapped   = {}

  pool = mp.Pool(processes=4)
  for t in iterread(samfile):
    pool.apply_async(_procread, t, callback=_collect)

  pool.close()
  print >>sys.stderr, "pool closed"
  pool.join()
  print >>sys.stderr, "pool join"
  return mapped

def sp_parse_reads(samfile, flookup, opts=None):
  _verbose = opts.verbose if opts is not None else True

  counts = {'unmapped':0, 'nofeat':0, 'mapped':0}
  mapped   = {}

  # Lookup reference name from reference ID
  refnames = dict(enumerate(samfile.references))

  for rname,segments in iterread(samfile):
    r = PSRead(rname,segments)
    if r.is_unmapped:
      counts['unmapped'] += 1
    else:
      r.assign_feats(refnames, flookup)
      if r.aligns_to_feat():
        r.assign_best()
        mapped[rname] = r
        counts['mapped'] += 1
      else:
        counts['nofeat'] += 1

    if _verbose and sum(counts.values()) % 100000 == 0:
      print >>sys.stderr, "...Processed %d fragments" % sum(counts.values())

  if _verbose:
    print >>sys.stderr, "Processed %d fragments" % sum(counts.values())
    print >>sys.stderr, "\t%d fragments were unmapped" % counts['unmapped']
    print >>sys.stderr, "\t%d fragments mapped to one or more positions on reference genome" % (counts['mapped'] + counts['nofeat'])
    print >>sys.stderr, "\t\t%d fragments mapped to reference but did not map to annotation" % counts['nofeat']
    print >>sys.stderr, "\t\t%d fragments have at least one alignment within annotation" % counts['mapped']

  return mapped

def parse_reads(samfile, flookup, opts=None):
  # Multiprocessing version does not perform well
  # Performance is worse than single processing version
  """
  try:
    import pathos.multiprocessing as mp
    print >>sys.stderr, "Using multiprocessing..."
    return mp_parse_reads(samfile, flookup, opts=None)
  except ImportError:
    print >>sys.stderr, "Using single processing..."
    return sp_parse_reads(samfile, flookup, opts=None)
  """
  return sp_parse_reads(samfile, flookup, opts=None)


def data_matrix(reads):
  '''
  Creates data structure used by pathoscope_em and other functions
  :param reads:
  :return:
  '''
  _unique = {}
  _repeat = {}
  maxscore = float('-inf')
  minscore = float('inf')
  gset = set()
  for rname,r in reads.iteritems():
    if r.is_unmapped: continue
    d = r.structured_data()
    # Ignore if read did not align to any annotations
    if all(g==PSRead.nofeature for g in d[0]):
      continue
    gset.update(set(d[0]))
    maxscore = max(d[3],maxscore)
    minscore = min(d[1]+[minscore])
    if r.unique_feat():
      _unique[rname] = d
    else:
      _repeat[rname] = d

  # Rescale alignment scores
  _ = samUtils.rescale_samscore(_unique, _repeat, maxscore, minscore)

  # Restructure unique to contain only genome index and score
  for rname in _unique.keys():
    _unique[rname] = [_unique[rname][0][0], _unique[rname][1][0]]

  # Normalize scores in repeat
  for rname in _repeat.keys():
    pScoreSum = sum(_repeat[rname][1])
    _repeat[rname][2] = [k/pScoreSum for k in _repeat[rname][1]]

  # Set genomes to integers
  _genomes = list(gset)
  gdict   = dict((v,i) for i,v in enumerate(_genomes))
  for rname,data in _unique.iteritems():
    data[0] = gdict[data[0]]

  for rname,data in _repeat.iteritems():
    data[0] = [gdict[g] for g in data[0]]

  # Set readnames to integers
  _reads = _unique.keys() + _repeat.keys()
  for i,rn in enumerate(_reads):
    if rn in _unique:
      _unique[i] = _unique.pop(rn)
    else:
      _repeat[i] = _repeat.pop(rn)

  return _unique, _repeat, _genomes, _reads

"""
Tags used in updated alignments:
ZF:Z - Name of annotation feature containing alignment
ZQ:i - Quality score
ZG
"""

def updated_alignments(psread, rdata, glookup, score_cutoff):
  '''

  :param psread:           PSRead to get alignments for
  :param rdata:            Matrix data for this read
  :param glookup:          Lookup table to get feature name
  :param score_cutoff:     Minimum score to be included in alignment
  :return:
  '''

  if len(rdata) == 2:                                                          # Initially mapped to only one feature
    gname = glookup[rdata[0]]
    pri_aln,alt_alns = psread.aligned_to_genome(gname)
    # Primary alignment was the one chosen by PSRead.assign_best and was used in PathoID
    pri_aln.set_tags({'ZP':'UP', 'ZQ':1.0,}).set_mapq(255)                     # Unique Primary, set mapq to 255
    if gname != PSRead.nofeature: pri_aln.set_tags('ZF',gname)
    # Alternate alignments are present if read mapped multiple times to same feature
    for a in alt_alns:
      a.set_tags({'ZP':'UA', 'ZQ':0.0, 'ZF':gname,}).set_mapq(0)               # Unique Alternate, set mapq to 0
      if gname != PSRead.nofeature: a.set_tags('ZF',gname)
    return [pri_aln] + alt_alns
  else:                                                                        # Initially mapped to more than one feature
    _updated = []
    top_pscore = max(rdata[2])
    top_genome = [glookup[rdata[0][i]] for i,s in enumerate(rdata[2]) if s==top_pscore]
    if len(top_genome)==1:
      top_genome = top_genome[0]
    else:
      top_genome = 'MULT'

    # Iterate over updated Pscores
    for i,upPscore in enumerate(rdata[2]):
      if upPscore > score_cutoff:
        gname = glookup[rdata[0][i]]
        pri_aln,alt_alns = psread.aligned_to_genome(gname)
        level = 'H' if upPscore >= 0.5 else 'L'
        pri_aln.set_tags({'ZP':'%sP' % level, 'ZQ':round(upPscore,3),}).set_mapq(phred(upPscore)).set_tags('ZG',top_genome)
        if gname != PSRead.nofeature: pri_aln.set_tags('ZF',gname)
        for a in alt_alns:
          a.set_tags({'ZP':'%sA' % level, 'ZQ':0.0,}).set_mapq(0)
          if gname != PSRead.nofeature: a.set_tags('ZF',gname)
        _updated += [pri_aln] + alt_alns
    return _updated

def write_tsv_report(genomes, initial_guess, final_guess, initial_report, final_report, nreads, opts):
  header1 = ['Total Number of Aligned Reads:', str(nreads), 'Total Number of Mapped Genomes:', str(len(genomes))]
  header2 = ['Genome', 'Final Guess', 'Final Best Hit', 'Final Best Hit Read Numbers',
             'Final High Confidence Hits', 'Final Low Confidence Hits', 'Initial Guess',
             'Initial Best Hit', 'Initial Best Hit Read Numbers',
             'Initial High Confidence Hits', 'Initial Low Confidence Hits']
  # The order in "report" must be the same as in "header2"
  report = zip(genomes,
               final_guess, final_report['bestHit'], final_report['bestHitReads'],
               final_report['level1'], final_report['level2'],
               initial_guess, initial_report['bestHit'], initial_report['bestHitReads'],
               initial_report['level1'], initial_report['level2']
              )

  # Sort report by final_guess
  report.sort(key=lambda x:x[1],reverse=True)

  # Only include genomes where Final guess is >= min_final_guess or have final hits > 0
  if not opts.report_all:
    report = [r for r in report if r[1] >= opts.min_final_guess or r[4] > 0 or r[5] > 0]

  with open(opts.generate_filename('report.tsv'),'w') as outh:
    print >>outh, '\t'.join(header1)
    print >>outh, '\t'.join(header2)
    for r in report:
      print >>outh, '\t'.join(str(_) for _ in r)

def write_abundance_report(genomes, initial_guess, final_guess, initial_report, final_report, nreads, opts,
                           genome_lengths, avg_read_len):
  '''
  :param genomes:
  :param initial_guess:
  :param final_guess:
  :param initial_report:
  :param final_report:
  :param nreads:
  :param opts:
  :param genome_lengths:
  :param avg_read_len:
  :param reportAll:
  :return:
  '''
  init_rpkm  = calculate_rpkm(genomes, genome_lengths, initial_report['bestHitReads'],nreads)
  init_tpm   = calculate_tpm(genomes, genome_lengths, initial_report['bestHitReads'], avg_read_len)
  final_rpkm = calculate_rpkm(genomes, genome_lengths, final_report['bestHitReads'],nreads)
  final_tpm  = calculate_tpm(genomes, genome_lengths, final_report['bestHitReads'], avg_read_len)
  glen_list  = [genome_lengths[g] for g in genomes]
  header1 = ['Total Number of Aligned Reads:', str(nreads), 'Total Number of Mapped Genomes:', str(len(genomes))]
  header2 = ['Genome', 'Length',
             'Final Best Hit Read Numbers',
             'Final RPKM','Final TPM',
             'Initial Best Hit Read Numbers',
             'Initial RPKM','Initial TPM',
             'Final Guess', 'Final Best Hit', 'Final High Confidence Hits', 'Final Low Confidence Hits',
             'Initial Guess','Initial Best Hit', 'Initial High Confidence Hits', 'Initial Low Confidence Hits',
             ]
  # The order in "report" must be the same as in "header2"
  report = zip(genomes, glen_list,
               final_report['bestHitReads'],
               final_rpkm, final_tpm,
               initial_report['bestHitReads'],
               init_rpkm, init_tpm,
               final_guess, final_report['bestHit'], final_report['level1'], final_report['level2'],
               initial_guess, initial_report['bestHit'], initial_report['level1'], initial_report['level2'],
              )

  # Sort report by final_rpkm
  report.sort(key=lambda x:x[3],reverse=True)

  # Optional filtering of genomes
  #if not reportAll:
  #  pass # filtered = [r for r in report if r[8] >= opts.score_cutoff or r[10] > 0 or r[11] > 0]

  with open(opts.generate_filename('abundance.tsv'),'w') as outh:
    print >>outh, '\t'.join(header1)
    print >>outh, '\t'.join(header2)
    for r in report:
      print >>outh, '\t'.join(str(_) for _ in r)

def average_read_length(rdict,rlist):
  numer = sum(max(a.query_length for a in rdict[r].alignments) for r in rlist)
  return int(round(float(numer) / len(rlist)))

def calculate_rpkm(features, feat_len, counts, nreads):
  rpkms = []
  for f,rc in zip(features,counts):
    rpkm = (rc * 1e9) / (feat_len[f] * nreads)
    rpkms.append(rpkm)
  return rpkms

def calculate_tpm(features, feat_len, counts, r_l):
  T = sum(((rc * r_l) / feat_len[f]) for f,rc in zip(features,counts))
  tpms = []
  for f,rc in zip(features,counts):
    tpm = (rc * r_l * 1e6) / (feat_len[f] * T)
    tpms.append(tpm)
  return tpms

def pathoscope_rna_reassign(opts):


  PSRead.nofeature = opts.no_feature_key

  flookup = AnnotationLookup(opts.gtf_file)
  samfile = pysam.AlignmentFile(opts.ali_file)

  if opts.verbose:
    print >>sys.stderr, "Loading alignment file (%s)" % opts.ali_file
    loadstart = time()

  allreads = parse_reads(samfile, flookup)

  if opts.verbose:
    print >>sys.stderr, "Time to load alignment:".ljust(40) + "%d seconds" % (time() - loadstart)

  U,NU,genomes,reads = data_matrix(allreads)

  if False: # Make a deep copy of the unique and non-unique data structures
    import copy
    initU = copy.deepcopy(U)
    initNU = copy.deepcopy(NU)

  initial_report = wrap_computeBestHit(U, NU, genomes, reads)

  if opts.verbose:
    print >>sys.stderr, "EM iteration..."
    print >>sys.stderr, "(Genomes,Reads)=%dx%d" % (len(genomes),len(reads))
    print >>sys.stderr, "Delta Change:"
    emtime = time()

  (initPi, pi, _, NU) = wrap_pathoscope_em(U, NU, genomes, opts)

  if opts.verbose:
    print >>sys.stderr, "Time for EM iteration:".ljust(40) +  "%d seconds" % (time() - emtime)

  if opts.out_matrix:
    with open(opts.generate_filename('genomeId.txt'),'w') as outh:
      print >>outh, '\n'.join(genomes)
    with open(opts.generate_filename('readId.txt'),'w') as outh:
      print >>outh, '\n'.join(reads)
    with open(opts.generate_filename('initGuess.txt'),'w') as outh:
      for p,g in sorted(zip(initPi,genomes), key=lambda x:x[0], reverse=True):
        print >>outh, '%.7g\t%s' % (p,g)
    with open(opts.generate_filename('finGuess.txt'),'w') as outh:
      for p,g in sorted(zip(pi,genomes), key=lambda x:x[0], reverse=True):
        print >>outh, '%.7g\t%s' % (p,g)

  final_report = wrap_computeBestHit(U, NU, genomes, reads)

  write_tsv_report(genomes, initPi, pi, initial_report, final_report, len(reads), opts)

  if opts.out_abundance: # Output abundance measures
    feature_lengths = flookup.feature_length()
    # Calculate the length of the "no feature"
    # nofeat_length = sum(samfile.lengths) - sum(feature_lengths.values())
    feature_lengths[opts.no_feature_key] = sum(samfile.lengths) - sum(feature_lengths.values()) #nofeat_length
    # Calculate average read length
    avg_rlen = average_read_length(allreads,reads)
    # abundance_list = calculate_abundance(genomes, feature_lengths, final_report['bestHitReads'], len(reads), avg_rlen)
    write_abundance_report(genomes, initPi, pi, initial_report, final_report, len(reads), opts, feature_lengths, avg_rlen)

  # Write the updated sam file
  if not opts.no_updated_sam:
    updated_samfile = pysam.AlignmentFile(opts.out_samfile, 'wh', header=samfile.header)
    glookup = dict(enumerate(genomes)) # Lookup genome name by index
    for ridx,rname in enumerate(reads):
      rdata = U[ridx] if ridx in U else NU[ridx]
      u_alns = updated_alignments(allreads[rname],rdata,glookup,opts.score_cutoff)
      for aln in u_alns:
        aln.write_samfile(updated_samfile)

  samfile.close()
  if not opts.no_updated_sam: updated_samfile.close()
  return
