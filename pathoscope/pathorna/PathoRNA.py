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
import pysam
from pathoscope.pathorna.utils import PSRead, FeatureLookup
from pathoscope.pathorna.utils import iterread
from pathoscope.utils import samUtils

class PathoRNAOptions:
  option_fields = ['verbose','score_cutoff','out_matrix','no_updated_sam','emEpsilon','maxIter','piPrior','thetaPrior',
                   'exp_tag','outdir','ali_format','ali_file','gtf_file', 'out_samfile']

  def __init__(self,alnfile,gtffile,**kwargs):
    # Default option values
    self.verbose         = False
    self.score_cutoff    = 0.01
    self.out_matrix      = False
    self.no_updated_sam  = False
    self.emEpsilon       = 1e-7
    self.maxIter         = 50
    self.piPrior         = 0
    self.thetaPrior      = 0
    self.exp_tag         = 'pathorna'
    self.outdir          = '.'
    self.ali_format      = 'sam'
    self.out_samfile     = None
    self.ali_file        = os.path.abspath(alnfile)
    self.gtf_file        = os.path.abspath(gtffile)

    # Set any other options passed in kwargs
    for k,v in kwargs.iteritems():
      if v is None: continue
      assert hasattr(self,k)
      if getattr(self,k) != v: setattr(self,k,v)

    self.outdir = os.path.abspath(self.outdir)

    if self.out_samfile is None:
      self.out_samfile = self.generate_filename('updated.sam')

    """
    # If outfile was not specified, set path for outfile
    if self.out_samfile is None:
      basename = os.path.split(self.ali_file)[1]
      fn,ext = os.path.splitext(basename)
      if ext == '.sam':
        self.out_samfile = os.path.join(self.outdir,'%s.updated.sam' % fn)
      else:
        self.out_samfile = os.path.join(self.outdir,'%s.updated.sam' % fn)
    """

  def generate_filename(self,suffix):
    basename = '%s-%s-%s' % (self.exp_tag, self.ali_format, suffix)
    return os.path.join(self.outdir,basename)

  def __str__(self):
    _ret = "PathoRNAOptions:\n"
    _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
    return _ret

def wrap_pathoscope_em(unique,repeat,genomes,opts):
  from pathoscope.pathoid import PathoID
  return PathoID.pathoscope_em(unique, repeat, genomes, opts.maxIter, opts.emEpsilon, opts.verbose, opts.piPrior, opts.thetaPrior)

def wrap_computeBestHit(unique,repeat,genomes,reads):
  from pathoscope.pathoreport import PathoReportA
  tup = PathoReportA.computeBestHit(unique, repeat, genomes, reads)
  labels = ['bestHitReads','bestHit','level1','level2']
  return dict(zip(labels,tup))

def parse_reads(samfile, flookup):
  from pathoscope.pathorna.utils import PSRead
  from pathoscope.pathorna.utils import iterread

  allreads = {}

  # Lookup reference name from reference ID
  refnames = dict(enumerate(samfile.references))

  for rname,segments in iterread(samfile):
    allreads[rname] = PSRead(rname,segments)
    allreads[rname].assign_feats(refnames, flookup)
    allreads[rname].assign_best()
    # if len(allreads) >= 2000: break

  return allreads

def data_matrix(reads):
  from pathoscope.utils import samUtils
  _unique = {}
  _repeat = {}
  maxscore = float('-inf')
  minscore = float('inf')
  gset = set()
  for rname,r in reads.iteritems():
    if r.is_unmapped: continue
    d = r.structured_data()
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

def updated_alignments(psread,rdata,glookup,score_cutoff):
  from pathoscope.pathorna.utils import phred

  if len(rdata) == 2: # This is a uniquely mapped read
    gname = glookup[rdata[0]]
    pri_aln,alt_alns = psread.aligned_to_genome(gname)
    pri_aln.set_tags('ZP','UP').set_tags('ZQ',255)     # Unique Primary
    for a in alt_alns:
      a.set_tags('ZP','UA').set_tags('ZQ',255)      # Unique Alternate
    return [pri_aln] + alt_alns
  else: # This is a non-uniquely mapped read
    _updated = []
    # Iterate over updated Pscores
    for i,upPscore in enumerate(rdata[2]):
      if upPscore > score_cutoff:
        gname = glookup[rdata[0][i]]
        pri_aln,alt_alns = psread.aligned_to_genome(gname)
        level = 'H' if upPscore >= 0.5 else 'L'
        pri_aln.set_tags('ZP','%sP' % level).set_tags('ZQ',phred(upPscore))
        for a in alt_alns:
          a.set_tags('ZP','%sA' % level).set_tags('ZQ',phred(upPscore))
        _updated += [pri_aln] + alt_alns
    return _updated

def write_tsv_report(genomes, initial_guess, final_guess, initial_report, final_report, nreads, opts):
  header1 = ['Total Number of Aligned Reads:', str(nreads), 'Total Number of Mapped Genomes:', str(len(genomes))]
  header2 = ['Genome', 'Final Guess', 'Final Best Hit', 'Final Best Hit Read Numbers',
             'Final High Confidence Hits', 'Final Low Confidence Hits', 'Initial Guess',
             'Initial Best Hit', 'Initial Best Hit Read Numbers',
             'Initial High Confidence Hits', 'Initial Low Confidence Hits']
  # Ensure that the order in "report" is the same as "header2"
  report = zip(genomes,
               final_guess, final_report['bestHit'], final_report['bestHitReads'],
               final_report['level1'], final_report['level2'],
               initial_guess, initial_report['bestHit'], initial_report['bestHitReads'],
               initial_report['level1'], initial_report['level2']
              )

  # Sort report by final_guess
  report.sort(key=lambda x:x[1],reverse=True)

  # Only include genomes where Final guess is >= score_cutoff or have final hits > 0
  filtered = [r for r in report if r[1] >= opts.score_cutoff or r[4] > 0 or r[4] > 0]

  with open(opts.generate_filename('report.tsv'),'w') as outh:
    print >>outh, '\t'.join(header1)
    print >>outh, '\t'.join(header2)
    for r in filtered:
      print >>outh, '\t'.join(str(_) for _ in r)



def pathoscope_rna_reassign(opts):
  import pysam
  from pathoscope.pathorna.utils import FeatureLookup
  from time import time

  flookup = FeatureLookup(opts.gtf_file)
  samfile = pysam.AlignmentFile(opts.ali_file)

  if opts.verbose:
    print >>sys.stderr, "Loading alignment file (%s)" % opts.ali_file
    loadstart = time()

  allreads = parse_reads(samfile, flookup)

  if opts.verbose:
    print >>sys.stderr, "Time to load alignment:".ljust(40) + "%d seconds" % (time() - loadstart)

  U,NU,genomes,reads = data_matrix(allreads)

  if False:
    import copy
    initU = copy.deepcopy(U)
    initNU = copy.deepcopy(NU)

  initial_report = wrap_computeBestHit(U, NU, genomes, reads)

  if opts.verbose:
    print >>sys.stderr, "EM iteration..."
    print >>sys.stderr, "(Genomes,Reads)=%dx%d" % (len(genomes),len(reads))
    print >>sys.stderr, "Delta Change:"
    emtime = time()\

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

  # Write the updated sam file
  if not opts.no_updated_sam:
    updated_samfile = pysam.AlignmentFile(opts.out_samfile,'wh',header=samfile.header)
    glookup = dict(enumerate(genomes))
    for ridx,rname in enumerate(reads):
      rdata = U[ridx] if ridx in U else NU[ridx]
      u_alns = updated_alignments(allreads[rname],rdata,glookup,opts.score_cutoff)
      for aln in u_alns:
        aln.write_samfile(updated_samfile)

  samfile.close()
  if not opts.no_updated_sam: updated_samfile.close()
  return
