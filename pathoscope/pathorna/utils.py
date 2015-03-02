__author__ = 'bendall'

import pysam
import math
from collections import defaultdict
import re

"""
SAM flags:
  0x1   template having multiple segments in sequencing
  0x2   each segment properly aligned according to the aligner
  0x4   segment unmapped
  0x8   next segment in the template unmapped
  0x10  SEQ being reverse complemented
  0x20  SEQ of the next segment in the template being reversed
  0x40  the first segment in the template
  0x80  the last segment in the template
  0x100 secondary alignment
  0x200 not passing quality controls
  0x400 PCR or optical duplicate
  0x800 supplementary alignment
"""

def phred(f):
  """ Calculate phred quality score for given error probability
  :param f: Error probability (float)
  :return:  Phred score (int)
  """
  return int(round(-10 * math.log10(1 - f))) if f < 1.0 else 255

def iterread(samfile):
  """ Iterate over samfile by query name (read ID)
      Each iteration returns all alignments that share the same read ID.
      Assumes that file is sorted by queryname (default output from bowtie2)
  :param samfile: Alignment file (pysam.AlignmentFile)
  :yield:         Read name, alignment list (str, list(pysam.AlignedSegment))
  """
  ralns = [samfile.next()]
  current = ralns[0].query_name
  for aln in samfile:
    if aln.query_name == current:
      ralns.append(aln)
    else:
      yield current,ralns
      ralns = [aln]
      current = aln.query_name
  yield current,ralns

def make_unmapped_mate(mate,template,add_tag=True):
  """ Create mate for read using sequence and quality from template
  :param mate:
  :param template:
  :param add_tag:
  :return:
  """
  a = pysam.AlignedSegment()
  a.query_name           = template.query_name
  a.flag                 = template.flag
  a.reference_id         = mate.reference_id
  a.reference_start      = mate.reference_start
  a.mapping_quality      = 0 # template.mapping_quality
  # a.cigar              = # Not set
  a.next_reference_id    = mate.reference_id
  a.next_reference_start = mate.reference_start
  # a.template_length    = # Not set
  a.query_sequence       = template.query_sequence
  a.query_qualities      = template.query_qualities
  a.tags                 = template.tags
  a.is_secondary         = mate.is_secondary
  a.is_paired      = True
  a.is_proper_pair = False
  a.is_unmapped    = True
  # This tag indicates the segment was created
  if add_tag: a.setTag('ZC',1)
  return a

"""
class Feature:
  def __init__(self, name, chrom, spos, epos, strand=None):
    self.name  = name
    self.chrom = chrom
    self.spos = spos
    self.epos = epos
    self.strand = strand
    self.length = epos - spos
"""

class FeatureLookup:
  def __init__(self,gtffile,attr_name="locus"):
    fh = open(gtffile,'rU') if isinstance(gtffile,str) else gtffile
    lines = (l.strip('\n').split('\t') for l in fh if not l.startswith('#'))
    # List of locus names
    self._locus = []
    # Dictionary with chromosome as key. List of tuples as value
    # Each tuple has a start coordinate, end coordinate, and index
    self._gdict = defaultdict(list)
    for i,l in enumerate(lines):
      attr = dict(re.search('(\S+)\s"(.+?)"',f.strip()).groups() for f in l[8].split(';') if f.strip())
      self._locus.append( attr[attr_name] if attr_name in attr else 'PSRE%04d' % i )
      self._gdict[l[0]].append((int(l[3]),int(l[4]),i))

  def lookup(self,ref,pos,get_index=False):
    ''' Return the feature for a given reference and position '''
    # Tests every locus in chromosome
    feats = [i for s,e,i in self._gdict[ref] if s <= pos <= e]
    if len(feats)==0:
      return None
    else:
      assert len(feats)==1
      if get_index: return feats[0]
      return self._locus[feats[0]]

  def lookup_interval(self,ref,spos,epos):
    ''' Resolve the feature that overlaps or contains the given interval
        NOTE: Only tests the start and end positions. This means that it does not handle
              cases where a feature lies completely within the interval. This is OK when the
              fragment length is expected to be smaller than the feature length. Also does
              not handle cases where the start position and end position lie in different
              features.
    '''
    featL = self.lookup(ref,spos)
    featR = self.lookup(ref,epos)
    if featL is None and featR is None:     # Neither start nor end is within a feature
      return None
    else:
      if featL is None or featR is None:    # One (and only one) end is within a feature
        return featL if featL is not None else featR
      elif featL == featR:                  # Both ends within the same feature
        return featL
      else:
        assert False, "Ends do not agree"

  def feature_length(self):
    _ret = {}
    for chr,ilist in self._gdict.iteritems():
      for spos,epos,locus_idx in ilist:
        _ret[self._locus[locus_idx]] = epos-spos
    return _ret

  def feature_name(self,id):
    return self._locus[id]

class PSAlignment:
  """
  One distinct alignment. For properly paired segments (paired end) this has two segments,
  otherwise has one segment
  """
  def __init__(self, s1, s2=None):
    self.seg1 = s1
    self.seg2 = s2
    if s2 is None:
      self.is_paired      = False
      self.is_proper_pair = False
      self.is_unmapped    = s1.is_unmapped
      self.is_secondary   = s1.is_secondary
      self.query_length   = s1.query_length if not s1.is_unmapped else 0
    else:
      self.is_paired         = True
      self.is_proper_pair = s1.is_proper_pair and s2.is_proper_pair
      self.is_unmapped    = s1.is_unmapped and s2.is_unmapped
      self.is_secondary   = s1.is_secondary and s2.is_secondary
      self.query_length   = sum(s.query_length if not s.is_unmapped else 0 for s in [s1,s2])
      assert s1.is_secondary == s2.is_secondary

    tags = [dict(s.tags) for s in [s1,s2] if s is not None and not s.is_unmapped]
    self.AS = sum(t['AS'] if 'AS' in t else 0 for t in tags)
    self.NM = sum(t['NM'] if 'AS' in t else 0 for t in tags)

  def set_tags(self,tag,value):
    self.seg1.setTag(tag,value)
    if self.is_paired:
      self.seg2.setTag(tag,value)
    # Return instance to allow chaining
    return self

  def set_mapq(self,value):
    # Ensure MAPQ in acceptable range
    if value > 255: value = 255
    if value < 0:   value = 0
    self.seg1.mapping_quality = 0 if self.seg1.is_unmapped else int(value)
    if self.is_paired:
      self.seg2.mapping_quality = 0 if self.seg2.is_unmapped else int(value)
    # Return instance to allow chaining
    return self

  def coordinates(self):
    if self.is_unmapped:
      return None
    if self.is_paired:
      assert self.seg1.reference_id == self.seg2.reference_id
    ref   = self.seg1.reference_id
    left  = min(s.reference_start for s in [self.seg1,self.seg2] if s is not None and not s.is_unmapped)
    right = max(s.reference_end for s in [self.seg1,self.seg2] if s is not None and not s.is_unmapped)
    return ref,left,right

  def write_samfile(self,outfile):
    retval = outfile.write(self.seg1)
    if self.is_paired:
      retval = outfile.write(self.seg2)

  def __str__(self):
    return '%s\n%s\n%s' % (self.AS, self.seg1, self.seg2)

class PSRead:
  """
  All alignments for one read
  """
  nofeature = '__nofeature__'

  def __init__(self,readname,segs):
    self.readname       = readname
    self.alignments     = []
    self.is_paired      = None
    self.is_proper_pair = None
    self.is_unmapped    = None    # Read is unmapped if all alignments are unmapped
    self.is_unique      = None    # Read is unique if only one mapped alignment

    # Assumptions about input file:
    #   - For a given read, segments are either ALL properly paired or ALL not properly paired
    #   - Properly paired segments have their mate on the adjacent line
    paired = segs[0].is_paired
    assert all(s.is_paired == paired for s in segs), "s.is_paired does not match! %s" % [s.is_paired for s in segs]
    if paired:
      self.is_paired = True
      proper = segs[0].is_proper_pair
      assert all(s.is_proper_pair==proper for s in segs), "s.is_proper_pair does not match! %s" % [s.is_proper_pair for s in segs]
      if proper:
        self.is_proper_pair = True
        self.is_unmapped    = False
        self.alignments     = [PSAlignment(segs[i],segs[i+1]) for i in range(0,len(segs),2)]
      else:
        self.is_proper_pair = False
        s0,s1 = segs[:2]
        assert not s0.is_secondary and not s1.is_secondary
        if s0.is_unmapped and s0.mate_is_unmapped:    # Both segments are unmapped
          assert len(segs)==2 and s1.is_unmapped
          self.is_unmapped  = True
          self.alignments   = [PSAlignment(s0,s1)]
        else:
          self.is_unmapped  = False
          if s1.is_unmapped:                          # One segment is unmapped
            assert sum(s.is_unmapped for s in segs) == 1
            self.alignments.append(PSAlignment(s0,s1))
          else:                                       # Both segments are mapped
            assert sum(s.is_unmapped for s in segs) == 0
            self.alignments.append(PSAlignment(s0, make_unmapped_mate(s0,template=s1)))
            s1.is_secondary = True
            self.alignments.append(PSAlignment(s1,make_unmapped_mate(s1,template=s0)))

          if len(segs) > 2:
            lookup_r1 = dict((True if _.is_read1 else False, _) for _ in [s0,s1])
            assert lookup_r1[True].is_read1
            assert lookup_r1[False].is_read2
            for s in segs[2:]:
              assert s.is_secondary
              m = make_unmapped_mate(s,template=lookup_r1[s.is_read2])
              self.alignments.append(PSAlignment(s,m))
    else:                                             # Read is not paired
      self.is_paired      = False
      self.is_proper_pair = False
      self.is_unmapped    = segs[0].is_unmapped
      self.alignments = [PSAlignment(s) for s in segs]

    # Read is unique if the number of mapped alignments == 1
    self.is_unique      = len([a.is_unmapped==False for a in self.alignments]) == 1
    # Initialize features array
    self.features = [None] * len(self.alignments)

  def assign_feats(self, ref_lookup, feature_lookup, use_chrom=False):
    for i,a in enumerate(self.alignments):
      if a.is_unmapped: continue
      refidx,spos,epos = a.coordinates()
      feat = feature_lookup.lookup_interval(ref_lookup[refidx], spos, epos)
      if feat is not None:
        self.features[i] = feat
      else:
        if use_chrom:
          self.features[i] = '%s.%s' % (self.nofeature,ref_lookup[refidx])
        else:
          self.features[i] = self.nofeature

  def assign_best(self):
    import random
    from collections import defaultdict

    # Group alignments according to feature
    tmp = defaultdict(list)
    for f,a in zip(self.features,self.alignments):
      if f is not None: # Alignments that did not get assigned to a feature will have None
        tmp[f].append(a)

    # Assign the best alignment within each feature
    self.feat_aln_map = {}
    for f,alist in tmp.iteritems():
      if len(alist)==1:
        self.feat_aln_map[f] = alist[0]
      else:
        # Sort by edit distance, then by alignment score
        alist.sort(key=lambda x:x.NM)
        alist.sort(key=lambda x:x.AS,reverse=True)
        # Choose randomly from best alignments
        self.feat_aln_map[f] = random.choice([a for a in alist if a.AS==alist[0].AS and a.NM==alist[0].NM])

  def unique_feat(self):
    ''' Returns True if read maps to exactly one feature '''
    return len(self.feat_aln_map)==1

  def aligned_to_genome(self,genome_name):
    primary = self.feat_aln_map[genome_name]
    alternates = [a for f,a in zip(self.features,self.alignments) if f==genome_name]
    alternates = [a for a in alternates if not a == primary and a.AS==primary.AS and a.NM==primary.NM]
    return primary, alternates

  def structured_data(self):
    _genomes,alns = zip(*self.feat_aln_map.iteritems())
    _scores = [a.AS + a.query_length for a in alns]
    return [list(_genomes),_scores,[float(_scores[0])],max(_scores)]


