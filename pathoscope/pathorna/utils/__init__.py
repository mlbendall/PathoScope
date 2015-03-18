__author__ = 'bendall'

import pysam
import math
from collections import defaultdict
import random

from AnnotationLookup import AnnotationLookup

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
  #a.tags                 = template.tags
  a.is_secondary         = mate.is_secondary
  a.is_paired      = True
  a.is_proper_pair = False
  a.is_unmapped    = True
  # This tag indicates the segment is a "mock pair"
  a.setTag('YT', mate.get_tag('YT'))
  if add_tag: a.setTag('ZT',"MP")
  return a

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
      self.is_paired       = True
      self.is_proper_pair = s1.is_proper_pair and s2.is_proper_pair
      self.is_unmapped    = s1.is_unmapped and s2.is_unmapped
      self.is_secondary   = s1.is_secondary and s2.is_secondary
      self.query_length   = sum(s.query_length if not s.is_unmapped else 0 for s in [s1,s2])
      assert s1.is_secondary == s2.is_secondary

    tags = [dict(s.tags) for s in [s1,s2] if s is not None and not s.is_unmapped]
    self.AS = sum(t['AS'] if 'AS' in t else 0 for t in tags)
    self.NM = sum(t['NM'] if 'AS' in t else 0 for t in tags)

  def set_tags(self, tag, value=None):
    ''' Set tags for alignment
    :param tag:
    :param value:
    :return:
    '''
    if value is None:
      assert type(tag) is dict, "ERROR: set_tags requires dict if single arg is provided"
      for k,v in tag.iteritems():
        self.seg1.setTag(k,v)
        if self.is_paired:
          self.seg2.setTag(k,v)
    else:
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

  def set_secondary(self,value):
    # Ensure MAPQ in acceptable range
    self.is_secondary      = value
    self.seg1.is_secondary = value
    if self.is_paired:
      self.seg1.is_secondary = value

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

  def __init__(self, readname, segs):
    self.readname       = readname
    self.is_paired      = None
    self.is_proper_pair = None
    self.is_unmapped    = None    # Read is unmapped if all alignments are unmapped
    self.is_unique      = None    # Read is unique if only one mapped alignment

    self.alignments     = []
    self.features       = []
    self.feat_aln_map   = {}

    self.bestAS         = float('-inf')

    self._make_alignments(segs)

  def _make_alignments(self, segs):
    """ Logic for making PSAlignments from segments
    :param segs:
    :return:
    """
    self.is_paired = segs[0].is_paired
    assert all(s.is_paired == self.is_paired for s in segs), "Expected all segments to be paired or unpaired!\n%s" % [str(s) for s in segs]
    if self.is_paired:
      self._make_paired(segs)
    else:
      self._make_unpaired(segs)

    # Set the primary alignment for this read
    self._set_primary()
    # Read is unique if the number of mapped alignments == 1
    self.is_unique       = len([a.is_unmapped==False for a in self.alignments]) == 1
    # Get the best alignments
    self.bestAS          = self.alignments[0].AS #max(a.AS for a in self.alignments)

  def _make_unpaired(self,segs):
    self.is_proper_pair = False
    self.is_unmapped    = segs[0].is_unmapped
    for s in segs:
      self.alignments.append(PSAlignment(s))

  def _make_paired(self,segs):
    self.is_proper_pair = segs[0].is_proper_pair
    assert all(s.is_proper_pair==self.is_proper_pair for s in segs), "Expected all segments to be proper or not!\n%s" % [str(s) for s in segs]
    # If properly paired, assume that mate is on the adjacent line
    if self.is_proper_pair:
      self.is_unmapped    = False
      for i in range(0,len(segs),2):
        self.alignments.append(PSAlignment(segs[i],segs[i+1]))
        #self.alignments     = [PSAlignment(segs[i],segs[i+1]) for i in range(0,len(segs),2)]
    else:
      ''' Bowtie handling of segements that are not properly paired
          There will be two primary segments (one for read1, one for read2)
          These might be unmapped
      '''
      primary   = []
      secondary = []
      for seg in segs:
        if not seg.is_secondary:
          primary.append(seg)
        else:
          secondary.append(seg)

      assert len(primary) == 2, "Expected exactly two primary segments!\n%s" % [str(s) for s in segs]

      # If either read1 or read2 has no mappings, its primary segment will be unmapped
      # Create a primary PSAlignment containing both segments
      if any(s.is_unmapped for s in primary):
          self.is_unmapped = all(s.is_unmapped for s in primary) # Unmapped if both segments are unmapped
          self.alignments.append(PSAlignment(*primary))

      # Both segments have at least one alignment (could be discordant or neither)
      # Create a separate PSAlignment for each segment, make the better one primary
      else:
        self.is_unmapped = False
        # Create mate for primary[0] using primary[1] as template
        a1 = PSAlignment(primary[0], make_unmapped_mate(primary[0], template=primary[1]))
        # Create mate for primary[1] using primary[0] as template
        a2 = PSAlignment(primary[1], make_unmapped_mate(primary[1], template=primary[0]))
        self.alignments.append(a1)
        self.alignments.append(a2)

      # Make secondary PSAlignments for remaining alignments
      if secondary:
        # Fast lookup of the primary template
        lookup_r1 = dict((True if _.is_read1 else False, _) for _ in primary)
        assert lookup_r1[True].is_read1
        assert lookup_r1[False].is_read2
        for s in secondary:
          assert s.is_secondary
          m = make_unmapped_mate(s,template=lookup_r1[s.is_read2])
          self.alignments.append(PSAlignment(s,m))

  def _set_primary(self):
    """ Logic for setting the primary alignment for PSRead
    :return:
    """
    self.alignments.sort(key=lambda x:x.AS, reverse=True)
    self.alignments[0].set_secondary(False)
    for a in self.alignments[1:]:
      a.set_secondary(True)

  def assign_feats(self, ref_lookup, feature_lookup, use_chrom=False):
    """ Assign each alignment to a feature
    :param ref_lookup:      Dictionary mapping reference index to reference name
    :param feature_lookup:  FeatureLookup object
    :param use_chrom:       No feature reads
    :return:
    """
    assert not self.features, "Features have already been assigned!"
    for i,a in enumerate(self.alignments):
      if a.is_unmapped:
        self.features.append(None)
        continue
      refidx,spos,epos = a.coordinates()
      feat = feature_lookup.lookup_interval(ref_lookup[refidx], spos, epos)
      if feat is not None:
        self.features.append(feat)
      else:
        if use_chrom:
          self.features.append('%s.%s' % (self.nofeature,ref_lookup[refidx]))
        else:
          self.features.append(self.nofeature)

  def assign_best(self):
    """ One alignment is selected for features with multiple alignments
    :return:
    """
    # One or more features has multiple alignments
    # Group alignments according to feature
    tmp = defaultdict(list)
    for f,a in zip(self.features,self.alignments):
      if f is not None: # Alignments that did not get assigned to a feature will have None
        tmp[f].append(a)

    while tmp:
      f,alist = tmp.popitem()
      if len(alist) == 1:
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

  def aligns_to_feat(self):
    ''' Returns True if read aligns to an annotated feature '''
    return any(f != self.nofeature for f in self.features)

  def aligned_to_genome(self,genome_name):
    ''' Returns all alignments to the given genome name
    :param genome_name: Name of genome
    :return:
    '''
    primary = self.feat_aln_map[genome_name]
    alternates = [a for f,a in zip(self.features,self.alignments) if f==genome_name]
    alternates = [a for a in alternates if not a == primary and a.AS==primary.AS and a.NM==primary.NM]
    return primary, alternates

  def structured_data(self):
    '''
    Creates data structure for read used by pathoscope_em and other functions
    :return: [list(genomes), list(scores), list(score), max score]
    '''
    _genomes,alns = zip(*self.feat_aln_map.iteritems())
    _scores = [a.AS + a.query_length for a in alns]
    return [list(_genomes),_scores,[float(_scores[0])],max(_scores)]


