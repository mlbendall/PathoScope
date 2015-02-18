__author__ = 'bendall'

import pysam

NO_FEATURE_KEY = '__nofeature__'

class FeatureLookup:
  def __init__(self,gtffile,attr_name="locus"):
    from collections import defaultdict
    import re
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

  def feature_name(self,id):
    return self._locus[id]



def iterread(samfile):
  ''' Each iteration returns all alignments that have the same read ID
      The file is expected to be sorted by queryname (default output from bowtie2)
      Parameters:
        samfile - pysam.AlignmentFile
  '''
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

class PSAlignment:
  """
  One distinct alignment. For properly paired segments (paired end) this has two segments,
  otherwise has one segment
  """
  def __init__(self, s1, s2=None):
    self.seg1 = s1
    self.seg2 = s2
    if s2 is None:
      assert not s1.is_proper_pair, "ERROR: s1 was properly paired but s2 is None"
      # Get flag values from s1
      self.is_paired      = s1.is_paired
      self.is_proper_pair = s1.is_proper_pair
      self.is_unmapped    = s1.is_unmapped
      self.secondary      = s1.is_secondary
      self.query_length   = s1.query_length
      tags = [dict(s1.tags)]
      # self.AS = tags1['AS'] if 'AS' in tags1 else 0
      # self.NM = tags1['NM'] if 'NM' in tags1 else 0
    else:
      self.is_paired      = s1.is_paired and s2.is_paired
      self.is_proper_pair = s1.is_proper_pair and s2.is_proper_pair
      self.is_unmapped    = s1.is_unmapped or s2.is_unmapped
      self.is_secondary   = s1.is_secondary or s2.is_secondary
      self.query_length   = s1.query_length + s2.query_length
      tags = [dict(s1.tags), dict(s2.tags)]

    self.AS = sum(t['AS'] if 'AS' in t else 0 for t in tags)
    self.NM = sum(t['NM'] if 'NM' in t else 0 for t in tags)

  def set_tags(self,tag,value):
    self.seg1.setTag(tag,value)
    if self.seg2 is not None:
      self.seg2.setTag(tag,value)

  def coordinates(self):
    if self.is_unmapped:
      return None
    tmp = [self.seg1.reference_start, self.seg1.reference_end]
    ref = self.seg1.reference_id
    if self.seg2 is not None:
      assert ref == self.seg2.reference_id
      tmp += [self.seg2.reference_start, self.seg2.reference_end]
    return ref,min(tmp),max(tmp)

  def write_samfile(self,outfile):
    retval = outfile.write(self.seg1)
    if self.seg2 is not None:
      retval = outfile.write(self.seg2)

  def __str__(self):
    return '%s\n%s\n%s' % (self.AS,self.seg1,self.seg2)

class PSRead:
  """
  All alignments for one read
  """
  nofeature = NO_FEATURE_KEY
  def __init__(self,readname,segs):
    self.readname = readname
    if segs[0].is_proper_pair:
      self.alignments = [PSAlignment(segs[i],segs[i+1]) for i in range(0,len(segs),2)]
    else:
      self.alignments = [PSAlignment(s) for s in segs]

    # Read is unmapped if all alignments are unmapped
    self.is_unmapped    = all(a.is_unmapped for a in self.alignments)
    # Read is uniquely if the number of mapped alignments == 1
    self.is_unique      = sum(not a.is_unmapped for a in self.alignments) == 1
    # Read is properly paired if all alignments are properly paired
    self.is_proper_pair = all(a.is_proper_pair for a in self.alignments)
    # features
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

  def best_for_feature(self,feat):
    import random
    alist = [a for f,a in zip(self.features,self.alignments) if f==feat]
    if len(alist)==1:
      return alist[0]
    else:
      alist.sort(key=lambda x:x.NM)
      alist.sort(key=lambda x:x.AS,reverse=True)
      poss = [a for a in alist if a.AS==alist[0].AS and a.NM==alist[0].NM]
      if len(poss)==1:
        return poss[0]
      else:
        return random.choice(poss)

  def best_alignments(self):
    from collections import defaultdict
    from random import choice

    bpf = []
    # group alignments according to feature
    byfeat = defaultdict(list)
    for f,a in zip(self.features,self.alignments):
      if f is not None:
        byfeat[f].append(a)

    # get the best alignment for each feature
    for f,alist in byfeat.iteritems():
      if len(alist)==1:            # Only one alignment to this feature
        bpf.append((f,alist[0]))
      else:                        # Multiple alignments to this feature
        alist.sort(key=lambda x:x.NM)
        alist.sort(key=lambda x:x.AS,reverse=True)
        poss = [a for a in alist if a.AS==alist[0].AS and a.NM==alist[0].NM]
        if len(poss) == 1:
          bpf.append((f,poss[0]))
        else:                      # Choose randomly from best alignments
          bpf.append((f,choice(poss)))
    return bpf

  def structured_data(self):
    baln = self.best_alignments()
    baln.sort(key=lambda x:x[1].AS+x[1].query_length, reverse=True)
    _genomes = [_[0] for _ in baln]
    _scores  = [(_[1].AS + _[1].query_length) for _ in baln]
    _norm    = [float(_scores[0])]
    _best    = _scores[0]
    return self.readname,[_genomes,_scores,_norm,_best]
