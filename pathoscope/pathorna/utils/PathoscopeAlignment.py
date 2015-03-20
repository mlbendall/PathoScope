__author__ = 'bendall'

class PathoscopeAlignment:
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

