#! /usr/bin/env python

from copy import deepcopy
import math

class RSMatrix:
  def __init__(self, data, row_idx, col_idx):
    self.ridx = row_idx
    self.rownames = [k for k,v in sorted(self.ridx.iteritems(),key=lambda x:x[1])]   
    self.cidx = col_idx
    self.colnames = [k for k,v in sorted(self.cidx.iteritems(),key=lambda x:x[1])]   
    self.shape = (len(self.ridx),len(self.cidx))
    
    self.rows = [dict() for _ in self.ridx]
    for i,j,v in data:
      if j in self.rows[i]:
        print "WARNING: already have value"
      self.rows[i][j] = v

  def _func_row(self,func):
    ret = []
    for d in self.rows:
      vals = d.values()
      vals += [] if len(d) == self.shape[1] else [0]
      ret.append(func(vals))
    return ret
  
  def count(self):
    return self._func_row(lambda x:sum(v>0 for v in x))
    
  def max(self):
    return self._func_row(max)
  
  def min(self):
    return self._func_row(min)
  
  def sum(self):
    return self._func_row(sum)

  def copy(self):
    return deepcopy(self)

  """
  def mult_row(self, v, inplace=False):
    ''' Multiply each row by corresponding value in v '''
    other = self if inplace else deepcopy(self)
    assert isinstance(v, list), "v must be a list"
    assert len(v) == self.shape[0], "Shape is incorrect %d" % len(v)
    for i,r in enumerate(other.rows):
      for c in r.keys():
        r[c] *= v[i]
    return other
  
  def mult_item(self, v, inplace=False):
    ''' Multiply each item by v '''
    other = self if inplace else deepcopy(self)
    assert isinstance(v, (int, float)), "v must be int or float"
    for r in other.rows:
      for c in r.keys():
        r[c] *= v
    return other
  """
  def multiply(self, v, inplace=False):
    ''' Determine whether to multiply item or row based on v '''
    other = self if inplace else deepcopy(self)

    if type(v) is list:
      assert len(v) == self.shape[0], "Shape is incorrect %d" % len(v)
      for i,r in enumerate(other.rows):
        for c in r.keys():
          r[c] *= v[i]    
    else:
      for r in other.rows:
        for c in r.keys():
          r[c] *= v
    return other

  def add(self, v, inplace=False):
    ''' Determine whether to multiply item or row based on v '''
    other = self if inplace else deepcopy(self)

    if type(v) is list:
      assert len(v) == self.shape[0], "Shape is incorrect %d" % len(v)
      for i,r in enumerate(other.rows):
        for c in r.keys():
          r[c] += v[i]
    else:
      for r in other.rows:
        for c in r.keys():
          r[c] += v
    return other

  def exp(self, inplace=False):
    other = self if inplace else deepcopy(self)
    for r in other.rows:
      for c in r.keys():
        r[c] = math.exp(r[c])
    return other

  def expm1(self, inplace=False):
    other = self if inplace else deepcopy(self)
    for r in other.rows:
      for c in r.keys():
        r[c] = math.expm1(r[c])
    return other

  def normalize(self, inplace=False):
    ''' Normalize matrix by dividing elements in each row by the row sum'''
    other = self if inplace else deepcopy(self)
    _inv = [1./v if v != 0 else 1 for v in other.sum()]
    return other.multiply(_inv, inplace)

  def rescale(self, inplace=False):
    ''' Rescale matrix '''
    other = self if inplace else deepcopy(self)
    # Get max and min values of matrix
    _mmax = max(other.max())
    _mmin = min(other.min())
    if _mmin < 0:
      _scaling_factor = 100.0 / (_mmax - _mmin)
      other = other.add(abs(_mmin), inplace)
    else:
      _scaling_factor = 100.0 / _mmax
    return other.multiply(_scaling_factor).exp()

  def tocounts(self, inplace=False):
    '''
    :param inplace:
    :return:
    '''
    other = self if inplace else deepcopy(self)
    for r in other.rows:
      for c in r.keys():
        if r[c] != 0: r[c] = 1
    return other

  def col_sum(self):
    ret = [0] * self.shape[1]
    for r in self.rows:
      for j,v in r.iteritems():
        ret[j] += v
    return ret

  def col_count(self):
    ret = [0] * self.shape[1]
    for r in self.rows:
      for j,v in r.iteritems():
        ret[j] += 1
    return ret

  def __unicode__(self):
    ret = u'\t%s\n' % u'\t'.join(self.colnames)
    for r,row in enumerate(self.rows):
      ret += u'%s\t' % self.rownames[r]
      ret += u'%s\n' % '\t'.join(str(row[c]) if c in row else '0' for c in range(len(self.colnames)))
    return ret

  def __str__(self):
    return unicode(self).encode('utf-8')

#######################################################
# Utility functions
########################################################
def calculate_fractional_counts(mat):
  ''' Calculates the number of fractional reads for each reference
        Fractional count for read i in genome j:
          fc_ij = 0 if read i does not align to j
          fc_ij = 1. / (# of genomes read i aligns to)
  '''
  return mat.tocounts().normalize().col_sum()
  #ret = [0.] * mat.shape[1]
  #row_counts = mat.count()
  #for i,r in enumerate(mat.rows):
  #  for j,v in r.iteritems():
  #    ret[j] += 1. / row_counts[i]
  #return ret


def old_calculate_fractional_counts(mat):
  ret = [0.] * mat.shape[1]
  row_counts = mat.count()
  for i,r in enumerate(mat.rows):
    for j,v in r.iteritems():
      ret[j] += 1. / row_counts[i]
  return ret

def calculate_weighted_counts(mat):
  ''' Calculates the weighted number f reads for each reference
        Weighted count for read i in genome j:
          fc_ij = 0 if read i does not align to j
          fc_ij = weight_ij / sum(weight_iK for all genomes K)
  '''
  return mat.normalize().col_sum()

def old_calculate_weighted_counts(mat):
  ''' Calculates the number of weighted reads for each reference
        Fractional count for read i in genome j:
          fc_ij = 0 if read i does not align to j
          fc_ij = 1. / (# of genomes read i aligns to)
  '''
  ret = [0.] * mat.shape[1]
  row_sums = mat.sum()
  for i,r in enumerate(mat.rows):
    for j,v in r.iteritems():
      ret[j] += float(v) / row_sums[i]
  return ret

def calculate_unique_counts(mat):
  ''' Calculates the number of unique reads for each reference '''
  ret = [0] * mat.shape[1]
  row_counts = mat.count()
  for i,r in enumerate(mat.rows):
    if row_counts[i]==1:
      ret[r.keys()[0]] += 1
  return ret



