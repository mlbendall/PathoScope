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
  
  def __str__(self):
    ret = '\t%s\n' % '\t'.join(self.colnames)
    for r,row in enumerate(self.rows):
      ret += '%s\t' % self.rownames[r]
      ret += '%s\n' % '\t'.join(str(row[c]) if c in row else '0' for c in range(len(self.colnames)))
    return ret



