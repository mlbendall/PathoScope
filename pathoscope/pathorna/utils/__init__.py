__author__ = 'bendall'

import math
from AnnotationLookup import AnnotationLookup
from PathoscopeRead import PathoscopeRead as PSRead
from RSMatrix import *

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

def array_to_prop(a):
  asum = sum(a)
  return [float(v) / asum for v in a]