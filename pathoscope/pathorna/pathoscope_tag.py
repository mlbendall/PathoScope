#! /usr/bin/env python

__author__ = 'bendall'

import sys, os
import pysam

# Imports
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir)
from pathoscope.pathorna.utils import iterread
from pathoscope.pathorna.utils import PSRead\
from pathoscope.pathorna.utils import AnnotationLookup

PALETTE = {'orange':             (230, 159, 0),
           'sky_blue':           (86, 180, 233),
           'bluish_green':       (0, 158, 115), # 009E73
           'yellow':             (240, 228, 66),
           'blue':               (0, 114, 178),
           'vermillion':         (213, 94, 0),
           'reddish_purple':     (204, 121, 167),
           'plum':               (142, 69, 133),
          }

blues = [(235,247,244),
         (215,240,233),
         (196,232,222),
         (176,225,211),
        ]

def lighten_color(rgb,percent):
  rF = max(min(255,int(rgb[0] * percent)),0)
  gF = max(min(255,int(rgb[1] * percent)),0)
  bF = max(min(255,int(rgb[2] * percent)),0)
  return rF,gF,bF

def c2str(rgb):
  return '%d,%d,%d' % rgb

def tag_color(psread):
  """ Set the YC tag (Viewers use this for the color)
  :param psread:
  :return:
  """
  if psread.is_unmapped:
    return
  if psread.is_unique:
    _ = psread.alignments[0].set_tags('YC', c2str(PALETTE['vermillion']))
  else:
    # bestAS   = max(a.AS for a in psread.alignments)
    # num_best =
    for a in psread.alignments:
      if a.AS == psread.bestAS:
        _ = a.set_tags('YC', c2str(PALETTE['bluish_green']))
      else:
        pct = float(a.AS) / psread.bestAS
        cat = int(pct*4)
        _ = a.set_tags('YC', c2str(blues[cat]))

def tag_nonunique(psread,has_features=True):
  if psread.is_unmapped or psread.is_unique:
    return
  else:
    # bestAS    = max(a.AS for a in psread.alignments)
    num_best  = sum(a.AS==psread.bestAS for a in psread.alignments)
    if has_features:
      bestfeats = set([f for a,f in zip(psread.alignments, psread.features) if a.AS == psread.bestAS])
      bestfeat_str = ','.join(sorted(bestfeats))
    for a in psread.alignments:
      _ = a.set_tags('ZN',len(psread.alignments)).set_tags('ZB',num_best).set_tags('ZS',psread.bestAS)
      if has_features:
        _ = a.set_tags('ZF',bestfeat_str)

def main(parser):
  args = parser.parse_args()

  samfile  = pysam.AlignmentFile(args.samfile)
  refnames = dict(enumerate(samfile.references))

  has_features = args.gtffile is not None
  if has_features:
    flookup = AnnotationLookup(args.gtffile)

  outfile = pysam.AlignmentFile(args.outfile, 'wh', header=samfile.header)

  for rname,segments in iterread(samfile):
    r = PSRead(rname,segments)
    if has_features: r.assign_feats(refnames, flookup)
    tag_color(r)
    tag_nonunique(r,has_features)
    for a in r.alignments:
      a.write_samfile(outfile)

  samfile.close()
  outfile.close()

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(prog="pathoscope_tag.py",
                                   description="Pathoscope TAG. Adds tags to SAM file",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  )
  parser.add_argument('-v','--version', action='version', version='%(prog)s 0.9b')
  parser.add_argument('--ali_format', default='sam', help='Alignment Format. Only sam is supported.')
  parser.add_argument('--gtffile', help='Path to annotation file (GTF format)')

  parser.add_argument('samfile', nargs="?", default="-", help='Path to alignment file (default is STDIN)')
  parser.add_argument('outfile', nargs="?", default="-", help='Output file (default is STDOUT)')
  main(parser)



