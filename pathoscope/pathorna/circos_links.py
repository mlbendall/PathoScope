#! /usr/bin/env python

__author__ = 'bendall'

import sys, os
import pysam

# Imports
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,pathoscopedir)
from pathoscope.pathorna.utils import iterread
from pathoscope.pathorna.utils import PSRead
from pathoscope.pathorna.utils import AnnotationLookup

"""
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
"""

'''
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
'''

def main(parser):
  args = parser.parse_args()
    
  samfile  = pysam.AlignmentFile(args.samfile)
  refnames = dict(enumerate(samfile.references))
  out1 = open('%s.links_feat_new.txt' % args.samfile,'w')
  #out2 = open('%s.links_nofeat.txt' % args.samfile,'w')  

  if args.gtffile is None:
    sys.exit("GTF file is required")

  flookup = AnnotationLookup(args.gtffile)

  #outfile = pysam.AlignmentFile(args.outfile, 'wh', header=samfile.header)
  target = 'HML-2_3q12.3'
  nbest_target = 0
  for rname,segments in iterread(samfile):
    r = PSRead(rname,segments)
    if r.is_unmapped:
      continue
    r.assign_feats(refnames, flookup)
    r.assign_best()
    if target in r.feat_aln_map:                    # r has alignment to target
      nbest_target += 1            
      if len(r.feat_aln_map) > 1:                   # r aligns to multiple features
        if True: #r.feat_aln_map[target].AS == r.bestAS:   # AS to target is best overall
          fname1 = target
          aln1 = r.feat_aln_map[target]
          rid1, s1, e1 = aln1.coordinates()
          ref1 = 'hs%s' % refnames[rid1][3:]
          for fname2,aln2 in zip(r.features,r.alignments):
            if aln1 == aln2: continue
            rid2, s2, e2 = aln2.coordinates()
            ref2 = 'hs%s' % refnames[rid2][3:]
            if fname2!='__nofeature__':
              addtext = 'color=lred_a5'
            else:
              addtext = 'color=lred_a5'
            print >>out1, '%s\t%d\t%d\t%s\t%d\t%d\t%s' % (ref1, s1, e1, ref2, s2, e2, fname2)
            #if fname2!='__nofeature__':
            #  print >>out1, '%s\t%d\t%d\t%s\t%d\t%d' % (ref1, s1, e1, ref2, s2, e2)
            #else:
            #  print >>out2, '%s\t%d\t%d\t%s\t%d\t%d' % (ref1, s1, e1, ref2, s2, e2)

  samfile.close()
  out1.close()
  #out2.close()
  print >>sys.stderr, 'reads: %d' % nbest_target

"""


      
    
    for a in psread.alignments:
      if a.AS == psread.bestAS:
        _ = a.set_tags('YC', c2str(PALETTE['bluish_green']))
    if r.aligns_to_feat():
      if len(r.feat_aln_map) == 77:
        counter += 1
        mycolor = '%s-%s' % (colprefix,counter)
        for fname1,aln1 in r.feat_aln_map.iteritems():
          rid1, s1, e1 = aln1.coordinates()
          ref1 = 'hs%s' % refnames[rid1][3:]
          for fname2,aln2 in r.feat_aln_map.iteritems():
            if fname1==fname2: continue
            rid2, s2, e2 = aln2.coordinates()
            ref2 = 'hs%s' % refnames[rid2][3:]
            print '%s\t%d\t%d\t%s\t%d\t%d\tcolor=%s' % (ref1, s1, e1, ref2, s2, e2, mycolor)

        
      maxfeat = max(maxfeat,len(r.feat_aln_map))
      continue
      if len(r.feat_aln_map) > 31 and filterby in r.feat_aln_map:
        maxfeat = max(maxfeat,len(r.feat_aln_map))
        counter += 1
        mycolor = '%s-%s' % (colprefix,counter)
        fname1 = filterby
        aln1 = r.feat_aln_map[fname1]
        #for fname1,aln1 in r.feat_aln_map.iteritems():
        rid1, s1, e1 = aln1.coordinates()
        ref1 = 'hs%s' % refnames[rid1][3:]
        for fname2,aln2 in r.feat_aln_map.iteritems():
            if fname1==fname2: continue
            rid2, s2, e2 = aln2.coordinates()
            ref2 = 'hs%s' % refnames[rid2][3:]
            print '%s\t%d\t%d\t%s\t%d\t%d\tcolor=%s' % (ref1, s1, e1, ref2, s2, e2, mycolor)
        
        #if counter >= 10:
        #  print >>sys.stderr, "exceeded %d" % counter
        #  break

"""

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
  # parser.add_argument('outfile', nargs="?", default="-", help='Output file (default is STDOUT)')
  main(parser)



