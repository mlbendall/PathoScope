#! /usr/bin/env python

__author__ = 'bendall'

#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Gets alignment file (currently support sam or BLAST-m8 format (.bl8)) and runs EM algorithm.
# Outputs the pathogens rank in the sample as a report that can be opened in Excel.
# Optionally outputs an updated alignment file (sam/bl8)

#usage information: pathoscope.py -h

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

def main(parser):
  # Imports
  from time import time
  pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
  sys.path.insert(0,pathoscopedir)
  from pathoscope.pathorna import PathoRNA

  # Process arguments
  args = parser.parse_args()
  adict = vars(args)
  annotation_file = adict.pop('gtffile')
  alignment_file  = adict.pop('samfile')
  assert os.path.exists(annotation_file) and os.path.isfile(annotation_file), 'ERROR: invalid annotation file [%s]' % annotation_file
  assert os.path.exists(alignment_file) and os.path.isfile(alignment_file), 'ERROR: invalid alignment file [%s]' % alignment_file
  opts = PathoRNA.PathoRNAOptions(alignment_file, annotation_file, **adict)
  if opts.verbose:
    print >>sys.stderr, opts

  start = time()
  PathoRNA.pathoscope_rna_reassign(opts)
  elapsed = time() - start
  if opts.verbose:
    print >>sys.stderr, "EM Elapsed Time: %d" % (elapsed)

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description="Pathoscope RNA")
  parser.add_argument('--verbose', action='store_true', help='Prints verbose text while running')
  parser.add_argument('--score_cutoff', type=float, default=0.01, help='Minimum final probability score for alignment')
  parser.add_argument('--out_matrix', action='store_true', dest='out_matrix', help='Output alignment matrix')
  parser.add_argument('--no_updated_sam', action='store_true', dest='no_updated_sam', help='Do not generate an updated alignment file')
  parser.add_argument('--emEpsilon', type=float, default=1e-7, help='EM Algorithm Epsilon cutoff')
  parser.add_argument('--maxIter', type=int, default=50, help='EM Algorithm maximum iterations')
  parser.add_argument('--piPrior', type=int, default=0, help='prior on pi')
  parser.add_argument('--thetaPrior', type=int, default=0, help='prior on theta')
  parser.add_argument('--exp_tag', default="pathorna", help='Experiment tag')
  parser.add_argument('--outdir', default=".", help='Output Directory')
  parser.add_argument('--version', action='version', version='%(prog)s 1.0')
  parser.add_argument('--ali_format', default='sam', help='Alignment Format: sam')
  parser.add_argument('--out_samfile', help='Name for updated alignment file')
  parser.add_argument('samfile', help='Alignment file path')
  parser.add_argument('gtffile', help='Annotation file (GTF format)')
  main(parser)
