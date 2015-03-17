#! /usr/bin/env python

__author__ = 'bendall'

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
  assert os.path.exists(annotation_file) and os.path.isfile(annotation_file), 'ERROR: annotation file [%s] does not exist' % annotation_file
  assert os.path.exists(alignment_file) and os.path.isfile(alignment_file), 'ERROR: alignment file [%s] does not exist' % alignment_file
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
  parser = argparse.ArgumentParser(prog="pathoscope_rna.py",
                                   description="Pathoscope RNA",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  )
  parser.add_argument('-v','--version', action='version', version='%(prog)s 0.9b')

  inputopts = parser.add_argument_group('input', 'Input options')
  inputopts.add_argument('--ali_format', default='sam', help='Alignment Format. Only sam is supported.')
  inputopts.add_argument('samfile', help='Path to alignment file')
  inputopts.add_argument('gtffile', help='Path to annotation file (GTF format)')
  inputopts.add_argument('--no_feature_key',
                         help='Feature name for unassigned reads. Must not match any other feature name')

  outputopts = parser.add_argument_group('output', 'Output options')
  outputopts.add_argument('--verbose', action='store_true',
                          help='Prints verbose text while running')
  outputopts.add_argument('--outdir', default=".",
                          help='Output Directory')
  outputopts.add_argument('--out_samfile',
                          help='Name for updated alignment file')
  outputopts.add_argument('--exp_tag', default="pathorna",
                          help='Experiment tag')
  outputopts.add_argument('--report_all', action='store_true',
                          help='Include all genomes in report')
  outputopts.add_argument('--min_final_guess', type=float, default=0.01,
                          help='Minimum final guess for genome to appear in report. Genomes with one or more final hits will always be included.')
  outputopts.add_argument('--out_matrix', action='store_true',
                          help='Output alignment matrix')
  outputopts.add_argument('--out_abundance', action='store_true',
                          help='Output abundances (FPKM, TPI)')
  outputopts.add_argument('--no_updated_sam', action='store_true', dest='no_updated_sam',
                          help='Do not generate an updated alignment file')

  modelopts = parser.add_argument_group('model', 'Model parameters')
  modelopts.add_argument('--piPrior', type=int, default=0,
                         help='Pi Prior equivalent to adding n unique reads')
  modelopts.add_argument('--thetaPrior', type=int, default=0,
                         help='Theta Prior equivalent to adding n non-unique reads')
  modelopts.add_argument('--score_cutoff', type=float, default=0.01,
                         help='Minimum final probability score for alignment')

  emopts = parser.add_argument_group('em', 'EM parameters')
  emopts.add_argument('--emEpsilon', type=float, default=1e-7,
                      help='EM Algorithm Epsilon cutoff')
  emopts.add_argument('--maxIter', type=int, default=100,
                      help='EM Algorithm maximum iterations')

  #  parser.add_argument('--out_samfile', help='Name for updated alignment file')
  main(parser)
