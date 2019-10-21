#!/usr/bin/env python
'''
  find reads matching a barcode
'''

import argparse
import collections
import gzip
import logging
import sys

def main(umi, r1, r2, barcode):
  logging.info('starting...')
  total = 0
  sys.stdout.write('R1\tR2\n')
  for idx, (line, r1_line, r2_line) in enumerate(zip(gzip.open(umi, 'r'), gzip.open(r1), gzip.open(r2))):
    if idx % 4 == 1:
      index = line.strip().decode('UTF-8')
      total += 1
      if index == highlight:
        sys.stdout.write('{}\t{}\n'.format(r1_line.strip().decode('UTF-8'), r2_line.strip().decode('UTF-8')))

    if idx % 100000 == 0:
      logging.debug('processed %i...', idx) 

  logging.info('done. %i reads', len(total))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Finds reads matching a barcode')
  parser.add_argument('--umi', required=True, help='umi file')
  parser.add_argument('--r1', required=False, help='fq r1 file')
  parser.add_argument('--r2', required=False, help='fq r2 file')
  parser.add_argument('--barcode', required=False, help='barcode to find')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.umi, args.r1, args.r2, args.barcode)
