#!/usr/bin/env python
'''
  find barcodes matching a position or region
'''

import argparse
import collections
import gzip
import logging
import sys

import pysam

def main(umi, bam, positions):
  logging.info('starting...')
  reads_list = {}
  samfile = pysam.AlignmentFile(bam, "rb" )
  for position in positions:
    logging.info('analysing %s', position)
    chrom, span = position.split(':')
    start, finish = [int(x) for x in span.split('-')]
    reads = set()
    for column in samfile.pileup(chrom, start, finish):
      for read in column.pileups:
        reads.add(read.alignment.query_name)
    logging.info('%i reads overlapping %s', len(reads), position)
    reads_list[position] = reads

  
  # now find barcodes 
  logging.info('reading %s...', umi)
  barcodes = {}
  matches = {}
  for idx, line in enumerate(gzip.open(umi, 'r')):
    if idx % 4 == 0:
      name = line.strip().decode('UTF-8')[1:].split(' ')[0]
      for position in positions:
        if name in reads_list[position]:
          matches[position] = True
        else:
          matches[position] = False
    elif idx % 4 == 1:
      for position in positions:
        if matches[position]:
          #logging.debug('adding %s', line.strip().decode('UTF-8'))
          barcode = line.strip().decode('UTF-8')
          if position not in barcodes:
            barcodes[position] = collections.defaultdict(int)
          barcodes[position][barcode] += 1
    if idx % 1000000 == 0:
      logging.debug('%i umi lines with %s barcodes found', idx, ','.join(['{}'.format(len(barcodes[x])) for x in barcodes]))
 
  for position in positions:
    sys.stdout.write('Position: {}\n'.format(position))
    hist = collections.defaultdict(int)
    for barcode in barcodes[position]:
      hist[barcodes[position][barcode]] += 1
    sys.stdout.write('BarcodeCount\tTimesSeen\n')
    for h in hist:
      sys.stdout.write('{}\t{}\n'.format(h, hist[h]))
    logging.info('matched %i distinct barcodes for %i reads', len(barcodes[position]), sum([barcodes[position][x] for x in barcodes]))

  logging.info('done.')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Finds reads matching a barcode')
  parser.add_argument('--umi', required=True, help='umi file')
  parser.add_argument('--bam', required=True, help='alignment file')
  parser.add_argument('--positions', required=True, nargs='+', help='positions of interest')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.umi, args.bam, args.positions)
