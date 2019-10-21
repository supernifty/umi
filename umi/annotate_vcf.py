#!/usr/bin/env python
'''
  find barcodes matching a position or region
'''

import argparse
import collections
import gzip
import logging
import sys

import cyvcf2
import pysam

def main(umi, bam):
  # find barcodes 
  logging.info('reading %s...', umi)
  barcodes = {} # readname -> barcode
  seen = set()
  idx = 0
  for idx, line in enumerate(gzip.open(umi, 'r')):
    if idx % 4 == 0:
      name = line.strip().decode('UTF-8')[1:].split(' ')[0] # read name
    elif idx % 4 == 1:
      barcode = line.strip().decode('UTF-8')
      barcodes[name] = barcode
      seen.add(barcode)
    if idx % 1000000 == 0:
      logging.debug('%i umi lines processed with %i distinct barcodes...', idx, len(seen))
  logging.info('%i umi lines processed with %i distinct barcodes...', idx, len(seen))

  logging.info('reading %s and stdin (vcf)...', bam)
  samfile = pysam.AlignmentFile(bam, 'rb')
  vcf_in = cyvcf2.VCF('-')
  vcf_in.add_info_to_header({'ID': 'umi', 'Description': 'Number of distinct UMI barcodes and number of reads overlapping the variant', 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)

  idx = 0
  worst = (1, None)
  best = (0, None)
  for idx, variant in enumerate(vcf_in):
    position = '{}:{}'.format(variant.CHROM, variant.POS)
    logging.debug('analysing %s', position)
    overlapping_barcodes = set()
    reads = set()
    for column in samfile.pileup(variant.CHROM, variant.POS, variant.POS + 1):
      for read in column.pileups:
        overlapping_barcodes.add(barcodes[read.alignment.query_name])
        reads.add(read.alignment.query_name)
    logging.debug('%i reads overlapping %s with %i distinct barcodes', len(reads), position, len(overlapping_barcodes))
    variant.INFO['umi'] = '{}/{}'.format(len(overlapping_barcodes), len(reads))

    proportion = len(overlapping_barcodes) / len(reads)
    if proportion < worst[0]:
      worst = (proportion, position)
    if proportion > best[0]:
      best = (proportion, position)

    sys.stdout.write(str(variant))
    
    if idx < 5 or idx % 1000 == 0:
      logging.debug('read %i variants...', idx)

  logging.info('done writing %i variants. best proportion %s worst proportion %s', idx, best, worst)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Finds reads matching a barcode')
  parser.add_argument('--umi', required=True, help='umi file')
  parser.add_argument('--bam', required=True, help='alignment file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.umi, args.bam)

