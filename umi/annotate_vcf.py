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
  barcodes = {}
  for idx, line in enumerate(gzip.open(umi, 'r')):
    if idx % 4 == 0:
      name = line.strip().decode('UTF-8')[1:].split(' ')[0] # read name
    elif idx % 4 == 1:
      barcode = line.strip().decode('UTF-8')
      barcodes[name] = barcode
    if idx % 1000000 == 0:
      logging.info('%i umi lines processed...', idx)

  logging.info('reading %s and stdin (vcf)...', bam)
  samfile = pysam.AlignmentFile(bam, "rb" )
  vcf_in = cyvcf2.VCF('-')
  vcf_in.add_info_to_header({'ID': 'umi', 'Description': 'Number of distinct UMI barcodes and number of reads overlapping the variant', 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)

  for idx, variant in enumerate(vcf_in):
    logging.debug('analysing %s:%i', variant.CHROM, variant.POS)
    overlapping_barcodes = set()
    reads = 0
    for column in samfile.pileup(variant.CHROM, variant.POS, variant.POS + 1):
      for read in column.pileups:
        overlapping_barcodes.add(read.alignment.query_name)
        reads += 1
    logging.debug('%i reads overlapping %s:%i with %i distinct barcodes', reads, variant.CHROM, variant.POS, len(overlapping_barcodes))
    variant.INFO['umi'] = '{}/{}'.format(len(overlapping_barcodes), reads)
    sys.stdout.write(str(variant))
    
    if idx < 5 or idx % 1000 == 0:
      logging.info('read %i variants...', idx)

  logging.info('done.')

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

