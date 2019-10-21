#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import collections
import gzip
import logging
import sys

def main(umi):
  logging.info('starting...')
  counts = collections.defaultdict(int)
  total = 0
  for idx, line in enumerate(gzip.open(umi, 'r'), gzip.open(r1), gzip.open(r2)):
    if idx % 4 == 1:
      index = line.strip().decode('UTF-8')
      counts[index] += 1
      total += 1

    if idx % 100000 == 0:
      logging.debug('processed %i: %i found...', idx, len(counts)) 

  sys.stdout.write('Barcode\tCount\t%\n')
  for k in sorted(counts.keys()):
    sys.stdout.write('{}\t{}\t{:.2f}\n'.format(str(k), counts[k], counts[k] / total * 100))

  logging.info('done. %i tags', len(counts))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Count UMI tags')
  parser.add_argument('--umi', required=True, help='umi file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.umi)
