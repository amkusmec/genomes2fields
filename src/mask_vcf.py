import argparse
import textwrap
import timeit
import os
from random import sample
from math import ceil

def warning(*objs):
  print("WARNING: ", *objs, end = '\n', file = sys.stderr)
  sys.exit()

def get_parser():
  parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)
  
  parser.add_argument('-p', '--path', help = 'Path of the input file', \
                      nargs = '?', default = os.getcwd())
  parser.add_argument('-i', '--input', help = 'Input file', type = str)
  parser.add_argument('-o', '--output', help = 'Output STEM', type = str)
  parser.add_argument('-m', '--mask', help = 'Proportion to mask', \
                      type = float, default = 0.1)
  
  return parser

if __name__ == '__main__':
  parser = get_parser()
  args = vars(parser.parse_args())
  
  if args['path'] is not None:
    os.chdir(args['path'])
  if args['input'] is None:
    warning("No input file.")
  if args['output'] is None:
    warning("No output stem.")
  
  st = timeit.default_timer()
  
  infile = open(args['input'], 'r')
  outfile = open(args['output'] + '.vcf', 'w')
  idxfile = open(args['output'] + '_index.txt', 'w')
  
  for line in infile:
    line = line.strip()
    if line[0] == "#":
      outfile.write(line + '\n')
      continue
    
    line = line.split()
    valid = [i for i in range(9, len(line)) if line[i][:3] != './.']
    if ceil((len(line) - 9)*args['mask']) > len(valid):
      warning("Not enough non-missing samples at " + line[2] + ".")
    idx = sample(valid, ceil((len(line) - 9)*args['mask']))
    idxfile.write(line[2] + ',' + ','.join([str(i) for i in idx]) + '\n')
    
    for i in idx: line[i] = './.' + line[i][3:]
    outfile.write('\t'.join(line) + '\n')
  
  infile.close()
  outfile.close()
  idxfile.close()
  
  et = timeit.default_timer()
  
  print("Masking finished.")
  print("Time: %.2f min." % ((et - st)/60))
