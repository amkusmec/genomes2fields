import argparse
import textwrap
import sys

###############################################################################
def warning(*objs):
  print("WARNING: ", *objs, end = '\n', file = sys.stderr)
  sys.exit()

################################################################################
def version():
  v0 = """
  ##############################################################################
  recode_indels
  (c) 2018 Aaron Kusmec
  
  Replace +/- insertion/deletion in VCF files with A/T coding.
  
  ##############################################################################
  """
  
  return v0

################################################################################
def get_parser():
  parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = textwrap.dedent(version()))
  
  parser.add_argument('-i', '--input', help = 'Input VCF file', type = str)
  parser.add_argument('-o', '--output', help = 'Output file stem', type = str)
  
  return parser

################################################################################
if __name__ == "__main__":
  parser = get_parser()
  args = vars(parser.parse_args())
  
  if args['input'] is None:
    warning("No input file.")
  if args['output'] is None:
    warning("No output file stem.")
  
  print(version())
  
  swaps = []
  counter, replaced = 0, 0
  convert = {'-' : 'A', '+' : 'T'}
  
  with open(args['input'], 'r') as infile:
    with open(args['output'] + '.vcf', 'w') as outfile:
      for line in infile:
        if line[0] == "#":
          outfile.write(line)
          continue
        
        line = line.strip().split()
        if '+' in line[3] or '-' in line[4]:
          line[3] = convert['+']
          line[4] = convert['-']
          replaced += 1
          swaps.append(line[2])
        elif '-' in line[3] or '+' in line[4]:
          line[3] = convert['-']
          line[4] = convert['+']
          replaced += 1
          swaps.append(line[2])
        
        outfile.write('\t'.join(line) + '\n')
        counter += 1
        
        if counter % 25000 == 0:
          print("Processed [", str(counter), "] SNPs.", end = '\r')
  
  print()
  with open(args['output'] + '_replaced.txt', 'w') as outfile:
    outfile.write('\n'.join(swaps) + '\n')
  
  print("Finished.")
