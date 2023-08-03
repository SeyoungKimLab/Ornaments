# SCRIPT TO TAKE ALL THE TRANSCRIPTOME SNPS, AND CONSOLIDATE INTO A FILE FOR A SUBSET OF THE SAMPLES

import argparse
import gzip
import sys

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input-dir', required=True, type=str)
    parser.add_argument('-o','--outfile', type=str)
    parser.add_argument('-s','--sample-file', required=True, type=str)
    return parser.parse_args()

# Arguments: first is location of directory, second is out file name, third is sample file or sample name
def main():
  args = parse_arguments()

  samples_file = args.sample_file
  samples = set([])
  try:
    with open(samples_file) as f:
      for line in f:
        if line.strip():
          samples.add(line.strip())
    samples = sorted(list(samples))
  except:
    samples = [samples_file]

  print('Extracting relevant SNPs for ornament transcriptome over: ', samples)
  i_dir = args.input_dir
  o_file = args.outfile
  with open(o_file, 'w+') as final_vcf:
    header_written = False
    for i in [*range(1, 23), 'X']:
      print('Processing chromosome ', i)
      #filename = 'chr' + str(i) + '.transcriptome.vcf'
      filename = 'chr' + str(i) + '.transcriptome.vcf'
      # filename = 'phase_chr' + str(i) + '.vcf.gz'
      num_processed = 0
      tot_lines = 0
      with open(i_dir + '/' + filename) as curr_chr_tvcf:
        start = False
        inds = []
        for line in curr_chr_tvcf:
          if 'CHROM' in line:
            start = True
            toks = [t.strip() for t in line.split('\t')]
            if len(samples) == 0:
              samples = toks[9:]
              samples = list(set(samples))
            inds = [toks[9:].index(sample) for sample in samples]
            list_of_samples = [toks[9:][ind] for ind in inds]
            if not header_written:
              final_vcf.write('\t'.join(toks[:9] + list_of_samples) + '\n')
              header_written = True
            continue

          if not start:
            continue

          toks = line.split('\t')
          gts = toks[9:]
          toks[7] = '.'
          current_snp_genotypes = [gts[ind].split(':')[0].strip() for ind in inds]
          write_record = False
          for gt in current_snp_genotypes:
            if gt == '1|0' or gt == '0|1' or gt == '1|1' or gt == '1':
              write_record = True
          if True:
            num_processed += 1
            final_vcf.write('\t'.join(toks[:9] +  current_snp_genotypes) + '\n')
          tot_lines += 1
        print('Processed ', num_processed, ' lines for samples ' + str(samples) + ' at indices', str(inds))

      

if __name__ == '__main__':
  main()
