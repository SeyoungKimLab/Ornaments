# SCRIPT TO CONVERT THE GENOMIC VCF FILE INTO A TRANSCRIPTOMIC VCF FILE FOR ALL THE SAMPLES
# CAN OPTIONALLY INTRODUCE FAKE SWITCH ERRORS INTO THE GENOME

import argparse
import numpy as np
import time
import gzip
import sys

from Bio.Seq import Seq
from gtfparse import read_gtf

########################################################################################################################
#######  Helper functions to perform the conversion of a genomic snp to a transcriptomic snp ###########
########################################################################################################################

# A genome position can be contained by many transcripts due to alternative splicing. Get all of them
# Get the transcripts that can be contained in this position
def get_transcript_containing_pos(chr, pos, chr_to_transcript, n=0):
  viable_transcripts = []
  if 'chr' in chr:
    chr_to_transcript_obj = chr_to_transcript[chr]
  elif (str(chr) == 'MT'):
    chr_to_transcript_obj = chr_to_transcript['chrM']
  else:
    chr_to_transcript_obj = chr_to_transcript['chr' + str(chr)]
  
  max_end_so_far = -1
  for i in range(max(0, n), len(chr_to_transcript_obj)):
    transcript = chr_to_transcript_obj[i]
    
    max_end_so_far = max(max_end_so_far, transcript['range'][1])
    if transcript['range'][1] < pos:
      if transcript['range'][1] >= max_end_so_far:
        n = i
    
    if transcript['range'][0] > pos:
      break
      
    # Transcript range contains the specified position.
    if pos >= transcript['range'][0] and pos <= transcript['range'][1]:
      # Now, check if the position is within an exon region
      for exon in transcript['exons']:
        if transcript['exons'][exon][0] <= pos and transcript['exons'][exon][1] >= pos:
          viable_transcripts.append(transcript)
          break
    
          
  return n, viable_transcripts

def complement_sequence(seq):
  return str(Seq(seq).reverse_complement())

def complement_nucleotide(nucleotide):
  if nucleotide == 'A':
    return 'T'
  elif nucleotide == 'T':
    return 'A'
  elif nucleotide == 'G':
    return 'C'
  elif nucleotide == 'C':
    return 'G'

# For positively stranded transcripts, we find the exon that it falls under, 
def convert_bp_positive(transcript, pos):
  exons = transcript['exons']
  # Iterating through exons already addresses the strandedness. In this function we assume this is positive
  curr_len = 0
  for i in exons:
    start, end = exons[i]
    # The SNP is contained in this exon
    if start <= pos and end >= pos:
      curr_len += (pos - start) + 1
      return curr_len
    else:
      curr_len += (end - start) + 1

def convert_bp_negative(transcript, pos):
  exons = transcript['exons']
  # Iterating through exons already addresses the strandedness, i.e. we are going in reverse order
  curr_len = 0
  tot_len = 0
  to_ret = 0
  for i in exons:
    start, end = exons[i]
    if start <= pos and end >= pos:
      curr_len += (pos - start) + 1
      to_ret = curr_len
    else:
      curr_len += (end - start) + 1
    tot_len += (end - start) + 1
  
  return tot_len - to_ret + 1

def transcript_generator(file):
  with open(file) as diploid_transcriptome:
    ref = ''
    for line in diploid_transcriptome:
      if '>' in line:
        # The current transcript has been finished is done
        if len(ref) > 0:
          yield (curr_tid, ref)

        curr_tid = line.split('>')[1].split(' ')[0].strip()
        ref = ''
      else:
        ref += line.strip()
        
    yield (curr_tid, ref)

def get_anchor_pos(transcript_obj, start_pos, ref, alt):
  prev = None
#     print(transcript_obj)
  for i in range(len(transcript_obj['exons'])):
    e = transcript_obj['exons'][i + 1]
    # starts before beginning of the transcript
    if start_pos < e[0] and prev is None:
      return None
    if prev:
      if start_pos < e[0] and start_pos > prev[1]:
        assert(start_pos + (len(ref) - len(alt)) >= e[0])
        tpos = convert_bp_positive(transcript_obj, prev[1])
        # indel can also run off end of transcript, in which case we need to truncate it
#                 print(start_pos + (len(ref) - len(alt)))
        tpos_end = convert_bp_positive(transcript_obj, start_pos + (len(ref) - len(alt)))
        # edge case: when both start and end are introns and cover entire exon(s)
#                 print(tpos, tpos_end)
        return(tpos, tpos_end)
    prev = e

def get_anchor_neg(transcript_obj, end_pos, ref, alt):
  prev = None
#     print(transcript_obj)
  for i in range(len(transcript_obj['exons'])):
    e = transcript_obj['exons'][i + 1]
#         print(i, e)
    # starts before beginning of the transcript
    if end_pos > e[1] and prev is None:
      return None
    if prev:
      if end_pos > e[1] and end_pos < prev[0]:
        assert(end_pos - (len(ref) - len(alt)) <= e[1])
        tpos = convert_bp_negative(transcript_obj, prev[0])
        # indel can also run off end of transcript, in which case we need to truncate it
#                 print(start_pos + (len(ref) - len(alt)))
        tpos_end = convert_bp_negative(transcript_obj, end_pos - (len(ref) - len(alt)))
        # edge case: when both start and end are introns and cover entire exon(s)
#                 print(tpos, tpos_end)
        return(tpos, tpos_end)
    prev = e

def truncate_pos(transcript_obj, start_pos, ref, alt):
  end_pos = start_pos + (len(ref) - len(alt))
  prev = None
#     print(transcript_obj)
  start_looking = False
  for i in range(len(transcript_obj['exons'])):
    e = transcript_obj['exons'][i + 1]
    # starts before beginning of the transcript
    if start_pos >= e[0] and start_pos <= e[1]:
      # want to find tpos of truncated end
      start_looking = True
    if start_looking and prev is not None:
      if end_pos > prev[1] and end_pos < e[0]:
        return convert_bp_positive(transcript_obj, prev[1])
    prev = e
  return convert_bp_positive(transcript_obj, prev[1])

def truncate_neg(transcript_obj, start_pos, ref, alt):
  end_pos = start_pos + (len(ref) - len(alt))
  prev = None
  start_looking = False
  for i in range(len(transcript_obj['exons'])):
    e = transcript_obj['exons'][i + 1]
#         print(i, e)
    if end_pos >= e[0] and end_pos <= e[1]:
      start_looking = True
    if start_looking and prev is not None:
      if start_pos < prev[0] and start_pos > e[1]:
        return convert_bp_negative(transcript_obj, prev[0])
    prev = e
  return convert_bp_negative(transcript_obj, prev[0])

def main(genome_vcf_file, gtf_file, transcriptome_file, output_file, coordinate_conversion_file):
  ########################################################################################################################
  #######  Load in command line paramters ###########
  print('Using Genome VCF file:', genome_vcf_file)
  print('Using GTF file: ', gtf_file)
  print('Using transcriptome file:', transcriptome_file)
  print('Using output file:', output_file)
  print('Outputting coordinate conversion:', coordinate_conversion_file)
  coords_f = open(coordinate_conversion_file, 'w+')
  coords_f.write('\t'.join(['CHR', 'GENOME_POS', 'TRANSCRIPT', 'TRANSCRIPT_POS', 'GREF', 'GALT', 'TREF', 'TALT']) + '\n')
  
  ########################################################################################################################


  ########################################################################################################################
  #######  Load GTF file in ###########
  print('Loading GTF file')
  ########################################################################################################################
  gtf_df = read_gtf(gtf_file)

  gtf_df = gtf_df.drop(['source', 'score', 'frame', 'gene_id', 'gene_type', 'gene_name', \
        'havana_gene', 'transcript_type', 'transcript_support_level', \
        'tag', 'level', 'transcript_name', 'havana_transcript', 'exon_id', \
        'ont'], axis=1)
  gtf_df = gtf_df[(gtf_df.feature == 'exon') | (gtf_df.feature == 'transcript')]
  ########################################################################################################################
  #######  Create transcript to exon lookup, which is a map from a transcript to its constituent exons ###########
  print('Creating transcript to exon map')
  ########################################################################################################################
  transcript_exon_lookup = {}

  curr_tid = None
  sw = 0
  sorted_df = gtf_df.sort_values(by=['transcript_id', 'feature', 'start'], ascending=[True, False, True])
  for i, row in sorted_df.iterrows():
    if curr_tid != row.transcript_id:
      curr_tid = row.transcript_id
      sw += 1
      exon_num = 1
      # if sw % 10000 == 0:
      #   print(sw)
        
  #     transcript_exon_lookup[curr_tid.split('_')[0]] = {
      transcript_exon_lookup[curr_tid] = {
        'exons': {},
        'range': (row.start, row.end),
        'chromosome': row.seqname.split('_')[0],
        'strand': row.strand,
        'tname': curr_tid
      }
    else:
  #     transcript_exon_lookup[curr_tid.split('_')[0]]['exons'][int(row.exon_number)] = (row.start, row.end)
      transcript_exon_lookup[curr_tid]['exons'][int(row.exon_number)] = (row.start, row.end)
      exon_num += 1
  ########################################################################################################################
  #######  Create map from chromosome to all the transcripts that it contains ###########
  print('Creating chromosome to transcript list map')
  ########################################################################################################################
  chr_to_transcript = {}

  for tname in transcript_exon_lookup:
    chr = transcript_exon_lookup[tname]['chromosome']
    if chr not in chr_to_transcript:
      chr_to_transcript[chr] = []
    
    chr_to_transcript[chr].append(transcript_exon_lookup[tname])

  for chr in chr_to_transcript:
    chr_to_transcript[chr] = sorted(chr_to_transcript[chr], key=lambda x: x['range'][0])

  ########################################################################################################################
  #######  Load in reference transcriptome to translate indels and double check conversion ###########
  print('Loading in reference transcriptome')
  ########################################################################################################################
  transcriptToRef = {}
  for tname, ref in transcript_generator(transcriptome_file):
    tname = tname.split('|')[0]
    transcriptToRef[tname] = ref

  
  ########################################################################################################################
  #######  Iterate through specified genomic vcf file, convert to transcriptomic file as we go ###########
  print('Beginning conversion of genome vcf file')
  ########################################################################################################################

  start_time = time.time()
  num_inconsistent = 0 # number of times our conversion doesn't match with the reference transcriptome
  num_consistent = 0 # number of times our conversion matches with the reference transcriptome
  num_conversion_right = 0
  num_conversion_wrong = 0
  strange = set([])
  curr_start = 0
  num_processed = 0
  thresh = .00747
  kill = False
  # thresh = .01028
  with gzip.open(genome_vcf_file) as genome_vcf:
      with open(output_file, 'w+') as transcriptome_vcf:
          start = False
          switch_gt = False
          sample_ind1 = -1
          sample_ind2 = -1
          sample_ind3 = -1
          for line in genome_vcf:
              line = str(line, 'utf-8')
              if num_processed > 0 and num_processed % 500000 == 0:
                  print('Processed', str(num_processed), 'records')
              # line = str(line, 'utf-8')
              if 'CHROM' in line:
                  start = True
                  transcriptome_vcf.write(line)
                  toks = line.split('\t')
                  #sample_ind1 = toks[5:].index('HG00405')
                  #sample_ind2 = toks[5:].index('NA18484')
                  #sample_ind3 = toks[5:].index('NA19185')
                  # print(sample_ind, toks[5:][sample_ind], file=sys.stderr)
                  continue

              if not start:
                  if '##contig' not in line:
                      transcriptome_vcf.write(line)
                  continue

              toks = line.split('\t')
              toks[7] = '.'
              chr, gpos, id, ref, alt = toks[:5]
              sanity = id.strip().split(':')
              if len(sanity) == 4:
                  san_ref = sanity[2]
                  san_alt = sanity[3]
                  if ref != san_ref or alt != san_alt:
                      assert(len(ref) > len(alt) or len(alt) > len(ref))
              gpos = int(gpos)

              if '<' in alt or ',' in alt: # skip multiallelic and weird variants for now
                  continue

              # add option for snp only?
              # figure out if we are to introduce a switch error here

              # for this chromosome and genomic position, find all transcripts that have this snp
              transcript_list = []
              transcripts_ends = []
              # Need to add a reset button for curr start
              curr_start, transcripts = get_transcript_containing_pos(chr, gpos, chr_to_transcript, curr_start)
              if len(ref) > len(alt):
                  _, transcripts_ends = get_transcript_containing_pos(chr, gpos + (len(ref) - len(alt)), chr_to_transcript, curr_start)
                  # if this position is in transcript ends but not transcript, then it overlaps with the start 
              in_orig = {}
              for k in range(len(transcripts)):
                  if transcripts[k]['strand'] == '+':
                      tpos = convert_bp_positive(transcripts[k], gpos)
                  else:
                      tpos = convert_bp_negative(transcripts[k], gpos)

                  tname = transcripts[k]['tname']
                  transcript_list.append((tname, tpos, transcripts[k]['strand']))
                  in_orig[tname] = k

              in_new = {}
              if len(ref) > len(alt):
                  for k in range(len(transcripts_ends)):
                      tname = transcripts_ends[k]['tname']
                      in_new[tname] = k
                      if tname not in in_orig:
                          # This deletion runs through the start of the transcript. add tname, tpos, strand
                          if transcripts_ends[k]['strand'] == '+':
                              tpos = convert_bp_positive(transcripts_ends[k], gpos + (len(ref) - len(alt)))
                              tpos = tpos - (len(ref) - len(alt))
                              # tpos needs to be minus one-d below.
                              # keep in mind if this is zero then the deletion only affects the first char
                              transcript_list.append((tname, tpos, transcripts_ends[k]['strand']))
                          else:
                              tpos = convert_bp_negative(transcripts_ends[k], gpos + (len(ref) - len(alt)))
                              tpos = tpos + (len(ref) - len(alt))
                              transcript_list.append((tname, tpos, transcripts_ends[k]['strand']))
                              # want to translate this to a deletion of the last n characters
              
              # Write line into new file for each transcript that it affects
              for i in range(len(transcript_list)):
                  tname, tpos, strand = transcript_list[i]
                  # A bailout, in case things are unexpected
                  skipThis = False
                  if tname not in transcriptToRef: # skip variants for regions that aren't transcribed...
                      strange.add((chr, gpos, tname, tpos))
                      continue   

                  if strand == '-':
                      if len(ref) < len(alt): # this is an insertion
                          assert(len(ref) == 1)
                          tref = transcriptToRef[tname][tpos - 2]
                          talt =  tref + complement_sequence(alt[1:])
                          tpos = tpos - 2 + 1

                          # If this is on the edge of an intron, skip it
                          next_pos = convert_bp_negative(transcripts[in_orig[tname]], gpos + 1)
                          if next_pos == len(transcriptToRef[tname]) + 1:
                              skipThis = True

                      elif len(ref) > len(alt): # this is a deletion
                          assert(len(alt) == 1)
                          # this deletion starts in an intron (before the transcript) and ends in an exon
                          if tpos - 1 >= len(transcriptToRef[tname]):
                              tpos = tpos - 2 - (len(ref) - len(alt))
                              tref = transcriptToRef[tname][tpos:]
                              talt = tref[0]
                              tpos = tpos + 1
                          else:
                              # this deletion starts in an exon (end of transcript) and runs off end
                              if tpos - 2 - (len(ref) - len(alt)) < 0:
                                  tref = '.'
                                  talt = str(len(ref) - len(alt))
                                  tpos = tpos - 2 - (len(ref) - len(alt)) + 1
                              # this deletion starts in an intron and ends in an exon alternative splice
                              elif tname not in in_orig and tname in in_new:
                                  # truncate neg for this
                                  tpos = tpos - (len(ref) - len(alt)) - 1 # check for edge case
                                  tpos_end = truncate_neg(transcripts_ends[in_new[tname]], gpos, ref, alt)
                                  tref = transcriptToRef[tname][tpos - 1:tpos_end]
                                  talt = tref[0]
                              # this deletion starts in an exon and ends in an intron alternative splice
                              elif tname in in_orig and tname not in in_new:
                                  # set anchor neg and reverse complement and shit
                                  t_anchor = get_anchor_neg(transcripts[in_orig[tname]], gpos + (len(ref) - len(alt)), ref, alt)
                                  if t_anchor is None:
                                      skipThis = True
                                  else:
                                      tpos, tpos_end = t_anchor
                                      tref = transcriptToRef[tname][tpos-1:tpos_end-1]
                                      talt = tref[0]
                                      if len(tref) == len(talt):
                                          skipThis = True
                              elif tname in in_orig and tname in in_new:
                                  tpos_end = convert_bp_negative(transcripts[in_orig[tname]], gpos + (len(ref) - len(alt)))
                                  if tpos - tpos_end != len(ref) - len(alt):
                                      print('Deletion over alternative splicing site in ', tname)
                                      tref = transcriptToRef[tname][tpos_end - 2 : tpos - 1]
                                      talt = tref[0]
                                      tpos = tpos_end - 1
                                  else:
                                      talt = transcriptToRef[tname][tpos - 2 - (len(ref) - len(alt))]
                                      tref = talt + complement_sequence(ref[1:])
                                  
                                      if tpos - 1 < len(transcriptToRef[tname]):
                                          if transcriptToRef[tname][tpos - 1] == complement_nucleotide(ref[0]):
                  #                             print('good', chr, gpos, ref, alt, tname, tpos, strand)
                                              num_consistent += 1
                                          else:
                                              print('bad', chr, gpos, ref, alt, tname, tpos, strand)
                                              num_inconsistent += 1 
                                      tpos = tpos - 2 - (len(ref) - len(alt)) + 1

                      else:
                          tref = complement_nucleotide(ref)
                          talt = complement_nucleotide(alt)
                  else:
                      if len(alt) > len(ref): # this is an insertion
                          tref = ref
                          talt = alt
                          # If this is on the edge of an intron, skip it?
                          next_pos = convert_bp_positive(transcripts[in_orig[tname]], gpos + 1)
                          if next_pos is None:
                              skipThis = True
                      elif len(ref) > len(alt): # this is a deletion. make sure it doesn't run off the length of the transcript
                          if tpos <= 0: # this deletion cuts into the start of this transcript
                              tref = '.'
                              talt = str(len(ref) - len(alt))
                          elif tname not in in_orig and tname in in_new:
                              # Go to previous exon (if it exists), take last element there as anchor
                              # Count how many transcriptomic positions this intron would take
                              # Do the thing
                              t_anchor = get_anchor_pos(transcripts_ends[in_new[tname]], gpos, ref, alt)
                              tpos, tpos_end = t_anchor
                              tref = transcriptToRef[tname][tpos-1 : tpos_end]
                              talt = tref[0]
                          elif tname in in_orig and tname not in in_new:
                              # if the end of this runs off an exon, truncate it
                              tpos_end = truncate_pos(transcripts[in_orig[tname]], gpos, ref, alt)
                              tref = transcriptToRef[tname][tpos - 1: tpos_end]
                              talt = tref[0]
                              if len(tref) == len(talt):
                                  skipThis = True
  #                                 print(':O')
                          elif tname in in_orig and tname in in_new: # both start and end positions of deletion are in bounds
                              if tpos > 0:
                                  if transcriptToRef[tname][tpos - 1] == ref[0]:
                                      num_consistent += 1
                                  else:
                                      print('bad', chr, gpos, ref, alt, tname, tpos, strand)
                                      num_inconsistent += 1
                              if tpos - 1 + (len(ref) - len(alt)) < len(transcriptToRef[tname]):
                                  tpos_end = convert_bp_positive(transcripts[in_orig[tname]], gpos + (len(ref) - len(alt)))
                                  if tpos_end - tpos != len(ref) - len(alt):
                                      print('Deletion over alternative splicing site in ', tname)
                                      tref = transcriptToRef[tname][tpos - 1 : tpos_end]
                                      talt = tref[0]
                                  else:
                                      tref = ref
                                      talt = alt
                              else:
                                  tref = ref[:len(transcriptToRef[tname])-(tpos - 1)] # transcriptToRef[tpos - 1:]
                                  if tref == transcriptToRef[tname][tpos-1:]:
                                      num_conversion_right += 1
                                  else:
                                      num_conversion_wrong += 1
                                      tref = transcriptToRef[tname][tpos-1:]
                                  talt = tref[0]

                                  if len(tref) == len(talt): # the insertion goes off the end, and the remaining portion does nothing
                                      skipThis = True

                      else:
                          tref = ref
                          talt = alt

                  if tref is None or talt is None:
                      print('skipping: ', chr, gpos, ref, alt, "transcript: ", tref, talt)
                      skipThis = True
                  
                  if not skipThis:
                      transcriptome_vcf.write('\t'.join([tname, str(tpos), id + 'v' + str(i), tref, talt] + toks[5:]))
                      # ['CHR', 'GENOME_POS', 'TRANSCRIPT', 'TRANSCRIPT_POS', 'GREF', 'GALT', 'TREF', 'TALT']
                      coords_f.write('\t'.join([chr, str(gpos), tname, str(tpos), ref, alt, tref, talt]) + '\n')

              num_processed += 1

  print('Total time taken for conversion: ', time.time() - start_time)
  print('Number of SNPs consistent with conversion: ', num_consistent)
  print('Number of SNPs inconsistent with conversion: ', num_inconsistent)
  print('Number of SNPs where the conversion matches up', num_conversion_right)
  print('Number of SNPs where the conversion does not match up', num_conversion_wrong)
  print('Positions converted into transcripts that are not found in reference transcriptome:')
  print(strange)

def parse_arguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-v','--genome-vcf-file', required=True, type=str)
  parser.add_argument('-g','--gtf-file', required=True, type=str)
  parser.add_argument('-t','--transcriptome', required=True, type=str)
  parser.add_argument('-o','--output', required=True, type=str)
  parser.add_argument('-c','--coordinate-output', required=True, type=str)
  return parser.parse_args()

if __name__ == "__main__":
  args = parse_arguments()
  main(args.genome_vcf_file, args.gtf_file, args.transcriptome, args.output, args.coordinate_output)
