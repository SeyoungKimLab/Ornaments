import argparse
import gzip
import sys
import numpy as np
from itertools import chain, combinations
from enum import Enum

class Utils:
  @classmethod
  def is_het(cls, gt):
    return gt == '0|1' or gt == '1|0' or gt == '0/1' or gt == '1/0'

  @classmethod
  def powerset(cls, iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

  @classmethod
  def ntuples(cls, l, het_vars_inds):
    if l == 0:
      yield []
    else:
      # for i in range(len(het_vars_inds)):
      for i in het_vars_inds[l - 1]:
        for L in Utils.ntuples(l - 1, het_vars_inds):
          yield [i] + L

  @classmethod
  def ntuple_variants(cls, dict_):
    arrs = [dict_[k] for k in dict_]
    for ntuple in Utils.ntuples(len(arrs), arrs):
      yield {k: ntuple[-i-1] for i, k in enumerate(dict_)}

  @classmethod
  def set_overlap(cls, range1, range2):
    return max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]) + 1)

  @classmethod
  def transcript_generator(cls, file):
    with open(file) as transcriptome:
      ref = ''
      curr_tid = None
      for line in transcriptome:
        if '>' in line:
        # The current transcript has been finished is done
          if len(ref) > 0:
            yield (curr_tid, ref)

          curr_tid = line.split('>')[1].strip()
          ref = ''
        else:
          ref += line.strip()

      yield (curr_tid, ref)

def get_snps(curr_range, snp_index):
  if len(snp_index) == 0:
    return []
  snps_in_range = []
  for ind in snp_index:
    if ind >= curr_range[0] and ind <= curr_range[1]:
      snps_in_range.append(ind)
  return snps_in_range

# If there is an overlap, we have to extend the region that we are considering
def calculate_extract_for_transcript(snp_list, ins_dict={}, del_dict={}, k=31):
  i = 0
  curr_run = 1
  added_bp = 0
  indel_dict = {**ins_dict,  **del_dict}
  indel_list = list(indel_dict.keys())
  indel_positions = set(indel_list)
  snp_list = snp_list + indel_list
  snp_list = snp_list + [100000000]
  snp_list = sorted(snp_list)
  num_seq = 0
  snp_partitions = []
  ranges = []
  active_bonus = 0
  while i < len(snp_list) - 1:
    if snp_list[i] in indel_positions and type(indel_dict[snp_list[i]][-1]) is int:
      active_bonus = np.max(np.abs(indel_dict[snp_list[i]]))
    else:
      active_bonus = 0
  
    if snp_list[i+1] - snp_list[i] < k + active_bonus:
      curr_run += 1
      i += 1
      continue
    else:
      grouped_snps = snp_list[i-curr_run+1:i+1] 
      curr_run = 1
      i += 1
      snp_partitions.append(grouped_snps)
      ranges.append((grouped_snps[0] - k, grouped_snps[-1] + k + active_bonus))
      active_bonus = 0
  
  return (snp_partitions, ranges)

def contained_snps_and_indels(snps, ins, dels, range_, lbound):
  ret_snps = {}
  ret_ins = {}
  ret_dels = {}

  for s in snps:
    if s >= range_[0] and s <= range_[1]:
      ret_snps[s - lbound] = snps[s]

  for i in ins:
    if i >= range_[0] and i <= range_[1]:
      ret_ins[i - lbound] = ins[i]

  for d in dels:
    if d >= range_[0] and d <= range_[1]:
      ret_dels[d - lbound] = dels[d]
    else:
      l, r = d, d + np.max(np.abs(dels[d]))
      if r >= range_[0] and r <= range_[1]:
        ret_dels[d - lbound] = dels[d]

  return (ret_snps, ret_ins, ret_dels)

def contained_snps(snps, range_, indel_dict, bonus=0):
  to_ret = []
  for s in snps:
    if s not in indel_dict:
      if s >= range_[0] and s <= range_[1]:
        to_ret.append(s)
    else:
      if s >= range_[0] and s < range_[1]:
        to_ret.append(s)
      else:
        if type(indel_dict[s][0]) is int: # deletion
          l, r = s, s + np.max(np.abs(indel_dict[s]))
          tup = (l, r)
          # If there is an overlap in the range, then append this snp
          if Utils.set_overlap((l + 1, r), range_) > 0:
            to_ret.append(s)
  return to_ret
  # return [s for s in snps if s >= range_[0] and s <= range_[1]]

def get_insertion_adjustment(seq, pos, indel):
  rel_area = seq[pos:pos+len(indel)]
  num_shift = 0
  for i in range(len(rel_area)):
    if rel_area[i] == indel[i]:
      num_shift += 1
    else:
      break
  return num_shift

def get_insertion_end_adjustment(seq, pos, indel):
  rel_area = seq[pos-len(indel)+1:pos+1]
  num_shift = 0
  for i in range(len(rel_area)):
    if rel_area[len(rel_area) - i - 1] == indel[len(indel) - i - 1]:
      num_shift += 1
    else:
      break
  return num_shift

def get_deletion_adjustment(seq, pos, indel):
  rel_area = seq[pos+1:]
  del_area = seq[pos+indel+1:]
  num_shift = 0
  for i in range(min(len(rel_area), len(del_area))):
    if rel_area[i] == del_area[i]:
      num_shift += 1
    else:
      break
  return num_shift + 1

# If a deletion would result in a kmer that also occurs in the ref haplotype
# Shift the end to exclude this kmer
def get_deletion_end_adjustment(seq, pos, indel):
  rel_area = seq[:pos+indel+1]
  prev_area = seq[:pos+1]
  num_shift = 0
  for i in range(min(len(rel_area), len(prev_area))):
    if rel_area[len(rel_area) - i - 1] == prev_area[len(prev_area) - i - 1]:
      num_shift += 1
    else:
      break
  return num_shift

def get_bonus(current_snps, indel_positions, indel_dict, i):
  pre_active_bonus = 0
  active_indel_stack = []
  for curr_snp in current_snps:
    if curr_snp in indel_positions:
      # if type(indel_dict[curr_snp][0]) is int and curr_snp >= i: # we have a deletion
      if type(indel_dict[curr_snp][0]) is int: # we have a deletion
        # If this deletion is not within another deletion's active range
        curr_real_range = ((curr_snp, curr_snp + np.max(np.abs(indel_dict[curr_snp]))))
        if len(active_indel_stack) > 0:
          l, r = active_indel_stack[-1]
          ovlp = max(0, min(r, curr_real_range[1]) - max(l, curr_real_range[0]) + 1)
          if ovlp != 0:
            pre_active_bonus -= (ovlp - 1)

        pre_active_bonus += np.max(np.abs(indel_dict[curr_snp]))
        if curr_snp < i:
          pre_active_bonus -= (i - curr_snp - 1)
        active_indel_stack.append(curr_real_range)
  return pre_active_bonus

def calculate_subsequences_new(sequence, o_snps, ins_dict={}, del_dict={}, k=31):
  indel_dict = {**ins_dict,  **del_dict}
  indel_positions = set(list(indel_dict.keys()))
  snps = o_snps + list(indel_dict.keys())
  snps = sorted(snps)
  start_shift = 0
  # raise issue if variant doesn't match reference? just don't include it in vcf?
  if len(snps) == 1: # REWORK IN CASE THIS IS AN INDEL
    start_range, end_range = snps[0] - k + 1, snps[0] + k - 1
    if snps[0] in indel_dict and type(indel_dict[snps[0]][0]) is int: # deletion
      start_range += 1
      end_range += np.max(np.abs(indel_dict[snps[0]]))
    elif snps[0] in indel_dict:
      start_range += 1

    if start_range < 0:
      start_range = 0
    if end_range > len(sequence):
      end_range = len(sequence)

    return [(start_range, end_range)], \
           [[s for s in snps if s not in indel_dict]], \
           [{s: indel_dict[s] for s in snps if s in indel_dict}], \
           [sequence[start_range:end_range+1]]
  
  # By definition, adjacent snps will be less than k distance from each other
  # However, further adjacent snps can also be less than k distance from each other and require special treatment
  current_range = (snps[0] - k + 1, snps[0]) # inclusive range
  current_snps = contained_snps(snps, current_range, indel_dict)
  start_range = snps[0] - k + 1
  
  ranges = []
  snps_for_ranges = []
  indels_for_ranges = []
  # Only apply one bonus if deletions overlap
  # Actually, how do we deal with overlaps in general?
  for i in range(snps[0] - k + 2, snps[-1] + 1):
    pre_active_bonus = get_bonus(current_snps, indel_positions, indel_dict, i)
    test_snps = contained_snps(snps, (i, i + k - 1), indel_dict, 0)
    # If there are any new indels in here, and using the active bonus from them would include
    # even more snps, then those snps are contained within the deletion
    curr_active_bonus = get_bonus(test_snps, indel_positions, indel_dict, i)
    new_snps = contained_snps(snps, (i, i + k + curr_active_bonus - 1), indel_dict, curr_active_bonus)
    if new_snps != current_snps:
      end_range = (i + k - 2) + pre_active_bonus
    
      if start_range < 0:
        start_range = 0
      if end_range > len(sequence):
        end_range = len(sequence)
    
      all_in = True
      check_contained = contained_snps(snps, (start_range, end_range), indel_dict, pre_active_bonus)
      for curr_snp in current_snps:
        if curr_snp not in check_contained:
          all_in = False
    
      if all_in:
        if end_range - start_range >= k - 1 and len(current_snps) > 0:
          ranges.append((start_range, end_range))
          snps_for_ranges.append([c for c in current_snps if c not in indel_dict])
          indels_for_ranges.append({c: indel_dict[c] for c in current_snps if c in indel_dict})

      start_range = i
    current_snps = new_snps
  
  # Will this always be one thing? What if last two snps are right next to each other
  if start_range < 0:
    start_range = 0
  
  curr_active_bonus = get_bonus(current_snps, indel_dict, indel_dict, i)
  if i + k + 1 - start_range >= k - 1: # make sure the range is actually valid
    ranges.append((start_range, i + k + curr_active_bonus - 1))
    snps_for_ranges.append([c for c in current_snps if c not in indel_dict])
    indels_for_ranges.append({c: indel_dict[c] for c in current_snps if c in indel_dict})

  return consolidate_ranges(ranges, snps_for_ranges, indels_for_ranges, sequence)

def snps_inds_are_in_range(snps, indels, range_):
  to_ret = True
  for s in snps:
    if s < range_[0] or s > range_[1]:
      to_ret = False
  for i in indels:
    if i < range_[0] or i >= range_[1]:
      to_ret = False
  return to_ret

# Ranges are sorted
def consolidate_ranges(ranges, snps_for_ranges, indels_for_ranges, sequence):
  if len(ranges) == 1:
    return ranges, snps_for_ranges, indels_for_ranges, [sequence[l:r+1] for (l, r) in ranges]

  new_ranges = []
  new_snps = []
  new_indels = []
  new_seqs = []
  i = 0
  cons_last = False
  while i < len(ranges) - 1:
    cl, cr = ranges[i]
    nl, nr = ranges[i+1]
    ns, ni = snps_for_ranges[i+1], indels_for_ranges[i+1]
    cs, ci = snps_for_ranges[i], indels_for_ranges[i]
    if snps_inds_are_in_range(ns, ni, (cl, cr)) and snps_inds_are_in_range(cs, ci, (nl, nr)):
      new_ranges.append((min(nl, cl), max(nr, cr)))
      new_snps.append(list(set(cs) | set(ns)))
      new_indels.append({**ci, **ni})
      if i == len(ranges) - 2:
        cons_last = False
      else:
        cons_last = True
      i += 2
    else:
      new_ranges.append(ranges[i])
      new_snps.append(snps_for_ranges[i])
      new_indels.append(indels_for_ranges[i])
      cons_last = True
      i += 1
  if cons_last:
    new_ranges.append(ranges[-1])
    new_snps.append(snps_for_ranges[-1])
    new_indels.append(indels_for_ranges[-1])

  return new_ranges, new_snps, new_indels, [sequence[l:r+1] for (l, r) in new_ranges]



def calculate_subsequences(sequence, o_snps, indel_dict={}, k=31):
  indel_positions = set(list(indel_dict.keys()))
  snps = o_snps + list(indel_dict.keys())
  snps = sorted(snps)
  if len(snps) == 1: # REWORK IN CASE THIS IS AN INDEL
    start_range, end_range = snps[0] - k + 1, snps[0] + k - 1
    start_shift = 0
    end_shift = 0
    if snps[0] in indel_dict and type(indel_dict[snps[0]]) is int: # deletion
        start_shift = get_deletion_adjustment(sequence, snps[0], abs(indel_dict[snps[0]]))
        end_shift = get_deletion_end_adjustment(sequence, snps[0], abs(indel_dict[snps[0]]))
        start_range += start_shift
        end_range += abs(indel_dict[snps[0]])
        end_range -= end_shift
    
    elif snps[0] in indel_dict and type(indel_dict[snps[0]]) is not int: # insertion
      # Shift the range by however much this matches the reference sequence at that position
      start_shift = get_insertion_adjustment(sequence, snps[0], indel_dict[snps[0]])
      end_shift = get_insertion_end_adjustment(sequence, snps[0], indel_dict[snps[0]])
      start_range += start_shift
      end_range -= end_shift

    if start_range < 0:
      start_range = 0
    if end_range > len(sequence):
      end_range = len(sequence)

    return [(start_range, end_range)], \
           [[s for s in snps if s not in indel_dict]], \
           [{s - start_shift: indel_dict[s] for s in snps if s in indel_dict}], \
           [{s: indel_dict[s] for s in snps if s in indel_dict}], \
           [sequence[start_range:end_range+1]]
  
  # By definition, adjacent snps will be less than k distance from each other
  # However, further adjacent snps can also be less than k distance from each other and require special treatment
  # Ignore indels that are within k distance of SNPs
  current_range = (snps[0] - k + 1, snps[0]) # inclusive range
  current_snps = contained_snps(snps, current_range)
  start_range = snps[0] - k + 1
  
  ranges = []
  snps_for_ranges = []
  indels_for_ranges = []
  duds = []
  for i in range(snps[0] - k + 2, snps[-1] + 1):
    pre_active_bonus = 0
    for curr_snp in current_snps:
      if curr_snp in indel_positions:
        if type(indel_dict[curr_snp]) is int and indel_dict[curr_snp] < 0: # we have a deletion
          pre_active_bonus += abs(indel_dict[curr_snp])
        
    new_snps = contained_snps(snps, (i, i + k + pre_active_bonus - 1), pre_active_bonus)
    if new_snps != current_snps:
      end_range = i + k - 2
      # calculate sum total of deletions in the current range
      active_bonus = 0
      for curr_snp in current_snps:
        if curr_snp in indel_positions:
          if type(indel_dict[curr_snp]) is int and indel_dict[curr_snp] < 0: # we have a deletion
            active_bonus += abs(indel_dict[curr_snp])
      end_range += active_bonus
    
      if start_range < 0:
        start_range = 0
      if end_range > len(sequence):
        end_range = len(sequence)
    
      all_in = True
      for curr_snp in current_snps:
        if curr_snp < start_range or curr_snp > end_range:
          all_in = False
    
      if all_in:
        if end_range - start_range >= k - 1 and len(current_snps) > 0:
          ranges.append((start_range, end_range))
          snps_for_ranges.append([c for c in current_snps if c not in indel_dict])
          indels_for_ranges.append({c: indel_dict[c] for c in current_snps if c in indel_dict})
          duds.append({})

      start_range = i
    current_snps = new_snps
  
  # Will this always be one thing? What if last two snps are right next to each other
  if start_range < 0:
    start_range = 0
  
  active_bonus = 0
  for curr_snp in current_snps:
    if curr_snp in indel_positions:
      if type(indel_dict[curr_snp]) is int and indel_dict[curr_snp] < 0: # we have a deletion
        active_bonus += abs(indel_dict[curr_snp])
  
  if i + k + 1 - start_range >= k - 1: # make sure the range is actually valid
    ranges.append((start_range, i + k + active_bonus - 1))
    snps_for_ranges.append([c for c in current_snps if c not in indel_dict])
    indels_for_ranges.append({c: indel_dict[c] for c in current_snps if c in indel_dict})
    duds.append({})

  s = [sequence[l:r+1] for (l, r) in ranges]
  return ranges, snps_for_ranges, indels_for_ranges, duds, s

def is_homo(gt):
  return gt == '0|0' or gt == '1|1' or gt == '0' or gt == '1' or gt == '0/0' or gt == '1/1'

def get_haplotype_for_sample_p(tsnps, snps, samples):
  curr_combos = set([])
  print(tsnps)
  for s in snps:
    tsnps[s] = tsnps[s].split('\t')
  for i in samples:
    left_genos = []
    right_genos = []
    for s in snps:
      gt = tsnps[s][i]
      if len(gt) == 3:
        if gt[0] == '1':
          left_genos.append(s)
        if gt[-1] == '1':
          right_genos.append(s)
      else:
        assert(len(gt) == 1)
        if gt[0] == '1':
          left_genos.append(s)
          right_genos.append(s)
    curr_combos.add(frozenset(left_genos))
    curr_combos.add(frozenset(right_genos))
  for s in snps:
    tsnps[s] = '\t'.join(tsnps[s])
  return list(curr_combos)

def get_haplotype_for_sample_up(tsnps, snps, samples):
  curr_combos = set([])
  for s in snps:
    tsnps[s] = tsnps[s].split('\t')
  for i in samples:
    is_curr_het = []
    is_curr_homo_alt = []
    for s in snps:
      gt = tsnps[s][i]
      if not is_homo(gt):
        is_curr_het.append(s)
      else:
        if '1' in gt:
          is_curr_homo_alt.append(s)
    pset = list(Utils.powerset(is_curr_het))
    for arr in pset:
      curr_combos.add(frozenset(set(arr) | set(is_curr_homo_alt)))
  for s in snps:
    tsnps[s] = '\t'.join(tsnps[s])
  return list(curr_combos)

def apply_some_snps_and_indels_and_edges(
    seq, snpDict, insDict, delDict, edgeDict,
    otherLeftSnps, otherLeftIns, otherLeftDels):
  to_ret = seq
  transformed_coords_ins = {i: i for i in insDict}
  squelched_snps = set([])
  squelched_ins = set([])
  squelched_dels = set([])
  squelched_other_snps = set([])
  squelched_other_ins = set([])
  squelched_other_dels = set([])

  for e in edgeDict:
    del_num = e + edgeDict[e]
    if del_num >= 0:
      to_ret = '-' * (del_num + 1) + to_ret[del_num+1:]

  for d in delDict:
    len_del = abs(delDict[d])
    ind = d
    if to_ret[ind] == '-' or to_ret[ind + len_del] == '-':
      squelched_dels.add(d)
      continue
    else:
      to_ret = to_ret[:ind+1] + '-' * len_del + to_ret[ind+len_del+1:]

  for s in snpDict:
    ind = s
    if to_ret[ind] == '-':
      squelched_snps.add(ind)
    else:
      to_ret = to_ret[:ind] + snpDict[ind] + to_ret[(ind+1):]

  for d in otherDels:
    len_del - abs(otherDels[d])
    ind = d
    if to_ret[ind] == '-' or to_ret[ind + len_del] == '-':
      squelched_other_dels.add(d)

  for s in otherSnps:
    ind = s
    if to_ret[ind] == '-':
      squelched_other_snps.add(ind)

  for i in otherIns:
    ind = i
    if to_ret[ind] == '-':
      squelched_other_ins.add(i)

  for i in insDict:
    ind = transformed_coords_ins[i]
    len_ins = len(insDict[i]) - 1
    if to_ret[ind] == '-':
      squelched_ins.add(i)
      continue
    else:
      to_ret = to_ret[:ind+1] + insDict[i][1:] + to_ret[ind+1:]
      for j in insDict:
        if transformed_coords_ins[j] > ind:
          transformed_coords_ins[j] += len_ins

  return (to_ret, \
      {s for s in otherSnps if s not in squelched_other_snps}, \
      {i for i in otherIns if i not in squelched_other_ins}, \
      {d for d in otherDels if d not in squelched_other_dels})

def apply_dem_snps_and_indels_and_edges(seq, snpDict, insDict, delDict, edgeDict):
  to_ret = seq
  transformed_coords_snps = {s: s for s in snpDict}
  transformed_coords_ins = {i: i for i in insDict}
  transformed_coords_ins_real = {i: i for i in insDict}
  transformed_coords_dels = {d: d for d in delDict}
  list_pos_transform = np.arange(0, len(seq), 1)

  squelched_snps = set([])
  squelched_ins = set([])
  squelched_dels = set([])
  squelched_edges = set([])

  for e in edgeDict:
    del_num = e + edgeDict[e]
    if to_ret[0] == '-':
      squelched_edges.add(e)
      continue
    if del_num >= 0:
      to_ret = '-' * (del_num + 1) + to_ret[del_num+1:]
      list_pos_transform[del_num:] -= del_num

  for d in delDict:
    len_del = abs(delDict[d])
    ind = d
    if to_ret[ind] == '-' or (ind + len_del < len(to_ret) and to_ret[ind + len_del] == '-'):
      squelched_dels.add(d)
      continue
    if len_del == 0:
      continue
    else:
      to_ret = to_ret[:ind+1] + '-' * len_del + to_ret[ind+len_del+1:]
      list_pos_transform[ind]
      for a in transformed_coords_snps:
        trans_ind = transformed_coords_snps[a]
        if trans_ind > d + len_del:
          transformed_coords_snps[a] -= len_del
      for a in transformed_coords_ins:
        trans_ind = transformed_coords_ins[a]
        if trans_ind > d + len_del:
          transformed_coords_ins[a] -= len_del
      for b in transformed_coords_dels:
        trans_ind = transformed_coords_dels[b]
        if trans_ind > d + len_del:
          transformed_coords_dels[b] -= len_del

  for s in snpDict:
    ind = s
    if ind >= len(to_ret):
      continue
    if to_ret[ind] == '-':
      squelched_snps.add(s)
    else:
      to_ret = to_ret[:ind] + snpDict[ind] + to_ret[(ind+1):]

  for i in insDict:
    ind = transformed_coords_ins_real[i]
    len_ins = len(insDict[i]) - 1
    if to_ret[ind] == '-':
      squelched_ins.add(i)
      continue
    if len_ins == 0:
      continue
    else:
      to_ret = to_ret[:ind+1] + insDict[i][1:] + to_ret[ind+1:]
      for j in insDict:
        if transformed_coords_ins_real[j] > ind:
          transformed_coords_ins_real[j] += len_ins

      for a in transformed_coords_snps:
        trans_ind = transformed_coords_snps[a]
        if trans_ind > ind:
          transformed_coords_snps[a] += len_ins

      for a in transformed_coords_ins:
        trans_ind = transformed_coords_ins[a]
        if trans_ind > ind:
          transformed_coords_ins[a] += len_ins

      for b in transformed_coords_dels:
        trans_ind = transformed_coords_dels[b]
        if trans_ind > ind:
          transformed_coords_dels[b] += len_ins

  return (to_ret.replace('-',''),
      {s: transformed_coords_snps[s] for s in snpDict if s not in squelched_snps}, 
      {i: transformed_coords_ins[i] for i in insDict if i not in squelched_ins},
      {d: transformed_coords_dels[d] for d in delDict if d not in squelched_dels},
      squelched_snps,
      squelched_ins,
      squelched_dels,
      squelched_edges)

# TODO: Need to remove squelched shit from this
def get_ornament_name(indexer, combo_snps, combo_ins, combo_dels, combo_edges, lbound, homos, squelched):
  homo_snps, homo_ins, homo_dels, homo_edges = homos[0], homos[1], homos[2], homos[3]
  squelched_snps, squelched_ins, squelched_dels, squelched_edges = \
    squelched[0], squelched[1], squelched[2], squelched[3]
  to_ret = '_refsnp_('
  grouper = {'r': [], 'a_0': []}
  num_variants = 0
  for s in combo_snps:
    if s in homo_snps or s in squelched_snps:
      continue
    t = indexer.get_snp(s + lbound, combo_snps[s])
    if t not in grouper:
      grouper[t] = []
    grouper[t].append(s + lbound)
    num_variants += 1
  for i in combo_ins:
    if i in squelched_ins:
      continue
    t = indexer.get_ins(i + lbound, combo_ins[i])
    if t not in grouper:
      grouper[t] = []
    grouper[t].append(i + lbound)
    num_variants += 1
  for d in combo_dels:
    if d in squelched_dels:
      continue
    t = indexer.get_del(d + lbound, combo_dels[d])
    if t not in grouper:
      grouper[t] = []
    grouper[t].append(d + lbound)
    num_variants += 1

  for e in combo_edges:
    if combo_edges[e] != 0: # only include ref for edges, since alt won't generate A.S. k-mers
      continue
    t = indexer.get_del(e + lbound, combo_edges[e])
    if t not in grouper:
      grouper[t] = []
    grouper[t].append(e + lbound)
    num_variants += 1

  for s in grouper['r']:
    to_ret += (str(s) + ',')
  if to_ret[-1] == ',':
    to_ret = to_ret[:-1]
  to_ret += ')|altsnp0_('
  for s in grouper['a_0']:
    to_ret += (str(s) + ',')
  if to_ret[-1] == ',':
    to_ret = to_ret[:-1]
  to_ret += ')|'

  for k in sorted(grouper):
    if k == 'r' or k == 'a_0':
      continue
    num = k[-1]
    to_ret += 'altsnp' + str(num) + '_('
    for s in grouper[k]:
      to_ret += (str(s) + ',')
    if to_ret[-1] == ',':
      to_ret = to_ret[:-1]
    to_ret += ')|'

  return to_ret, num_variants

# Need to deal with multiallelic sites though, or sites with same position but different allele
def apply_snps_indels(
    kmer, k, indexer,
    snps, ins, dels, edges,
    homo_snps, homo_ins, homo_dels, homo_edges, lbound=0, phased=False): 
  # Need to consider both combos for homo alt indels. Otherwise, for het indels, we are set already?
  het_snps = {s: snps[s] for s in snps if s not in homo_snps}
  # het_edges = {e: edges[e] for e in edges if e not in homo_edges}
  if not phased:
    het_ins = {i: ins[i] for i in ins}
    het_dels = {d: dels[d] for d in dels}
    for i in het_ins:
      if len(het_ins[i]) >= 1:
        if het_ins[i][0][0] not in het_ins[i]:
          het_ins[i] = [het_ins[i][0][0]] + het_ins[i]

    for d in het_dels:
      if len(het_dels[d]) >= 1:
        if 0 not in het_dels[d]:
          het_dels[d] = [0] + het_dels[d]
  else:
    het_ins = {i: ins[i] for i in ins if i not in homo_ins}
    het_dels = {d: dels[d] for d in dels if d not in homo_dels}

  het_edges = {d: het_dels[d] for d in het_dels if d < 0}
  het_dels = {d: het_dels[d] for d in het_dels if d >= 0}

  for combo_snps in Utils.ntuple_variants(het_snps):
    for combo_ins in Utils.ntuple_variants(het_ins):
      for combo_dels in Utils.ntuple_variants(het_dels):
        for combo_edges in Utils.ntuple_variants(het_edges):
          true_combo_snps = {c: combo_snps[c] for c in combo_snps}
          true_combo_ins = {c: combo_ins[c] for c in combo_ins}
          true_combo_dels = {c: combo_dels[c] for c in combo_dels}

          for s in homo_snps:
            true_combo_snps[s] = snps[s][0]

          if phased:
            for i in homo_ins:
              true_combo_ins[i] = ins[i][0]
            for d in homo_dels:
              true_combo_dels[d] = dels[d][0]


          # presence in set means alternative, absence means reference
          to_ret, snpmap, insmap, delmap, ss, si, sd, se = apply_dem_snps_and_indels_and_edges(
              kmer,
              true_combo_snps,
              true_combo_ins,
              true_combo_dels,
              combo_edges)

          # Double check all variants qualify for this range and trim if not
          ltrim = 0
          rtrim = len(to_ret)
          ltrim_snps = set()
          rtrim_snps = set()
          ltrim_ins = set()
          rtrim_ins = set()
          ltrim_dels = set()
          rtrim_dels = set()

          for s in true_combo_snps:
            if s not in snpmap or s in homo_snps: # was squelched
              continue

            ind = snpmap[s]
            if ind + k < len(to_ret):
              rtrim = min(rtrim, ind + k)
              rtrim_snps.add(s)
            if ind + 1 - k > 0:
              ltrim = max(ltrim, ind + 1 - k)
              ltrim_snps.add(s)

          for i in true_combo_ins:
            if i not in insmap: # was squelched
              continue
            if phased and i in homo_ins:
              continue

            ind = insmap[i]
            if ind + len(true_combo_ins[i]) - 1 + k < len(to_ret):
              rtrim = min(rtrim, ind + len(true_combo_ins[i]) - 1 + k)
              rtrim_ins.add(i)
            if ind + 2 - k > 0:
              ltrim = max(ltrim, ind + 2 - k)
              ltrim_ins.add(i)

          for d in true_combo_dels:
            if d not in delmap: # was squelched
              continue
            if phased and d in homo_dels:
              continue

            ind = delmap[d]
            find = ind
            if true_combo_dels[d] == 0: # deletion is ref, not active
              rems = set(np.abs(het_dels[d])) - set([0])
              if len(rems) > 0:
                find += min(rems)

            if find + k < len(to_ret):
              rtrim = min(rtrim, find + k)
              rtrim_dels.add(d)
            if ind + 2 - k > 0:
              ltrim = max(ltrim, ind + 2 - k)
              ltrim_dels.add(d)

          ornament_name_preview, num_variants = get_ornament_name(
            indexer,
            combo_snps,
            combo_ins,
            combo_dels,
            combo_edges,
            lbound,
            (homo_snps, set(), set(), set()),
            (ss, si, sd, se))

          ornament_position = lbound + ltrim
          ornament_name = ornament_name_preview + str(ornament_position)
          if num_variants > 0:
            # print(lbound, combo_snps, combo_ins, combo_dels, ornament_name_preview, to_ret[ltrim:rtrim])
            yield ornament_name, to_ret[ltrim:rtrim]

          # Start trimming towards the left and see if there are any ornaments that we missed
          if ltrim > k:
            untrimmed_s = {s: combo_snps[s] for s in combo_snps if s not in ltrim_snps}
            untrimmed_i= {i: combo_ins[i] for i in combo_ins if i not in ltrim_ins}
            untrimmed_d = {d: combo_dels[d] for d in combo_dels if d not in ltrim_dels}
            untrimmed_e = {e: combo_edges[e] - ltrim for e in combo_edges if e + combo_edges[e] >= ltrim}
            ornament_name_preview, num_variants = get_ornament_name(
              indexer,
              untrimmed_s,
              untrimmed_i,
              untrimmed_d,
              untrimmed_e,
              lbound,
              (homo_snps, set(), set(), set()),
              (ss, si, sd, se))

            if num_variants > 0:
              ornament_position = lbound
              ornament_name = ornament_name_preview + str(ornament_position)
              yield ornament_name, to_ret[:ltrim]

          if len(to_ret) - rtrim > k:
            untrimmed_s = {s: combo_snps[s] for s in combo_snps if s not in rtrim_snps}
            untrimmed_i= {i: combo_ins[i] for i in combo_ins if i not in rtrim_ins}
            untrimmed_d = {d: combo_dels[d] for d in combo_dels if d not in rtrim_dels}
            ornament_name_preview, num_variants = get_ornament_name(
              indexer,
              untrimmed_s,
              untrimmed_i,
              untrimmed_d,
              {},
              lbound,
              (homo_snps, set(), set(), set()),
              (ss, si, sd, se))
            if num_variants > 0:
              ornament_position = lbound + rtrim
              ornament_name = ornament_name_preview + str(ornament_position)
              yield ornament_name, to_ret[rtrim:]


# Absence in the set right now can either mean squelched or refsnp. Separate these two out,
# and we have a complete implementation that applies variants (that are compatible) together

def apply_snps(kmer, tsnps, samples, phased, list_of_snps, list_of_indels={}):
  # print(kmer)
  kmers = []
  snp_list = list(list_of_snps.keys())
  if len(list_of_indels) > 0:
    list_of_insertions = {k: list_of_indels[k] for k in list_of_indels if type(list_of_indels[k]) is not int}
    list_of_deletions = {k: list_of_indels[k] for k in list_of_indels if type(list_of_indels[k]) is int}
  else:
    list_of_insertions = {}
    list_of_deletions = {}
  snp_list += list(list_of_insertions.keys())
  snp_list += list(list_of_deletions.keys())
  indel_positions = set(list(list_of_insertions.keys()))
  indel_positions |= set(list(list_of_deletions.keys()))

  to_remove = set([])
  # current behavior: if indels overlap, ignore them
  for p1 in indel_positions:
    for p2 in indel_positions:
      if p1 == p2:
        continue
      if p1 in list_of_insertions:
        eff1 = (p1, p1+len(list_of_insertions[p1]))
      elif p1 in list_of_deletions:
        eff1 = (p1, p1 - list_of_deletions[p1])
      if p2 in list_of_insertions:
        eff2 = (p2, p2+len(list_of_insertions[p2]))
      elif p2 in list_of_deletions:
        eff2 = (p2, p2 - list_of_deletions[p2])
      
      if (eff1[1] >= eff2[0] and eff1[1] <= eff2[1]) or (eff1[0] >= eff2[0] and eff1[0] <= eff2[1]):
        to_remove.add(p1)
        to_remove.add(p2)

  snp_list = [p for p in snp_list if p not in to_remove]

  # generate observable haplotypes for each sample
  if phased:
    pset = get_haplotype_for_sample_p(tsnps, snp_list, samples)
  else:
    pset = get_haplotype_for_sample_up(tsnps, snp_list, samples)
  
  for curr_snps in pset:
    new_kmer = kmer
    new_kmer = apply_all_snps(new_kmer, {s: list_of_snps[s] for s in list_of_snps if s in curr_snps})
    new_kmer = apply_all_indels(new_kmer, {s: list_of_insertions[s] for s in curr_snps if s in list_of_insertions}, {s: list_of_deletions[s] for s in curr_snps if s in list_of_deletions})
    kmers.append(new_kmer)
  return pset, kmers         

def pretty_print(arr):
  if len(arr) == 0:
    return "()"
  ret = "("
  for element in arr:
    ret += str(element) + ","
  return ret[:-1] + ")"
def apply_all_snps(ref, snpDict):
  to_ret = ref
  for s in snpDict:
    to_ret = to_ret[:s] + snpDict[s] + to_ret[(s+1):]
  return to_ret

def apply_all_insertions(ref, insDict):
  to_ret = ref
  transformed_coords = {i: i for i in insDict}
  for i in insDict:
    ind = transformed_coords[i]
    len_ins = len(insDict[i][1]) - len(insDict[i][0])
    to_ret = to_ret[:ind+1] + insDict[i][1][1:] + to_ret[ind+1:]
    for j in insDict:
      if transformed_coords[j] > ind:
        transformed_coords[j] += len_ins
  return to_ret

# Apply all deletions to a string by padding with -, then removing all -
def apply_all_deletions(ref, delDict):
  to_ret = ref
  for d in delDict:
    len_del = len(delDict[d][0]) - len(delDict[d][1])
    to_ret = to_ret[:d+1] + '-' * len_del + to_ret[d+len_del+1:]
  return to_ret.replace('-','')

def apply_all_indels(ref, insDict, delDict):
  to_ret = ref
  transformed_coords_ins = {i: i for i in insDict}
  transformed_coords_del = {d: d for d in delDict}
  for i in insDict:
    ind = transformed_coords_ins[i]
    len_ins = len(insDict[i]) - 1
    to_ret = to_ret[:ind+1] + insDict[i][1:] + to_ret[ind+1:]
    for j in insDict:
      if transformed_coords_ins[j] > ind:
        transformed_coords_ins[j] += len_ins
    for j in delDict:
      if transformed_coords_del[j] > ind:
        transformed_coords_del[j] += len_ins
  for d in delDict:
    len_del = abs(delDict[d])
    ind = transformed_coords_del[d]
    to_ret = to_ret[:ind+1] + '-' * len_del + to_ret[ind+len_del+1:]
  return to_ret.replace('-','')

def apply_all_delins(ref, insDict, delDict):
  to_ret = ref
  transformed_coords_ins = {i: i for i in insDict}
  transformed_coords_del = {d: d for d in delDict}
  skip = set([])
  skip_i = set([])
  for d in list(delDict.keys()):
    len_del = abs(delDict[d])
    ind = transformed_coords_del[d]
    to_ret = to_ret[:ind+1] + '-' * len_del + to_ret[ind+len_del+1:]
    # Remove all variants that are contained within this deletion
    # SNPs will be deleted anyways, so only indels that extend out of the range matter
    for j in list(insDict.keys()):
      if transformed_coords_ins[j] > ind:
        if transformed_coords_ins[j] > ind + len_del:
          transformed_coords_ins[j] -= len_del
        else:
          # del transformed_coords_ins[j]
          # del insDict[j]
          skip_i.add(j)
    for j in delDict:
      if transformed_coords_del[j] > ind:
        if transformed_coords_del[j] > ind + len_del:
          transformed_coords_del[j] -= len_del
        else:
          skip.add(j)
      elif transformed_coords_del[j] < ind and transformed_coords_del[j] + delDict[j] > ind:
        skip.add(j)
    to_ret = to_ret.replace('-','')

  for i in insDict:
    ind = transformed_coords_ins[i]
    len_ins = len(insDict[i]) - 1
    to_ret = to_ret[:ind+1] + insDict[i][1:] + to_ret[ind+1:]
    for j in insDict:
      if transformed_coords_ins[j] > ind:
        transformed_coords_ins[j] += len_ins
  return to_ret

def shift(ind, shift_amt, snps, ins, dels):
  ret_snps = {}
  ret_ins = {}
  ret_dels = {}
  for s in snps:
    if s > ind:
      ret_snps[s - shift_amt] = snps[s]
  for i in ins:
    if i > ind:
      ret_ins[i - shift_amt] = ins[i]
  for d in dels:
    if d > ind:
      ret_dels[d - shift_amt] = dels[d]
  return ret_snps, ret_ins, ret_dels

# Pre: Segregated variants, applying some of the variants should never overlap with others
# These variants actually can overlap with others. What do we do about it? Delete out?
# Post: Translate all the variants that aren't being applied. Impute in homo alt
def apply_some_delins(
    ref, insDict, delDict, \
    otherLeftSnps={}, otherLeftIns={}, otherLeftDels={}, \
    otherRightSnps={}, otherRightIns={}, otherRightDels={}):
  to_ret = ref
  transformed_coords_ins = {i: i for i in insDict}
  transformed_coords_del = {d: d for d in delDict}
  for d in list(delDict.keys()):
    len_del = abs(delDict[d])
    ind = transformed_coords_del[d]
    to_ret = to_ret[:ind+1] + '-' * len_del + to_ret[ind+len_del+1:]
    # Remove all variants that are contained within this deletion
    # SNPs will be deleted anyways, so only indels that extend out of the range matter
    for j in list(insDict.keys()):
      if transformed_coords_ins[j] > ind:
        if transformed_coords_ins[j] > ind + len_del:
          transformed_coords_ins[j] -= len_del
    for j in delDict:
      if transformed_coords_del[j] > ind:
        if transformed_coords_del[j] > ind + len_del:
          transformed_coords_del[j] -= len_del

    otherLeftSnps, otherLeftIns, otherLeftDels = \
        shift(ind, len_del, otherLeftSnps, otherLeftIns, otherLeftDels)
    otherRightSnps, otherRightIns, otherRightDels = \
        shift(ind, len_del, otherRightSnps, otherRightIns, otherRightDels)
    to_ret = to_ret.replace('-','')

  for i in insDict:
    ind = transformed_coords_ins[i]
    len_ins = len(insDict[i]) - 1
    to_ret = to_ret[:ind+1] + insDict[i][1:] + to_ret[ind+1:]
    for j in insDict:
      if transformed_coords_ins[j] > ind:
        transformed_coords_ins[j] += len_ins
    otherLeftSnps, otherLeftIns, otherLeftDels = \
        shift(ind, -1 * len_ins, otherLeftSnps, otherLeftIns, otherLeftDels)
    otherRightSnps, otherRightIns, otherRightDels = \
        shift(ind, -1 * len_ins, otherRightSnps, otherRightIns, otherRightDels)

  return to_ret, otherLeftSnps, otherLeftIns, otherLeftDels, otherRightSnps, otherRightIns, otherRightDels 


def apply_edges(ref, edges):
  assert(len(edges) == 1)
  k = list(edges.keys())[0]
  assert(k <= 0)
  del_num = k + edges[k]
  if del_num >= 0:
    return ref[del_num+1:], k, edges[k]
  else:
    return ref, k, edges[k]

def get_edge_info(edges):
  assert(len(edges) == 1)
  k = list(edges.keys())[0]
  assert(k <= 0)
  return k, edges[k]

# Load in all the SNPs for the transcriptome
# This function is for ornament transcriptomes
def load_transcriptome_snps_diploid(input_tvcf, sample):
  transcriptome_snps = {}
  num_skipped = 0
  num_processed = 0
  # load in the left and the right information for the sample that we are interested in
  with open(input_tvcf) as tvcf:
    start = False
    ind = []
    for line in tvcf:
      if '#' in line:
        if 'CHROM' in line:
          toks = [tok.strip() for tok in line.split('\t')]
          ind = toks.index(sample)
          print('Using ind ', ind, ' for sample ', sample)
        continue
    
      tokens = line.split('\t')
      tid, pos, _, ref, alt = tokens[:5]
      pos = int(pos)
      num_ppl = 0
      for genotype in tokens[9:]:
        if '1' in genotype:
          num_ppl += 1
        
      genotype = tokens[ind].strip()
      if len(genotype) == 3 and (genotype == '0|0' or genotype == '0/0'):
          continue
      elif len(genotype) == 1 and genotype == '0':
          continue
      if len(genotype) == 3:
        left = ref if genotype[0] == '0' else alt
        right = ref if genotype[2] == '0' else alt
      else:
        assert len(genotype) == 1
        left = ref if genotype[0] == '0' else alt
        right = ref if genotype[0] == '0' else alt

      if tid not in transcriptome_snps:
        transcriptome_snps[tid] = {}
      
      if pos - 1 not in transcriptome_snps[tid]:
        transcriptome_snps[tid][pos - 1] = []

      transcriptome_snps[tid][pos - 1].append(
          (ref, alt, left, right, len(genotype), len(ref) == len(alt) or alt.isdigit()))
  return transcriptome_snps, num_processed, num_skipped

# Load in all the SNPs for the transcriptome
# This function is for ornament transcriptomes
def load_transcriptome_snps_ornaments(input_tvcf, samples=None):
  transcriptome_snps = {}
  num_skipped = 0
  num_processed = 0
  # load in the left and the right information for the sample that we are interested in
  inds = []
  with open(input_tvcf) as tvcf:
    start = False
    for line in tvcf:
      if '#' in line:
        toks = [tok.strip() for tok in line.split('\t')]
        if samples is None:
          inds = list(range(9, len(toks)))
          print(inds)
        else:
          inds = [toks.index(sample) for sample in samples]
        print('Using inds ', inds, ' for samples ', samples)
        continue
    
      tokens = line.split('\t')
      tid, pos, _, ref, alt = tokens[:5]
      
      # if len(ref) != len(alt): # no longer skipping indels!
      #   num_skipped += 1
      #   continue
      # else:
      #   num_processed += 1
      
      pos = int(pos)
      
      num_ppl = 0
      for genotype in tokens[9:]:
        if '1' in genotype:
          num_ppl += 1
        
      gts = [tokens[ind].strip() for ind in inds]
      is_homo_alt = True
      is_homo_ref = True
      for gt in gts:
        if gt != '1|1' and gt != '1' and gt != '1/1':
          is_homo_alt = False
        if gt != '0|0' and gt != '0' and gt != '0/0':
          is_homo_ref = False
      
      if tid not in transcriptome_snps:
        transcriptome_snps[tid] = {}
      
      # THESE ARE ZERO INDEXED
      transcriptome_snps[tid][pos-1] = (ref, alt, is_homo_alt, is_homo_ref, "\t".join(gts), len(ref) == len(alt))
  return transcriptome_snps, num_processed, num_skipped, inds

def load_transcriptome_snps_phased_ornament(input_tvcf, samples=None):
  transcriptome_snps = {}
  num_skipped = 0
  num_processed = 0
  num_snps, num_insertions, num_deletions = 0, 0, 0
  # load in the left and the right information for the sample that we are interested in
  with open(input_tvcf) as tvcf:
    start = False
    ind = []
    for line in tvcf:
      if '#' in line:
        if 'CHROM' in line:
          toks = [tok.strip() for tok in line.split('\t')]
          inds = [toks.index(sample) for sample in samples]
          print('Using ind ', inds, ' for samples ', samples)
        continue
    
      tokens = line.split('\t')
      tid, pos, _, ref, alt = tokens[:5]
      pos = int(pos)
      num_ppl = 0
      for genotype in tokens[9:]:
        if '1' in genotype:
          num_ppl += 1
        
      gts = [tokens[ind].strip() for ind in inds]
      is_homo_alt = True
      is_homo_ref = True
      for gt in gts:
        if gt != '1|1' and gt != '1' and gt != '1/1':
          is_homo_alt = False
        if gt != '0|0' and gt != '0' and gt != '0/0':
          is_homo_ref = False

      if is_homo_ref:
        continue
      if tid not in transcriptome_snps:
        transcriptome_snps[tid] = TranscriptVariantList(True, tid)
      
      # If this is present in at least one of the samples
      transcriptome_snps[tid].add(
          pos - 1,
          ref,
          alt,
          is_homo_alt,
          is_homo_ref,
          gts,
          len(ref) == len(alt) or alt.isdigit())

      if len(ref) == len(alt):
          num_snps += 1
      elif len(ref) > len(alt):
          num_deletions += 1
      elif len(ref) < len(alt):
          num_insertions += 1

    print(f'Using {num_snps} heterozygous SNP(s), {num_insertions} heterozygous insertion(s), and {num_deletions} heterozygous deletion(s)')
    return transcriptome_snps, num_processed, num_skipped

class Printer:
  def __init__(self, f, c=None):
    self.f = open(f, 'w+')
    if c:
      self.c = open(c, 'w+')
    else:
      self.c = None

  def write_sequence(self, tname, seq):
    self.f.write('>' + tname + '\n')
    for i in range((len(seq) + 59) // 60):
      self.f.write(seq[i*60:(i+1)*60] + '\n')

  def write_coordinate(self, tname, haplotype, t, v, opos, npos):
    if self.c:
      self.c.write('\t'.join([tname, haplotype, t, str(v), str(opos), str(npos)]) + '\n')

  def close(self):
    self.f.close()
    if self.c:
      self.c.close()

class VariantToVCF:
  def __init__(self):
    self.snps = {}
    self.target_to_snp = {}
    self.ins = {}
    self.target_to_ins = {}
    self.dels = {}
    self.target_to_dels = {}
    self.edges = {}
    self.target_to_edges= {}

  def add_snp(self, pos, val, target):
    if pos not in self.snps:
      self.snps[pos] = {}
      self.target_to_snp[pos] = {}
    self.snps[pos][val] = target
    self.target_to_snp[pos][target] = val

  def get_snp(self, pos, val):
    return self.snps[pos][val]

  def get_snp_target(self, pos, target):
    if target not in self.target_to_snp[pos]:
      return self.target_to_snp[pos]['a_0']
    else:
      return self.target_to_snp[pos][target]

  def add_ins(self, pos, val, target):
    if pos not in self.ins:
      self.ins[pos] = {}
      self.target_to_ins[pos] = {}
    self.ins[pos][val] = target
    self.target_to_ins[pos][target] = val

  def get_ins(self, pos, val):
    return self.ins[pos][val]

  def get_ins_target(self, pos, target):
    if target not in self.target_to_ins[pos]:
      return self.target_to_ins[pos]['a_0']
    else:
      return self.target_to_ins[pos][target]

  def add_del(self, pos, val, target):
    if pos not in self.dels:
      self.dels[pos] = {}
      self.target_to_dels[pos] = {}
    self.dels[pos][val] = target
    self.target_to_dels[pos][target] = val

  def get_del(self, pos, val):
    return self.dels[pos][val]

  def get_del_target(self, pos, target):
    if target not in self.target_to_dels[pos]:
      return self.target_to_dels[pos]['a_0']
    else:
      return self.target_to_dels[pos][target]

class Haplotype(Enum):
  left = 0
  right = 1

class OrnamentVariant:
  def __init__(self, ref, alt, is_homo_alt, is_homo_ref, gts, is_snp):
    self.ref = ref
    self.alt = alt
    self.is_homo_alt = is_homo_alt
    self.is_homo_ref = is_homo_ref
    self.gt_string = "\t".join(gts)
    self.indel = not is_snp

  def is_snp(self):
    return not self.indel

  def is_indel(self):
    return self.indel

  def get_genotypes(self):
    return self.gt_string.split('\t')

  def get_genotype(self, ind):
    return self.gt_string.split('\t')[ind]

  def get_phased_alleles(self, ind):
    gt = self.get_genotype(ind)
    if '|' in gt:
      l, r = gt.split('|')
    elif '/' in gt:
      l, r = gt.split('/')
    else:
      assert(len(gt) == 1)
      l, r = gt, gt
    l = self.ref if l == '0' else self.alt
    r = self.ref if r == '0' else self.alt
    return l, r

  def __repr__(self):
      return f'{self.ref}, {self.alt}, {self.is_homo_alt}, {self.gt_string}, {self.indel})'

class TranscriptVariantList:
  def __init__(self, phased, tname):
    self.variants = {}
    self.phased = phased
    self.tname = tname

  def add(self, pos, ref, alt, is_homo_alt, is_homo_ref, gts, is_indel):
    if pos not in self.variants:
      self.variants[pos] = []

    self.variants[pos].append(
      OrnamentVariant(ref, alt, is_homo_alt, is_homo_ref, gts, is_indel))

  # For the phased case: Take the final set of snps before writing and applying 
  # And make ornaments out of them? Same logic as before, but impute in homo alts

  def get_homo_alt_variants(self):
    homo_alt_edges = {}
    homo_alt_snps = {}
    homo_alt_ins = {}
    homo_alt_dels = {}

    for pos in self.variants:
      for var in self.variants[pos]:
        # If the current variant is a SNP, or an edge deletion
        if len(var.ref) == len(var.alt) and var.is_snp():
          if is_homo_alt:
            if pos < 0:
              homo_alt_edges[pos] = int(var.alt)
            else:
              homo_alt_snps[pos] = var.alt
        else:
          # If the current variant is a deletion
          if is_homo_alt:
            if len(var.ref) > len(var.alt):
              homo_alt_dels[pos] = len(var.ref) - len(var.alt)
            # If the current variant is aninsertion
            else:
              homo_alt_ins[pos] = var.alt

    return homo_alt_edges, homo_alt_snps, homo_alt_ins, homo_alt_dels

  def get_het_variants(self, sample, ind):
    het_edges = {}
    het_snps = {}
    het_ins = {}
    het_dels = {}

    for pos in self.variants:
      for var in self.variants[pos]:
        if not Utils.is_het(var.get_genotype(ind)):
          continue

        # Need to make this work with multiallelic sites?
        # Need to make actual ornaments script work with multiallelic sites?
        # Include the allelic value of snp in the ornament?
        if len(var.ref) == len(var.alt) and var.is_snp():
          if pos < 0:
            het_edges[pos] = int(var.alt)
          else:
            het_snps[pos] = var.alt

  def what_am_i(self, pos, val):
    a_num = 0
    for var in self.variants[pos]:
      if var.is_snp():
        if pos < 0:
          if val == '.':
            return 'r'
          else:
            if val == pos + int(var.alt):
              return 'a_' + str(a_num)
        else:
          if val == var.ref:
            return 'r'
          elif val == var.alt:
            return 'a_' + str(a_num)
      else:
        if len(var.ref) > len(var.alt):
          if val == 0:
            return 'r'
          else:
            if val == len(var.ref) - len(var.alt):
              return 'a_' + str(a_num)
        else:
          if val == var.ref:
            return 'r'
          elif val == var.alt:
            return 'a_' + str(a_num)
      a_num += 1
      
  def get_phased_ornament_variants(self, ind, hap):
    snps = {}
    ins = {}
    dels = {}
    edges = {}
    homo_alt_snps = set()
    homo_alt_ins = set()
    homo_alt_dels = set()
    homo_alt_edges = set()
    # This is used to map the variants and values back to what they are in the vcf file
    indexer = VariantToVCF()

    for pos in self.variants:
      a_num = 0
      for var in self.variants[pos]:
        ref_allele = var.ref
        alt_allele = var.alt
        l, r = var.get_phased_alleles(ind)
        al = l if hap == Haplotype.left else r
        if var.is_snp():
          if pos < 0:
            # Variant is present and homozygous
            # If pos isn't there yet, then add it to the list of variants at this position
            # This is phased -- only add one per sample! Overwrite with alt if it ever shows up
            if al != '.':
              # if pos not in edges:
              #   edges[pos] = []
              if al == ref_allele:
                if pos not in edges:
                  edges[pos] = [ref_allele]
                # if ref_allele not in edges[pos]:
                  # edges[pos].append(ref_allele)
              else:
                assert(al == alt_allele)
                # TODO: Inconsistent vcf file will make this cry 
                edges[pos] = [int(alt_allele)]
                # if alt_allele not in edges[pos]:
                #   edges[pos].append(int(alt_allele))
              if var.is_homo_alt:
                assert(len(edges[pos]) == 1)
                homo_alt_edges.add(pos)
          else:
            indexer.add_snp(pos, ref_allele, 'r')
            indexer.add_snp(pos, alt_allele, 'a_' + str(a_num))

            if al == ref_allele:
              if pos not in snps:
                snps[pos] = [ref_allele]
            else:
              assert(al == alt_allele)
              snps[pos] = [alt_allele]
            if var.is_homo_alt:
              assert(len(snps[pos]) == 1)
              homo_alt_snps.add(pos)
        else:
          # If this is a deletion
          if len(ref_allele) > len(alt_allele):
            # If the deletion was on the left haplotype, add as deletion to left_indels
            # if pos not in right_ins and pos not in right_dels and r == alt_allele:
            indexer.add_del(pos, 0, 'r')
            indexer.add_del(pos, len(ref_allele) - len(alt_allele), 'a_' + str(a_num))
            if al == ref_allele:
              if pos not in dels:
                dels[pos] = [0]
            else:
              assert(al == alt_allele)
              dels[pos] = [len(ref_allele) - len(alt_allele)]
            if var.is_homo_alt:
              assert(len(dels[pos]) == 1)
              homo_alt_dels.add(pos)
          # If this is an insertion
          else:
            indexer.add_ins(pos, ref_allele, 'r')
            indexer.add_ins(pos, alt_allele, 'a_' + str(a_num))
            if al == ref_allele:
              if pos not in ins:
                ins[pos] = [ref_allele]
            else:
              assert(al == alt_allele)
              ins[pos] = [alt_allele]
            if var.is_homo_alt:
              assert(len(ins[pos]) == 1)
              homo_alt_ins.add(pos)

        a_num += 1

    return snps, ins, dels, edges, homo_alt_snps, homo_alt_ins, homo_alt_dels, homo_alt_edges, indexer

  # For each variant type, for each position, keep track of all alleles 
  def get_unphased_ornament_variants(self, ind):
    snps = {}
    ins = {}
    dels = {}
    edges = {}
    homo_alt_snps = set()
    homo_alt_ins = set()
    homo_alt_dels = set()
    homo_alt_edges = set()
    # This is used to map the variants and values back to what they are in the vcf file
    indexer = VariantToVCF()

    for pos in self.variants:
      a_num = 0
      for var in self.variants[pos]:
        ref_allele = var.ref
        alt_allele = var.alt
        l, r = var.get_phased_alleles(ind)
        if var.is_snp():
          if pos < 0:
            # Variant is present and homozygous
            # If pos isn't there yet, then add it to the list of variants at this position
            if l != '.' or r != '.':
              if pos not in edges:
                edges[pos] = []
              if l == ref_allele or r == ref_allele:
                if ref_allele not in edges[pos]:
                  edges[pos].append(ref_allele)
              if l == alt_allele or r == alt_allele:
                if alt_allele not in edges[pos]:
                  edges[pos].append(int(alt_allele))
              if var.is_homo_alt:
                assert(len(edges[pos]) == 1)
                homo_alt_edges.add(pos)
          else:
            indexer.add_snp(pos, ref_allele, 'r')
            indexer.add_snp(pos, alt_allele, 'a_' + str(a_num))
            if pos not in snps:
              snps[pos] = []
            if l == ref_allele or r == ref_allele:
              if ref_allele not in snps[pos]:
                snps[pos].append(ref_allele)
            if l == alt_allele or r == alt_allele:
              if alt_allele not in snps[pos]:
                snps[pos].append(alt_allele)
            if var.is_homo_alt:
              assert(len(snps[pos]) == 1)
              homo_alt_snps.add(pos)
        else:
          # If this is a deletion
          if len(ref_allele) > len(alt_allele):
            # If the deletion was on the left haplotype, add as deletion to left_indels
            # if pos not in right_ins and pos not in right_dels and r == alt_allele:
            indexer.add_del(pos, 0, 'r')
            indexer.add_del(pos, len(ref_allele) - len(alt_allele), 'a_' + str(a_num))
            if pos not in dels:
              dels[pos] = []
            if l == ref_allele or r == ref_allele:
              dels[pos].append(0)
            if l == alt_allele or r == alt_allele:
              dels[pos].append(len(ref_allele) - len(alt_allele))
            if var.is_homo_alt:
              # assert(len(dels[pos]) == 1)
              if len(dels[pos]) != 1:
                dels[pos] = dels[pos][-1:]
              homo_alt_dels.add(pos)
          # If this is an insertion
          else:
            indexer.add_ins(pos, ref_allele, 'r')
            indexer.add_ins(pos, alt_allele, 'a_' + str(a_num))
            if pos not in ins:
              ins[pos] = []
            if l == ref_allele or r == ref_allele:
              ins[pos].append(ref_allele)
            if l == alt_allele or r == alt_allele:
              ins[pos].append(alt_allele)
            if var.is_homo_alt:
              assert(len(ins[pos]) == 1)
              homo_alt_ins.add(pos)

        a_num += 1

    return snps, ins, dels, edges, homo_alt_snps, homo_alt_ins, homo_alt_dels, homo_alt_edges, indexer

  # If 0|1 and then 1|1, then variant is skipped on the right haplotype but not left
  # So, only have one ornament in this case? But what if there are multiple samples?
  def get_phased_variants(self, ind):
    left_snps = {}
    right_snps = {}
    left_dels = {}
    right_dels = {}
    left_ins = {}
    right_ins = {}
    left_edges = {}
    right_edges = {}
    for pos in self.variants:
      for var in self.variants[pos]:
        ref_allele = var.ref
        alt_allele = var.alt
        l, r = var.get_phased_alleles(ind)
        # print(l, r, var.is_snp())
        if var.is_snp():
          if pos < 0:
            if l != '.':
              left_edges[pos] = int(l)
            if r != '.':
              right_edges[pos] = int(r)
          else:
            if pos not in left_snps or l == alt_allele:
              left_snps[pos] = l
            if pos not in right_snps or r == alt_allele:
              right_snps[pos] = r
        else:
          # If this is a deletion
          if len(ref_allele) > len(alt_allele):
            # If the deletion was on the left haplotype, add as deletion to left_indels
            if pos not in right_ins and pos not in right_dels and r == alt_allele:
              right_dels[pos] = len(ref_allele) - len(alt_allele)
            if pos not in left_ins and pos not in left_dels and l == alt_allele:
              left_dels[pos] = len(ref_allele) - len(alt_allele)
          # If this is an insertion
          else:
            if pos not in right_ins and pos not in right_dels and r == alt_allele:
              right_ins[pos] = r
            if pos not in left_ins and pos not in left_dels and l == alt_allele:
              left_ins[pos] = l

    return left_edges, left_snps, left_ins, left_dels, right_edges, right_snps, right_ins, right_dels

  # Apply the homozygous alternative snps, and filter the remaining variants and their coordinates accordingly
  def apply_homo_alts(self, sequence, left_edges, left_snps, left_ins, left_dels, right_edges, right_snps, right_ins, right_dels):
    homo_alt_snps = {}
    homo_alt_ins = {}
    homo_alt_dels = {}
    homo_alt_edges = {}
    for s in left_snps:
      if s in right_snps:
        # Homo alt
        if left_snps[s] == right_snps[s]:
          homo_alt_snps[s] = left_snps[s]
        # Multiallelic site
        # else:
        #   do some other shit
    for s in homo_alt_snps:
      del left_snps[s]
      del right_snps[s]

    for d in left_dels:
      if d in right_dels:
        if left_dels[d] == right_dels[d]:
          homo_alt_dels[d] = left_dels[d]
    for d in homo_alt_dels:
      del left_dels[d]
      del right_dels[d]

    for i in left_ins:
      if i in right_ins:
        if left_ins[i] == right_ins[i]:
          homo_alt_ins[i] = left_ins[i]
    for i in homo_alt_ins:
      del left_ins[i]
      del right_ins[i]

    for e in left_edges:
      if e in right_edges:
        if left_edges[e] == right_edges[e]:
          homo_alt_edges[e] = left_edges[e]
    for e in homo_alt_edges:
      del left_edges[e]
      del right_edges[e]

def write_phased_ornament_transcriptome(
    transcript_to_sequence,
    transcriptome_snps,
    output_file,
    samples,
    k=31
):

  printer = Printer(output_file)
  for tname in transcript_to_sequence:
    sequence = transcript_to_sequence[tname]
    tname = tname.split('|')[0]

    if tname not in transcriptome_snps:
      printer.write_sequence(tname, sequence)
      continue

    # 1. tsnps is a list of variants. filter for set of het and homo alt for this sample
    write_base_seq = False
    encountered = {}
    for sind in range(len(samples)):
      for hap in [Haplotype.left, Haplotype.right]:
        tsnps = transcriptome_snps[tname]
        snps, ins, dels, edges, homo_alt_snps, homo_alt_ins, homo_alt_dels, homo_alt_edges, indexer = \
            tsnps.get_phased_ornament_variants(sind, hap)
        homo_alt_imputed_sequence, _, _, _, _, _, _, _ = apply_dem_snps_and_indels_and_edges(
            sequence,
            {s: snps[s][0] for s in snps if s in homo_alt_snps},
            {},
            {},
            {})

        if not write_base_seq:
          printer.write_sequence(tname, homo_alt_imputed_sequence)
          write_base_seq = True

        non_homoa_snps = {s: snps[s] for s in snps if s not in homo_alt_snps}
        non_homoa_ins = {i: ins[i] for i in ins if i not in homo_alt_ins}
        non_homoa_dels = {d: dels[d] for d in dels if d not in homo_alt_dels}

        homoa_snps = {s: snps[s] for s in snps if s in homo_alt_snps}
        homoa_ins = {i: ins[i] for i in ins if i in homo_alt_ins}
        homoa_dels = {d: dels[d] for d in dels if d in homo_alt_dels}

        snp_partitions, extract_regions = calculate_extract_for_transcript(
          sorted(non_homoa_snps.keys()), ins, dels, k=k
        )
        for csnps, erange in zip(snp_partitions, extract_regions):
          if erange[0] < 0:
            lbound = -1
          else:
            lbound = erange[0]

          if erange[1] > len(sequence):
            ubound = len(sequence)
          else:
            ubound = erange[1]

          lbound += 1
          current_snps = {s - lbound: snps[s] for s in csnps if s in snps}
          current_ins = {i - lbound: ins[i] for i in csnps if i in ins}
          current_dels = {d - lbound: dels[d] for d in csnps if d in dels}
          roi = sequence[lbound:ubound]

          ranges, snps_for_range, indels_for_range, subsequences = \
              calculate_subsequences_new(roi, list(current_snps.keys()), current_ins, current_dels, k=k)
          # split subindels into subins and subdels?
          for subrange, subsnps, subindels, subseq in zip(
              ranges, snps_for_range, indels_for_range, subsequences):

            # print(tname, lbound, subrange, subsnps, subindels, subseq)
            homoa_csnps, homoa_cins, homoa_cdels = contained_snps_and_indels(
                homoa_snps, {}, {}, # homoa_ins, homoa_dels,
                (subrange[0] + lbound, subrange[1] + lbound),
                lbound + subrange[0])
            # TODO: Need to subtract subrange bound from all of these and further filter.

            cc_snps = {s - subrange[0]: current_snps[s] 
                for s in subsnps}
            cc_ins = {i - subrange[0]: current_ins[i] 
                for i in subindels if type(subindels[i][0]) is not int}
            cc_dels = {d - subrange[0]: current_dels[d] 
                for d in subindels if type(subindels[d][0]) is int}

            for orn_name, orn_seq in apply_snps_indels(
                subseq, k, indexer, 
                {**cc_snps, **homoa_csnps}, {**cc_ins, **homoa_cins}, {**cc_dels, **homoa_cdels}, {},
                set(homoa_csnps.keys()), set(homoa_cins.keys()), set(homoa_cdels.keys()), set(),
                lbound + subrange[0], phased=True):

              if orn_name not in encountered:
                encountered[orn_name] = orn_seq
              else:
                if len(orn_seq) > len(encountered[orn_name]):
                  encountered[orn_name] = orn_seq

    for ok in encountered:
      if len(encountered[ok]) >= k:
        printer.write_sequence(tname + ok, encountered[ok])

# Write for a single sample first
def write_unphased_ornament_transcriptome(
    transcript_to_sequence,
    transcriptome_snps,
    output_file,
    samples,
    k=31
):

  # 1. Get a list of all heterozygous and homozygous alt variants for a sample
  # 2. Apply all the homozygous alt variants for a sample, and modify the coordinates of het loci
  # 3. Given the list of heterozygous variants for a sample, calculate grouped snps
  # 4. Calculate all possible combinations of these grouped snps. Retain overlapping variants
  #    in the same way that the diploid transcriptome retains variants
  # 5. Calculate the shaded k-mers, by applying snps and then delins, and then trimming the resulting
  #    sequences to make sure that k-mers are not shared across the ornaments
  printer = Printer(output_file)
  for tname in transcript_to_sequence:
    sequence = transcript_to_sequence[tname]
    tname = tname.split('|')[0]

    if tname not in transcriptome_snps:
      printer.write_sequence(tname, sequence)
      continue

    # print(tname)
    # 1. tsnps is a list of variants. filter for set of het and homo alt for this sample
    tsnps = transcriptome_snps[tname]
    snps, ins, dels, edges, homo_alt_snps, homo_alt_ins, homo_alt_dels, homo_alt_edges, indexer = \
        tsnps.get_unphased_ornament_variants(0)
    homo_alt_imputed_sequence, _, _, _, _, _, _, _ = apply_dem_snps_and_indels_and_edges(
        sequence,
        {s: snps[s][0] for s in snps if s in homo_alt_snps},
        {},
        {},
        {})

    printer.write_sequence(tname, homo_alt_imputed_sequence)

    encountered = {}
    # If a variant from left and a variant from right would overlap, don't include them in ornament
    # So, merge together indices. Break into combos, only make compatible variants within each range

    # Need to be able to label the variants
    # TODO: Make sure that we can make a phased transcriptome too, just skip over sites
    # that would have wrong phase
    non_homoa_snps = {s: snps[s] for s in snps if s not in homo_alt_snps}
    non_homoa_ins = {i: ins[i] for i in ins if i not in homo_alt_ins}
    non_homoa_dels = {d: dels[d] for d in dels if d not in homo_alt_dels}

    homoa_snps = {s: snps[s] for s in snps if s in homo_alt_snps}
    homoa_ins = {i: ins[i] for i in ins if i in homo_alt_ins}
    homoa_dels = {d: dels[d] for d in dels if d in homo_alt_dels}

    # Don't use homo alt for calculating extract
    # But, impute in homo alt for all ornaments around them, and impute into the base sequence
    # say we have a homo alt snp that would have been eaten by a homo alt edge
    # do we impute in homo alt snps? what if they would be deleted by a homo alt deletion
    snp_partitions, extract_regions = calculate_extract_for_transcript(
      sorted(non_homoa_snps.keys()), ins, dels, k=k
    )
    for csnps, erange in zip(snp_partitions, extract_regions):
      if erange[0] < 0:
        lbound = -1
      else:
        lbound = erange[0]

      if erange[1] > len(sequence):
        ubound = len(sequence)
      else:
        ubound = erange[1]

      lbound += 1
      # Need to also get homo alt snps, homo alt indels, homo alt dels in this region
      # Subsequence represents the range for these variants.
      # Need to apply them in a phase appropriate manner
      # Need to include edges here as well
      current_snps = {s - lbound: snps[s] for s in csnps if s in snps}
      current_ins = {i - lbound: ins[i] for i in csnps if i in ins}
      current_dels = {d - lbound: dels[d] for d in csnps if d in dels}
      roi = sequence[lbound:ubound]

      # Need to calculate subsequences without the homo alt snps
      # Then, add back in the ones that are in the current range to use for the apply_snps calculation
      # But, homo_alts should not be considered at all. 
      # If homo_alt is squelched, then _________
      # If anything else is squelched, then don't mention it at all
      # Different from non-membership, which is when something is ref_snp

      ranges, snps_for_range, indels_for_range, subsequences = \
          calculate_subsequences_new(roi, list(current_snps.keys()), current_ins, current_dels, k=k)
      # split subindels into subins and subdels?
      for subrange, subsnps, subindels, subseq in zip(
          ranges, snps_for_range, indels_for_range, subsequences):

        # print(tname, lbound, subrange, subsnps, subindels, subseq)
        homoa_csnps, _, _ = contained_snps_and_indels(homoa_snps, {}, {}, (subrange[0] + lbound, subrange[1] + lbound), lbound + subrange[0])
        # TODO: Need to subtract subrange bound from all of these and further filter.

        cc_snps = {s - subrange[0]: current_snps[s] for s in subsnps}
        cc_ins = {i - subrange[0]: current_ins[i] for i in subindels if type(subindels[i][0]) is not int}
        cc_dels = {d - subrange[0]: current_dels[d] for d in subindels if type(subindels[d][0]) is int}

        for orn_name, orn_seq in apply_snps_indels(
            subseq, k, indexer, 
            {**cc_snps, **homoa_csnps}, cc_ins, cc_dels, {}, 
            set(homoa_csnps.keys()), set(), set(), set(), lbound + subrange[0]):
          if orn_name not in encountered:
            encountered[orn_name] = orn_seq
          else:
            if len(orn_seq) > len(encountered[orn_name]):
              encountered[orn_name] = orn_seq

    for ok in encountered:
      if len(encountered[ok]) >= k:
        printer.write_sequence(tname + ok, encountered[ok])

        # Simplest way: undo that change, add trimmer in phased case
def write_diploid_transcriptome(
    transcript_to_sequence, transcriptome_snps, output_file, coordinate_file=None):
  printer = Printer(output_file, coordinate_file)
  for tname in transcript_to_sequence:
    ref = transcript_to_sequence[tname]
    short_tname = tname.split('|')[0]
    left_ref = ref
    right_ref = ref
    lsnpmap, linsmap, ldelmap = {}, {}, {}
    rsnpmap, rinsmap, rdelmap = {}, {}, {}

    if short_tname in transcriptome_snps: # this is a variant transcript, or homozygous alt
      tsnps = transcriptome_snps[short_tname]
      left_edges, left_snps, left_ins, left_dels, right_edges, right_snps, right_ins, right_dels = \
          tsnps.get_phased_variants(0)
          
      # sequence, coordinate maps, sets of squelched variants
      left_ref, lsnpmap, linsmap, ldelmap, lss, lsi, lsd, lse = apply_dem_snps_and_indels_and_edges(
          left_ref,
          left_snps,
          left_ins,
          left_dels,
          left_edges)

      right_ref, rsnpmap, rinsmap, rdelmap, rss, rsi, rsd, rse = apply_dem_snps_and_indels_and_edges(
          right_ref,
          right_snps,
          right_ins,
          right_dels,
          right_edges)


      
    # print left version of transcript and then right version of transcript
    printer.write_sequence(short_tname + '_L', left_ref)
    printer.write_sequence(short_tname + '_R', right_ref)
    # if user wants to output coordinates, then print these out
    for s in lsnpmap:
      printer.write_coordinate(short_tname, 'L', 's', left_snps[s], s, lsnpmap[s])
    for i in linsmap:
      printer.write_coordinate(short_tname, 'L', 'i', left_ins[i], i, linsmap[i])
    for d in ldelmap:
      printer.write_coordinate(short_tname, 'L', 'd', left_dels[d], d, ldelmap[d])
    for s in rsnpmap:
      printer.write_coordinate(short_tname, 'R', 's', right_snps[s], s, rsnpmap[s])
    for i in rinsmap:
      printer.write_coordinate(short_tname, 'R', 'i', right_ins[i], i, rinsmap[i])
    for d in rdelmap:
      printer.write_coordinate(short_tname, 'R', 'd', right_dels[d], d, rdelmap[d])

  printer.close()

# TODO: For bias correction, want to do orn check against what would have been a solo orn trans
def write_ornament_transcriptome(transcript_to_sequence, transcriptome_snps, output_file_prefix, bias_correction, phased, samples, inds):
  num_processed = 0
  bias_proportions = {}
  bias_files = []
  if bias_correction:
    for ind in range(len(samples)):
      bias_files.append(open(output_file_prefix + '.' + samples[ind] + '.bias.csv', 'w+'))

  with open(output_file_prefix + '.fa', 'w+') as otrans:
    for transcript in transcript_to_sequence:
      num_processed += 1
      if num_processed % 10000 == 0:
        print(num_processed)

      name = transcript
      sequence = transcript_to_sequence[transcript]
      name = name.split('|')[0]
      # print(name)
      # if name != 'ENST00000642116.1':
      #   continue
      sample_level_sequences = []
      not_imputed_in_alts = []
      not_imputed_in_refs = []

      if name in transcriptome_snps:
        tsnps = transcriptome_snps[name]
        # get a list of all snps for this transcript that are homozygous alternative across all samples
        homo_alt_snps = {}
        other_snps = {}
        homo_alt_ins = {}
        other_ins = {}
        homo_alt_dels = {}
        other_dels = {}

        other_indels = {}
        for pos in tsnps:
          if tsnps[pos][2] and tsnps[pos][5]:
            homo_alt_snps[pos] = tsnps[pos][1]
          elif tsnps[pos][2] and not tsnps[pos][5] and len(tsnps[pos][0]) > len(tsnps[pos][1]):
            homo_alt_dels[pos] = tsnps[pos][1]
          elif tsnps[pos][2] and not tsnps[pos][5] and len(tsnps[pos][0]) < len(tsnps[pos][1]):
            homo_alt_ins[pos] = tsnps[pos][1]
          elif pos not in homo_alt_snps and not tsnps[pos][3] and tsnps[pos][5]:
            other_snps[pos] = tsnps[pos]
          elif pos not in homo_alt_snps and not tsnps[pos][3] and not tsnps[pos][5] and len(tsnps[pos][0]) > len(tsnps[pos][1]):
            other_dels[pos] = len(tsnps[pos][1]) - len(tsnps[pos][0])
            other_indels[pos] = len(tsnps[pos][1]) - len(tsnps[pos][0])
          elif pos not in homo_alt_snps and not tsnps[pos][3] and not tsnps[pos][5] and len(tsnps[pos][0]) < len(tsnps[pos][1]):
            other_ins[pos] = tsnps[pos][1]
            other_indels[pos] = tsnps[pos][1]

        sequence = apply_all_snps(sequence, homo_alt_snps)
        # To get this to be fully functional, need to edit the map so that all het snps are updated
        # after these imputations
        # TODO: this
        # sequence = apply_all_indels(sequence, homo_alt_ins, homo_alt_dels)

        if bias_correction:
          for ind in range(len(samples)):
            additional_homo_snps = {}
            not_imputed_in_alts.append(set([]))
            not_imputed_in_refs.append(set([]))
            for pos in tsnps:
              gts = tsnps[pos][4].split('\t')
              if gts[ind] == '1|1' or gts[ind] == '1':
                additional_homo_snps[pos] = tsnps[pos][1]
                not_imputed_in_alts[-1].add(pos)
              elif gts[ind] == '0|0' or gts[ind] == '0':
                additional_homo_snps[pos] = tsnps[pos][0]
                not_imputed_in_refs[-1].add(pos)
            sample_level_sequences.append(apply_all_snps(sequence, additional_homo_snps))

        
      # First, print the current transcript into the file. Impute in any homozygous alt SNPs at this stage
      # After this point, we have dealt with all SNPs that are homo ref or homo alt in the samples
      otrans.write('>' + name + '\n')
      for i in range((len(sequence) + 59) // 60):
        otrans.write(sequence[i*60:(i+1)*60] + '\n')
      
      # Calculate where the SNPs are, and how much we have to extract for each
      if name not in transcriptome_snps:
        continue
        
      # print(name)
      snp_partitions, extract_regions = calculate_extract_for_transcript(
        sorted(list(other_snps.keys())), other_indels
      )
      print('SNP PARTITIONS: ', snp_partitions)
      print('EXTRACTED REGIONS: ', extract_regions)
      for snps, erange in zip(snp_partitions, extract_regions):
        # If the range goes below 0 or after the transcript sequence length, clip it
        if erange[0] < 0:
          lbound = -1
        else:
          lbound = erange[0]
          
        # Calculate the subsequence in the relevant range, and what the reference letters and snp letters should be
        lbound += 1
        subsequence = sequence[lbound:erange[1]]
        
        # For each SNP, calcualte the breakdown of the subsequence into its constituent parts.
        # Ranges are inclusive on both ends
        adjusted_snps = [s - lbound for s in snps if s < len(sequence) and s not in other_indels]
        adjusted_indels = {p - lbound: other_indels[p] for p in other_indels if p < len(sequence) and p in snps}
        if len(adjusted_indels) > 1 or len(adjusted_indels) == 1 and len(adjusted_snps) > 0:
          adjusted_indels = {}
        if len(adjusted_snps) == 0 and len(adjusted_indels) == 0:
          continue
        
        print(snps, other_indels)
        subseq_ranges, snps_in_ranges, indels_in_ranges, unshifted_indels_ir, subseq_seq \
            = calculate_subsequences(subsequence, adjusted_snps, adjusted_indels)
        # Iterate over the snps
        for subseq_range, relevant_snps, relevant_indels, unshifted_indels, subseq in \
            zip(subseq_ranges, snps_in_ranges, indels_in_ranges, unshifted_indels_ir, subseq_seq):
          ref_letters = list(map(lambda x: other_snps[x][0], [s + lbound for s in relevant_snps]))
          snp_letters = list(map(lambda x: other_snps[x][1], [s + lbound for s in relevant_snps]))
          pos_to_snp_let = {k: v for k, v in zip([s - subseq_range[0] for s in relevant_snps], snp_letters)}
          # need to subtract out relative offset from segmenting the snps
          # Multisequences now contains the relevant strings, need to form headers for them now
          print('INSIDE LOOP', pos_to_snp_let, subseq_range, relevant_snps, relevant_indels, unshifted_indels)
          pset, multisequences = apply_snps(
            subseq,
            {k - subseq_range[0]: tsnps[k+lbound][4] for k in [*relevant_snps, *relevant_indels]},
            inds,
            phased,
            pos_to_snp_let,
            relevant_indels)
          # Currently applied snps are snps that currently have the alternative value
          currently_alt_snps = [np.array(list(p)) + subseq_range[0] + lbound for p in pset]
          if type(unshifted_indels) is not int:
            all_snps = set([s + lbound for s in relevant_snps + list(unshifted_indels.keys())])
          else:
            all_snps = set([s + lbound for s in relevant_snps])

          # Currently ref snps are snps that currently have the reference value
          currently_ref_snps = [list(all_snps - set(curr_alt)) for curr_alt in currently_alt_snps]
          
          for a, r, ms in zip(currently_alt_snps, currently_ref_snps, multisequences): 
            # if len(ms) < 31:
            #   continue
            shade_str = name + '_' + 'refsnp' + '_' + pretty_print(r) + '|' + 'altsnp' + '_' + pretty_print(a)
            otrans.write('>' + shade_str + '\n')
            # Write it in the fasta format
            for i in range((len(ms) + 59) // 60):
              otrans.write(ms[i*60:(i+1)*60] + '\n')

            # if bias correction is on, check and see if there is a bias proportion for this transcript
            if bias_correction:
              # First construct that the sequence would have been for each individual, i.e. impute in snps where that ind is 1/1
              if len(ms) < 31:
                continue

              for ind in range(len(samples)):
                curr_sample = samples[ind]
                curr_sample_sequence = sample_level_sequences[ind]

                variant_real_estate = 0
                total_real_estate = 0
                rem_rs = set(r) - not_imputed_in_refs[ind]
                rem_as = set(a) - not_imputed_in_alts[ind]
                if len(rem_rs) == 0 and len(rem_as) == 0:
                  continue
                is_all_ref = len(rem_rs) > 0 and len(rem_as) == 0
                for i in range(len(ms) - 30):
                  kmer = ms[i:i+31]
                  variant_real_estate += len(kmer)
                  if is_all_ref:
                    c = curr_sample_sequence.count(kmer)
                    total_real_estate += (2*c-1) * len(kmer)
                  else:
                    c = curr_sample_sequence.count(kmer)
                    total_real_estate += (2*c+1) * len(kmer)

                if total_real_estate == 0:
                  continue

                variant_proportion = variant_real_estate / total_real_estate
                if variant_proportion == 1:
                  continue

                bias_files[ind].write(shade_str + '\t' + str(variant_proportion) + '\n')
    

# Command line arguments: transcriptome, transcriptome vcf file, type (diploid or ornaments), output file name prefix
def main(tfasta, tvcf, output_type, output_file, samples, coordinates, k):
  if len(sys.argv) == 2 and sys.argv[1] == '-h':
    print('USAGE: transcriptome fasta file, transcriptome vcf file, diploid/ornaments, output filename prefix, sample(s), bias correction (off by default)')
    return

  bias_correction = False
  phased = False

  if output_type != 'diploid' and output_type != 'ornament':
    print('Invalid output type. Please specify either "diploid" or "ornament"')
    print(sys.argv)
    return

  transcript_to_sequence = {}
  print('Loading in transcript sequences...')
  for tname, ref in Utils.transcript_generator(tfasta):
    tname = tname.split('|')[0]
    transcript_to_sequence[tname] = ref
  print('Finished loading in sequences')
  print('Loading in transcriptome variants...')

  if output_type == 'diploid':
    if ',' in samples:
        print('Can only output one diploid transcriptome at a time')
        return
    transcriptome_snps, _, _ = load_transcriptome_snps_phased_ornament(tvcf, [samples])
    print(len(transcriptome_snps))
  else:
    if samples.strip() == 'ALL':
      current_samples = None
    else:
      print(samples)
      current_samples = [t.strip() for t in samples.split(',')]
    transcriptome_snps, _, _ = load_transcriptome_snps_phased_ornament(tvcf, current_samples)
    inds = list(range(len(current_samples)))

  print('Finished loading in variants')
  print('Writing output ' + output_type + ' transcriptome')

  if output_type == 'diploid':
    coord_file = coordinates
    if coord_file is None:
        print('Need coordinate files for diploid fasta')
        return
    write_diploid_transcriptome(transcript_to_sequence, transcriptome_snps, output_file, coord_file)
  else:
    if phased:
      write_phased_ornament_transcriptome(
              transcript_to_sequence, transcriptome_snps, output_file, current_samples, k)
    else:
      write_unphased_ornament_transcriptome(
              transcript_to_sequence, transcriptome_snps, output_file, current_samples, k)

def parse_arguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f','--transcriptome-fasta', required=True, type=str)
  parser.add_argument('-v','--transcriptome-vcf', required=True, type=str)
  parser.add_argument('-t','--output-type', required=True, type=str, choices=['diploid', 'ornament'])
  parser.add_argument('-o','--output-file', required=True, type=str)
  parser.add_argument('-s','--samples', required=True, type=str, help='comma separated list, or "ALL"')
  parser.add_argument('-c', '--coordinates', type=str)
  parser.add_argument('-k', '--kmer-length', default=31, type=int)
  return parser.parse_args()

if __name__ == "__main__":
  args = parse_arguments()
  main(args.transcriptome_fasta, args.transcriptome_vcf, args.output_type, args.output_file, args.samples, args.coordinates, args.kmer_length)
