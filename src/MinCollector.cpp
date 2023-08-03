#include "MinCollector.h"
#include <algorithm>
#include <unordered_set>

// utility functions

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"}";
    return os;
}

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y) {
  std::vector<int> v;
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end()) {
    if (*a < *b) {
      ++a;
    } else if (*b < *a) {
      ++b;
    } else {
      v.push_back(*a);
      ++a;
      ++b;
    }
  }
  return v;
}

void MinCollector::init_mean_fl_trunc(double mean, double sd) {
  auto tmp_trunc_fl = trunc_gaussian_fld(0, MAX_FRAG_LEN, mean, sd);
  assert( tmp_trunc_fl.size() == mean_fl_trunc.size() );

  std::copy( tmp_trunc_fl.begin(), tmp_trunc_fl.end(), mean_fl_trunc.begin() );

  mean_fl = mean_fl_trunc[ MAX_FRAG_LEN - 1 ];

  has_mean_fl = true;
  has_mean_fl_trunc = true;
}
// v1 is vector of pairs of KmerEntries and their positions in the read
int MinCollector::intersectKmers(std::vector<std::pair<KmerEntry,int>>& v1,
                          std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired, std::vector<int> &u) const {

  // u1 is for single end reads, u2 is additionally used for paired end reads
  std::vector<int> u1 = intersectECs(v1);
  std::vector<int> u2 = intersectECs(v2);

  if (u1.empty() && u2.empty()) {
    return -1;
  }

  // non-strict intersection.
  if (u1.empty()) {
    if (v1.empty()) {
      u = u2;
    } else {
      return -1;
    }
  } else if (u2.empty()) {
    if (v2.empty()) {
      u = u1;
    } else {
      return -1;
    }
  } else {
    u = intersect(u1,u2);
  }

  if (u.empty()) {
    return -1;
  }
  return 1;
}

int MinCollector::collect(std::vector<std::pair<KmerEntry,int>>& v1,
                          std::vector<std::pair<KmerEntry,int>>& v2, bool nonpaired) {
  std::vector<int> u;
  int r = intersectKmers(v1, v2, nonpaired, u);
  if (r != -1) {
    return increaseCount(u);
  } else {
    return -1;
  }
}

// u is list of transcripts, v is list of effective parameters that the shades map to
// finds VC, or places it in and returns it if not found
// change this to not modify index within findVC
int MinCollector::findVC(const std::pair<std::vector<int>, std::vector<int>>& p) const {
  auto search = index.vcmapinv.find(p);
  if (search != index.vcmapinv.end()) {
    return search->second;
  } else {
    return -1;
  }
}

int MinCollector::increaseVCount(const variant_class& p) {
  int vc = findVC(p);
  if (vc != -1) {
    ++aseCounts[vc];
  } else {
    int size = index.vcmap.size();
    index.vcmapinv.insert({p, size});
    index.vcmap.push_back(p);
    aseCounts[size] = 1;
    vc = size;
  }
  return vc;
}

int MinCollector::findEC(const std::vector<int>& u) const {
  if (u.empty()) {
    return -1;
  }
  if (u.size() == 1) {
    return u[0];
  }
  auto search = index.ecmapinv.find(u);
  if (search != index.ecmapinv.end()) {
    return search ->second;
  } else {
    return -1;
  }
}

int MinCollector::increaseCount(const std::vector<int>& u) {
  int ec = findEC(u);

  if (u.empty()) {
    return -1;
  } else {
    if (ec != -1) {
      ++counts[ec];
      return ec;
    } else {
      auto necs = counts.size();
      //index.ecmap.insert({necs,u});
      index.ecmap.push_back(u);
      index.ecmapinv.insert({u,necs});
      counts.push_back(1);
      return necs;
    }
  }

  /* -- removable
  if (u.size() == 1) {
    int ec = u[0];
    ++counts[ec];
    return ec;
  }
  auto search = index.ecmapinv.find(u);
  if (search != index.ecmapinv.end()) {
    // ec class already exists, update count
    ++counts[search->second];
    return search->second;
  } else {
    // new ec class, update the index and count
    auto necs = counts.size();
    //index.ecmap.insert({necs,u});
    index.ecmap.push_back(u);
    index.ecmapinv.insert({u,necs});
    counts.push_back(1);
    return necs;
  }
  */
}

int MinCollector::decreaseCount(const int ec) {
  assert(ec >= 0 && ec <= index.ecmap.size());
  --counts[ec];
  return ec;
}

struct ComparePairsBySecond {
  bool operator()(std::pair<KmerEntry,int> a, std::pair<KmerEntry,int> b) {
    return a.second < b.second;
  }
};

int getPositionInTranscript(Contig& c, KmerEntry& ke, Kmer& km, ContigToTranscript& ctt) {
  int k = 31;
  bool fw = (km == km.rep());
  bool csense = (fw == ke.isFw());
  int trpos = ctt.pos;
  bool trsense = ctt.sense;
  if (trsense) {
    if (csense) {
      trpos = trpos + ke.getPos() + 1; // 1-based, case I
    } else {
      trpos = trpos + ke.getPos() + k ; // 1-based, case III
    }
  } else {
    if (csense) {
      trpos = trpos + (c.length - ke.getPos() -1) + k;  // 1-based, case IV
    } else {
      trpos =  trpos + (c.length - ke.getPos()); // 1-based, case II
    }
  }
  return trpos;
}

/* Return 1 if position for this shade is uniquely implied, 0 if compatible, -1 if uncompatible */
/* Can also cache result of this */
int MinCollector::isUniqueOrCompatible(
    std::vector<std::pair<KmerEntry,int>>& v,
    int layer,
    bool fw,
    bool strand,
    int color,
    int ornament_pos,
    int shade,
    const char* s) const {
  /* Scan in the forward direction */
  if (fw) {
    /* Go in the forward direction */
    if (layer < v.size() - 1) {
      Kmer km((s + v[layer + 1].second));
      KmerEntry ke = v[layer + 1].first;
      int seen_count = 0;
      bool is_compat = false;
      int add_or_subtract = (strand) ? 1 : -1;
      Contig& c = index.dbGraph.contigs[ke.contig];
      for (auto& ctt : c.transcripts) {
        if (ctt.trid == color) {
          seen_count += 1;
          if (ctt.trid == color && std::fabs(
                getPositionInTranscript(c, ke, km, ctt) -
                add_or_subtract * (v[layer + 1].second - v[layer].second) - 
                ornament_pos) < 3) {
            is_compat = true;
          }
        }
      }

      if (!is_compat) {
        return -1;
      }

      int curr_ret = (seen_count == 1) ? (is_compat ? 1 : -1) : (is_compat ? 0 : -1);
      if (curr_ret == 1) {
        return curr_ret;
      } else {
        int future = isUniqueOrCompatible(
            v,
            layer + 1,
            fw,
            strand,
            color,
            ornament_pos + add_or_subtract * (v[layer + 1].second - v[layer].second),
            shade,
            s);
        return future;
      }
    } else {
      return 0;
    }
  } else {
    /* Go in the backwards direction */
    if (layer > 0) {
      Kmer km((s + v[layer - 1].second));
      KmerEntry ke = v[layer - 1].first;
      int seen_count = 0;
      bool is_compat = false;
      int add_or_subtract = (strand) ? 1 : -1;
      Contig& c = index.dbGraph.contigs[ke.contig];
      for (auto& ctt : c.transcripts) {
        if (ctt.trid == color) {
          seen_count += 1;
          if (ctt.trid == color && std::fabs(
                getPositionInTranscript(c, ke, km, ctt) +
                add_or_subtract * (v[layer].second - v[layer - 1].second) - 
                ornament_pos) < 3) {
            is_compat = true;
          }
        }
      }

      if (!is_compat) {
        return -1;
      }

      int curr_ret = (seen_count == 1) ? (is_compat ? 1 : -1) : (is_compat ? 0 : -1);
      if (curr_ret == 1) {
        return curr_ret;
      } else {
        int past = isUniqueOrCompatible(
            v,
            layer - 1,
            fw,
            strand,
            color,
            ornament_pos - add_or_subtract * (v[layer].second - v[layer - 1].second),
            shade,
            s);
        return past;
      }
    } else {
      return 0;
    }
  }
}

enum RefInfo { ref, alt, alt2, nonvariant };

RefInfo refOrAlt(std::string t) {
  int refcount = 0;
  int altcount = 0;
  int alt2count = 0;
  for (int i = 0; i < t.size(); i++) {
    if (t[i] == 'r') {
      refcount += 1;
    } else if (t[i] == 'a') {
      altcount += 1;
    } else if (t[i] == 'y') {
      alt2count += 1;
    }
  }

  if (refcount > altcount && refcount > alt2count) {
    return RefInfo::ref;
  } else if (altcount > refcount && altcount > alt2count) {
    return RefInfo::alt;
  } else if (alt2count > refcount && alt2count > altcount) {
    return RefInfo::alt2;
  } else {
    return RefInfo::nonvariant;
  }
}

// return a dictionary mapping the transcript name to the shades for it
// need to optimize, but full run over a set of reads takes around 10 minutes
/* this usage of variant class is slightly different, keep track of rsnps (0) and asnps(1) */
/* TODO: Make this one map that is modified upon input, and cleared for every read in processread? */
void MinCollector::unionShades(
    std::vector<std::pair<KmerEntry,int>>& v,
    std::unordered_map<int, variant_info>& variantInfo,
    std::vector<int>& u,
    std::unordered_map<int, int>& variantU,
    const char* s) const {
  if (v.empty() || u.size() == 0) {
    return;
  }
  /* Sort the equivalence classes by read position */
  sort(v.begin(), v.end(), [&](std::pair<KmerEntry, int> a, std::pair<KmerEntry, int> b)
      {
      return a.second < b.second;
      });

  KmerEntry fi = v[0].first;
  Kmer kfi = Kmer((s + v[0].second));
  bool strand = (fi.isFw() == (kfi == kfi.rep())); // k-mer maps to fw strand?
  Contig& cfi = index.dbGraph.contigs[fi.contig];
  strand = (strand == cfi.transcripts[0].sense);
  int read_len = strlen(s);
  int k = 31;

  /* Assume that the strand is strand specific? */

    // shades will contain all the variant shades that this read pseudoaligned to
  std::vector<int> shades;
  /* Iterate through the equivalence classes, and flag the ones with repeat trouble */
  /* Make sure the ones we add are within the pseudoalignment */
  std::unordered_set<int> spurious_ornaments;
  std::unordered_set<int> longer_repeat_ornaments;
  std::unordered_set<int> calculated_ornaments;
  /* Take a pass through and see if there are any transcripts that appear more than once? */
  std::vector<std::vector<int>> shades_across_vs;
  std::vector<int> shade_positions;
  int ppos = -1;
  for (int i = 0; i < v.size(); i++) {
    int read_position = v[i].second;
    if (read_position == ppos || read_len - read_position < k) {
      continue;
    }
    KmerEntry ke = v[i].first;
    Kmer km((s + read_position));
    shades_across_vs.push_back({});
    Contig& c = index.dbGraph.contigs[ke.contig];
    /* Implementation is not safe against having multiple of the same shade within this array */
    for (ContigToTranscript& ctt : c.transcripts) {
      /* This is a list of all occurences of this kmer within the transcriptome */
      /* If this is a shade, then check for compatibility for current position */
      auto search = index.shadeToColorTranscriptMap.find(ctt.trid);
      if (search != index.shadeToColorTranscriptMap.end()) {
        int ornament_position = index.get_ornament_position(ctt.trid) + getPositionInTranscript(c, ke, km, ctt);
        int color = search->second;
        if (std::find(u.begin(), u.end(), color) == u.end()) {
          continue;
        }
        int num_color = 0;
        std::vector<int> obs;
        for (auto& rep : c.transcripts) {
          if (rep.trid == color && std::find(obs.begin(), obs.end(), rep.pos) == obs.end()) {
            obs.push_back(rep.pos);
            num_color++;
          }
        }
        /* This is a repeat region */
        if (num_color > 1) {
          /* This shade has already been cached. Find out if it's legit or not, and if so, add it */
          if (calculated_ornaments.find(ctt.trid) != calculated_ornaments.end()) {
            if (spurious_ornaments.find(ctt.trid) == spurious_ornaments.end() &&
                longer_repeat_ornaments.find(ctt.trid) == longer_repeat_ornaments.end()) {
              shades_across_vs[shades_across_vs.size() - 1].push_back(ctt.trid);
              if (shade_positions.size() == 0 || shade_positions[shade_positions.size() - 1] != v[i].second) {
                shade_positions.push_back(v[i].second);
              }
            }
            continue;
          }
          int forward = isUniqueOrCompatible(
              v, i, true, (strand), color, ornament_position, ctt.trid, s);
          int backward = isUniqueOrCompatible(
              v, i, false, (strand), color, ornament_position, ctt.trid, s);
          calculated_ornaments.insert(ctt.trid);
          /* If one of these are 1 and the other is zero, then valid */
          /* If one of these are 1 and the other is -1, then that's weird */
          /* If both of these are -1 then reject */
          /* If both of these are 0 then unknown */
          if (backward == -1 || forward == -1) {
            spurious_ornaments.insert(ctt.trid);
          } else if (backward != 1 && forward != 1) {
            longer_repeat_ornaments.insert(ctt.trid);
          } else {
            shades_across_vs[shades_across_vs.size() - 1].push_back(ctt.trid);
            if (shade_positions.size() == 0 || shade_positions[shade_positions.size() - 1] != v[i].second) {
              shade_positions.push_back(v[i].second);
            }
          }
        } else {
          /* This is the normal case -- no repeats, union across shades and intersect within variants */
          shades_across_vs[shades_across_vs.size() - 1].push_back(ctt.trid);
          if (shade_positions.size() == 0 || shade_positions[shade_positions.size() - 1] != v[i].second) {
            shade_positions.push_back(v[i].second);
          }
        }
      }
    }
    if (shades_across_vs[shades_across_vs.size() - 1].size() == 0) {
      shades_across_vs.pop_back();
    }
    ppos = read_position;
  }

  if (shades_across_vs.size() == 0) return;

  std::unordered_map<int, std::unordered_map<int, std::string>> variantsIntersected;
  std::vector<int> refsnpVals;
  std::vector<int> altsnpVals;
  std::vector<int> alt2snpVals;
  for (int level = 0; level < shades_across_vs.size(); level++) {
    int potential_weight = 1;
    if (level > 0) {
      potential_weight = shade_positions[level] - shade_positions[level - 1];
    }

    for (int shade : shades_across_vs[level]) {
      int color = index.shadeToColorTranscriptMap.at(shade);
      int start_pos;
      std::string ornament_name = index.target_names_[shade];
      auto& tmap = variantsIntersected[color];
      refsnpVals.clear();
      altsnpVals.clear();
      alt2snpVals.clear();
      index.convert_ornament_to_snps(ornament_name, refsnpVals, altsnpVals, alt2snpVals, &start_pos);

      for (int r : refsnpVals) {
        int pcount = 0;
        bool had_r = false;
        for (int q = tmap[r].size() - 1; q >= 0; q--) { 
          if (tmap[r][q] == 'r') {
            had_r = true;
          }  else if (tmap[r][q] == '.') {
            pcount += 1;
            if (pcount > 1) {
              break;
            }
          }
        }
        std::string to_add((had_r) ? potential_weight : 1, 'r');
        tmap[r] += to_add;
      }
      for (int a : altsnpVals) {
        int pcount = 0;
        bool had_a = false;
        for (int q = tmap[a].size() - 1; q >= 0; q--) { 
          if (tmap[a][q] == 'r') {
            had_a = true;
          }  else if (tmap[a][q] == '.') {
            pcount += 1;
            if (pcount > 1) {
              break;
            }
          }
        }
        std::string to_add((had_a) ? potential_weight : 1, 'a');
        tmap[a] += to_add;
      }
      for (int y : alt2snpVals) {
        int pcount = 0;
        bool had_y = false;
        for (int q = tmap[y].size() - 1; q >= 0; q--) { 
          if (tmap[y][q] == 'r') {
            had_y = true;
          }  else if (tmap[y][q] == '.') {
            pcount += 1;
            if (pcount > 1) {
              break;
            }
          }
        }
        std::string to_add((had_y) ? potential_weight : 1, 'y');
        tmap[y] += to_add;
      }
    }
    for (auto& tmap : variantsIntersected) {
      for (auto& k : tmap.second) {
        k.second += ".";
      }
    }
  }

  for (auto& tmap : variantsIntersected) {
    for (auto& k : tmap.second) {
      RefInfo curr = refOrAlt(k.second);
      if (curr == RefInfo::ref) {
        std::vector<int>& arr = std::get<0>(variantInfo[tmap.first]);
        arr.push_back(k.first);
      } else if (curr == RefInfo::alt) {
        std::vector<int>& arr = std::get<1>(variantInfo[tmap.first]);
        arr.push_back(k.first);
      } else if (curr == RefInfo::alt2) {
        std::vector<int>& arr = std::get<2>(variantInfo[tmap.first]);
        arr.push_back(k.first);
      }
    }
  }
}

std::vector<int> MinCollector::intersectECs(std::vector<std::pair<KmerEntry,int>>& v) const {
  if (v.empty()) {
    return {};
  }
  
  // Filter these here to ignore contigs with only one variant shade
  sort(v.begin(), v.end(), [&](std::pair<KmerEntry, int> a, std::pair<KmerEntry, int> b)
       {
         if (a.first.contig==b.first.contig) {
           return a.second < b.second;
         } else {
           return a.first.contig < b.first.contig;
         }
       }); // sort by contig, and then first position


  
  int ec = index.dbGraph.ecs[v[0].first.contig];
  int lastEC = ec;
  std::vector<int> u = index.ecmap[ec];

  for (int i = 1; i < v.size(); i++) {
    if (v[i].first.contig != v[i-1].first.contig) {
      ec = index.dbGraph.ecs[v[i].first.contig];
      if (ec != lastEC) {
        u = index.intersect(ec, u);
        lastEC = ec;
        if (u.empty()) {
          return u;
        }
      }
    }
  }

  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;

  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }

  if ((maxpos-minpos + k) < min_range) {
    return {};
  }

  return u;
}


void MinCollector::loadCounts(ProgramOptions& opt) {
  int num_ecs = counts.size();
  counts.clear();
  std::ifstream in((opt.output + "/counts.txt"));
  int i = 0;
  if (in.is_open()) {
    std::string line;
    while (getline(in, line)) {
      std::stringstream ss(line);
      int j,c;
      ss >> j;
      ss >> c;
      if (j != i) {
        std::cerr << "Error: equivalence class does not match index. Found "
                  << j << ", expected " << i << std::endl;
        exit(1);
      }
      counts.push_back(c);
      i++;
    }

    if (i != num_ecs) {
      std::cerr << "Error: number of equivalence classes does not match index. Found "
                << i << ", expected " << num_ecs << std::endl;
      exit(1);
    }
  } else {
    std::cerr << "Error: could not open file " << opt.output << "/counts.txt" << std::endl;
    exit(1);

  }

}

void MinCollector::write(const std::string& pseudoprefix) const {
  std::string ecfilename = pseudoprefix + ".ec";
  std::string countsfilename = pseudoprefix + ".tsv";

  std::ofstream ecof, countsof;
  ecof.open(ecfilename.c_str(), std::ios::out);
  // output equivalence classes in the form "EC TXLIST";
  for (int i = 0; i < index.ecmap.size(); i++) {
    ecof << i << "\t";
    // output the rest of the class
    const auto &v = index.ecmap[i];
    bool first = true;
    for (auto x : v) {
      if (!first) {
        ecof << ",";
      } else {
        first = false;
      }
      ecof << x;
    }
    ecof << "\n";
  }
  ecof.close();

  countsof.open(countsfilename.c_str(), std::ios::out);
  for (int i = 0; i < counts.size(); i++) {
    countsof << i << "\t" << counts[i] << "\n";
  }
  countsof.close();
}

double MinCollector::get_mean_frag_len(bool lenient) const {
  if (has_mean_fl) {
    return mean_fl;
  }

  auto total_counts = 0;
  double total_mass = 0.0;

  for ( size_t i = 0 ; i < flens.size(); ++i ) {
    total_counts += flens[i];
    total_mass += static_cast<double>(flens[i] * i);
  }

  if (total_counts == 0) {
    if (!lenient) {
      std::cerr << "Error: could not determine mean fragment length from paired end reads, no pairs mapped to a unique transcript." << std::endl
              << "       Run kallisto quant again with a pre-specified fragment length (option -l)." << std::endl;
      exit(1);
    } else {
      return std::numeric_limits<double>::max();
    }

  }

  // cache the value
  const_cast<double&>(mean_fl) = total_mass / static_cast<double>(total_counts);
  const_cast<bool&>(has_mean_fl) = true;
  return mean_fl;
}

double MinCollector::get_sd_frag_len() const {
  double tmp = get_mean_frag_len(true);
  double m = mean_fl;

  size_t total_counts = 0;
  double total_mass = 0.0;
  
  for (size_t i = 0; i < flens.size(); ++i) {
    total_counts += flens[i];
    total_mass += flens[i]*(i-m)*(i-m);
  }

  double sd_fl = std::sqrt(total_mass/total_counts);
  return sd_fl;
}



void MinCollector::compute_mean_frag_lens_trunc(bool verbose)  {

  std::vector<int> counts(MAX_FRAG_LEN, 0);
  std::vector<double> mass(MAX_FRAG_LEN, 0.0);

  counts[0] = flens[0];

  for (size_t i = 1; i < MAX_FRAG_LEN; ++i) {
    // mass and counts keep track of the mass/counts up to and including index i
    mass[i] = static_cast<double>( flens[i] * i) + mass[i-1];
    counts[i] = flens[i] + counts[i-1];
    if (counts[i] > 0) {
      mean_fl_trunc[i] = mass[i] / static_cast<double>(counts[i]);
    }
    // std::cerr << "--- " << i << '\t' << mean_fl_trunc[i] << std::endl;
  }

  has_mean_fl_trunc = true;

  if (verbose) {
    std::cerr << "[quant] estimated average fragment length: " <<
      mean_fl_trunc[MAX_FRAG_LEN - 1] << std::endl;
  }
}

int hexamerToInt(const char *s, bool revcomp) {
  int hex = 0;
  if (!revcomp) {
    for (int i = 0; i < 6; i++) {
      hex <<= 2;
      switch (*(s+i) & 0xDF) {
      case 'A': break;
      case 'C': hex += 1; break;
      case 'G': hex += 2; break;
      case 'T': hex += 3; break;
      default: return -1;
      }
    }
  } else {
    for (int i = 0; i < 6; i++) {
      switch (*(s+i) & 0xDF) {
      case 'A': hex += 3 << (2*i);break;
      case 'C': hex += 2 << (2*i); break;
      case 'G': hex += 1 << (2*i); break;
      case 'T': break;
      default: return -1;
      }
    }
  }
  return hex;
}

bool MinCollector::countBias(const char *s1, const char *s2, const std::vector<std::pair<KmerEntry,int>> v1, const std::vector<std::pair<KmerEntry,int>> v2, bool paired) {
  return countBias(s1,s2,v1,v2,paired,bias5);
}

bool MinCollector::countBias(const char *s1, const char *s2, const std::vector<std::pair<KmerEntry,int>> v1, const std::vector<std::pair<KmerEntry,int>> v2, bool paired, std::vector<int>& biasOut) const {

  const int pre = 2, post = 4;

  if (v1.empty() || (paired && v2.empty())) {
    return false;
  }

  

  auto getPreSeq = [&](const char *s, Kmer km, bool fw, bool csense,  KmerEntry val, int p) -> int {
    if (s==0) {
      return -1;
    }
    if ((csense && val.getPos() - p >= pre) || (!csense && (val.contig_length - 1 - val.getPos() - p) >= pre )) {
      const Contig &c = index.dbGraph.contigs[val.contig];
      bool sense = c.transcripts[0].sense;

      int hex = -1;
      //std::cout << "  " << s << "\n";
      if (csense) {
        hex = hexamerToInt(c.seq.c_str() + (val.getPos()-p - pre), true);
        //std::cout << c.seq.substr(val.getPos()- p - pre,6) << "\n";
      } else {
        int pos = (val.getPos() + p) + k - post;
        hex = hexamerToInt(c.seq.c_str() + (pos), false);
        //std::cout << revcomp(c.seq.substr(pos,6)) << "\n";
      }
      return hex;
    }

    return -1;
  };

  // find first contig of read
  KmerEntry val1 = v1[0].first;
  int p1 = v1[0].second;
  for (auto &x : v1) {
    if (x.second < p1) {
      val1 = x.first;
      p1 = x.second;
    }
  }

  Kmer km1 = Kmer((s1+p1));
  bool fw1 = (km1==km1.rep());
  bool csense1 = (fw1 == val1.isFw()); // is this in the direction of the contig?

  int hex5 = getPreSeq(s1, km1, fw1, csense1, val1, p1);

  /*
  int hex3 = -1;
  if (paired) {
    // do the same for the second read
    KmerEntry val2 = v2[0].first;
    int p2 = v2[0].second;
    for (auto &x : v2) {
      if (x.second < p2) {
        val2 = x.first;
        p2 = x.second;
      }
    }

    Kmer km2 = Kmer((s2+p2));
    bool fw2 = (km2==km2.rep());
    bool csense2 = (fw2 == val2.isFw()); // is this in the direction of the contig?

    hex3 = getPreSeq(s2, km2, fw2, csense2, val2, p2);
  }
  */

  if (hex5 >=0) { // && (!paired || hex3 >= 0)) {
    biasOut[hex5]++;
    //bias3[hex3]++;
  } else {
    return false;
  }

  return false;
}
