#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdint.h>
#include <ostream>
//#include <map>


#include "common.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"

#include "KmerHashTable.h"

#include "hash.hpp"

std::string revcomp(const std::string s);


struct TRInfo {
  int trid;
  int start;
  int stop; //exclusive [start,stop)
  bool sense; // true for sense, false for anti-sense
};

using EcMap = std::vector<std::vector<int>>; //std::unordered_map<int, std::vector<int>>;

struct SortedVectorHasher {
  size_t operator()(const std::vector<int>& v) const {
    uint64_t r = 0;
    int i=0;
    for (auto x : v) {
      uint64_t t;
      MurmurHash3_x64_64(&x,sizeof(x), 0,&t);
      t = (x>>i) | (x<<(64-i));
      r = r ^ t;
      i = (i+1)%64;
    }
    return r;
  }
};

struct SortedVectorPairHasher {
  size_t operator()(const std::pair<std::vector<int>, std::vector<int>>& p) const {
    uintmax_t hash = SortedVectorHasher{}(p.first);
    uintmax_t hash2 = SortedVectorHasher{}(p.second);
    hash ^= hash2 + 0x9e3779b9 + (hash<<6) + (hash>>2);
    return hash;
  }
};

struct KmerEntry {
  int32_t contig; // id of contig
  uint32_t _pos; // 0-based forward distance to EC-junction
  int32_t contig_length;

  KmerEntry() : contig(-1), _pos(0xFFFFFFF), contig_length(0) {}
  KmerEntry(int id, int length, int pos, bool isFw) : contig(id), contig_length(length) {
    setPos(pos);
    setDir(isFw);
  }

  inline int getPos() const {return (_pos & 0x0FFFFFFF);}
  inline int isFw() const  {return (_pos & 0xF0000000) == 0; }
  inline void setPos(int p) {_pos = (_pos & 0xF0000000) | (p & 0x0FFFFFFF);}
  inline void setDir(bool _isFw) {_pos = (_pos & 0x0FFFFFFF) | ((_isFw) ? 0 : 0xF0000000);}
  inline int getDist(bool fw) const {
    if (isFw() == fw) {
      return (contig_length - 1 - getPos());
    } else {
      return getPos();
    }
  }
};

struct ContigToTranscript {
  int trid;
  int pos; 
  bool sense; // true for sense, 
};

struct Contig {
  int id; // internal id
  int length; // number of k-mers
  int ec;
  std::string seq; // sequence
  std::vector<ContigToTranscript> transcripts;
};

struct DBGraph {
  std::vector<int> ecs; // contig id -> ec-id
  std::vector<Contig> contigs; // contig id -> contig
//  std::vector<pair<int, bool>> edges; // contig id -> edges
};

struct VariantParameterEntry {
    int trid;
    int pos;
    char orientation;
    std::string ref;
    std::string alt;
    bool is_first;


    VariantParameterEntry(int trid, int pos, char orientation, std::string ref, std::string alt) :
      trid(trid),
      pos(pos),
      orientation(orientation),
      ref(ref),
      alt(alt),
      is_first(true) {

    }

    bool is_indel() const {
      return (ref.size() != alt.size());
    }

    int delta_length() const {
      return (alt.size() - ref.size());
    }

    std::string toString() const {
        return 
          std::to_string(trid) + " " + 
          std::to_string(pos) +
          std::string(1, orientation) + " " +
          (is_indel() ? "indel" : "snp");
    }
};

struct KmerIndex {
  KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0), skip(opt.skip), target_seqs_loaded(false) {
    //LoadTranscripts(opt.transfasta);
  }

  ~KmerIndex() {}

  void match(const char *s, int l, std::vector<std::pair<KmerEntry, int>>& v) const;
//  bool matchEnd(const char *s, int l, std::vector<std::pair<int, int>>& v, int p) const;
  int mapPair(const char *s1, int l1, const char *s2, int l2, int ec) const;
  std::vector<int> intersect(int ec, const std::vector<int>& v) const;

  void convert_string_to_arr(std::string variant_info, std::vector<int>& arr) const;
  void convert_ornament_to_snps(std::string ornament_name, std::vector<int>& refsnp_vals, std::vector<int>& altsnp_vals, std::vector<int>& alt2snp_vals, int* start_pos) const;
  int get_ornament_position(int shade) const;
  std::vector<size_t> get_substring_positions(std::string t, std::string ornament, int orn_tid, int k);
  bool get_array_matches(int start_pos, std::string t, std::string ornament, int orn_tid, int k);
  bool is_problem_transcript(std::string t, std::string ornament, int orn_tid, int k);

  void BuildTranscripts(const ProgramOptions& opt);
  void BuildDeBruijnGraph(const ProgramOptions& opt, const std::vector<std::string>& seqs);
  void BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs);
  void FixSplitContigs(const ProgramOptions& opt, std::vector<std::vector<TRInfo>>& trinfos);
  bool fwStep(Kmer km, Kmer& end) const;

  // output methods
  void write(const std::string& index_out, bool writeKmerTable = true);
  void writePseudoBamHeader(std::ostream &o) const;
  
  // note opt is not const
  // load methods
  void load(ProgramOptions& opt, bool loadKmerTable = true);
  void loadTranscriptSequences() const;
  void initVariantParamEntry(std::string vcf_file, std::string sample_name, bool is_phased);
  void initRepeatTranscriptBiasCorrection(std::string bias_file);
  void filterVariants(const ProgramOptions& opt, bool filter_phased);
  void clear();

  // positional information
  std::pair<int, bool> findFirstOrnamentAgnosticPosition(int tr, Kmer km, KmerEntry& val, int p) const;
  std::pair<int,bool> findPosition(int tr, Kmer km, KmerEntry val, int p = 0) const;
  std::pair<int,bool> findPosition(int tr, Kmer km, int p) const;

  std::vector<std::string> getTranscriptNamesForKmer(const char* kmer);

  int k; // k-mer size used
  int num_trans; // number of targets
  int skip;

  KmerHashTable<KmerEntry, KmerHash> kmap;
  EcMap ecmap;
  DBGraph dbGraph;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;
  std::vector<std::pair<std::vector<int>, std::vector<int>>> vcmap;
  std::unordered_map<std::pair<std::vector<int>, std::vector<int>>, int, SortedVectorPairHasher> vcmapinv;
  std::unordered_map<int, int> shadeToColorTranscriptMap;

  std::unordered_map<int, std::vector<VariantParameterEntry>> effectiveASEParams;
  std::unordered_map<int, std::unordered_map<int, int>> transcriptToParam;

  std::unordered_map<int, std::unordered_set<int>> homo_alt_snps;
  const size_t INDEX_VERSION = 10; // increase this every time you change the fileformat

  std::vector<int> target_lens_;
  std::vector<int> target_left_lens_;
  std::vector<int> target_right_lens_;

  std::vector<std::string> target_names_;
  std::vector<std::string> target_seqs_; // populated on demand
  bool target_seqs_loaded;


};








#endif // KALLISTO_KMERINDEX_H
