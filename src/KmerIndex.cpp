#include "KmerIndex.h"
#include <algorithm>
#include <random>
#include <ctype.h>
#include <zlib.h>
#include <unordered_set>
#include "kseq.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"}";
    return os;
}

static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
    return !std::isspace(ch);
  }));
}

static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

// helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] == u[j]) {
      return false;
    }
  }
  return true;
}


std::vector<int> unique(const std::vector<int>& u) {
  std::vector<int> v;
  v.reserve(u.size());
  v.push_back(u[0]);
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] != u[j]) {
      v.push_back(u[j]);
    }
  }
  return v;
}

const char Dna(int i) {
  static const char *dna = "ACGT";
  return dna[i & 0x03];
}

int hamming(const char *a, const char *b) {
  int h = 0;
  while (*a != 0 && *b != 0) {
    if (*a != *b) {
      h++;
    }
    a++;
    b++;
  }
  return h;
}

std::string revcomp(const std::string s) {
  std::string r(s);
  std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
      switch(c) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      default: return 'N';
      }
      return 'N';
    });
  return r;
}

int num_eq(int start_pos, int num_comp, int k, std::string s1, std::string s2) {
  int num_same = 0;
  for (int i = start_pos; i < start_pos + num_comp; i++) {
    if (s1[i] == s2[i - start_pos]) {
      num_same += 1;
    }
  }
  return num_same;
}

bool KmerIndex::get_array_matches(int start_pos, std::string t, std::string ornament, int orn_tid, int k) {
  int first_match = 0;
  std::vector<int> matches;
  for (int i = 0; i < k; i++) {
    if (t[i + start_pos] == ornament[i]) {
      first_match += 1;
    }
  }
  bool potentially_problem = false;
  if (first_match == 31) {
    potentially_problem = true;
  }
  matches.push_back(first_match);

  for (int i = k; i < ornament.size(); i++) {
    if (t[i + start_pos - k] == ornament[i - k]) {
      first_match -= 1;
    }

    if (t[i + start_pos] == ornament[i]) {
      first_match += 1;
    }

    if (first_match == 31) {
      potentially_problem = true;
    }
    matches.push_back(first_match);
  }

  if (potentially_problem) {
    // If this is all equal to k, and if the start position is the same as orn tid's start position, then no
    size_t pos_in_t = get_ornament_position(orn_tid);
    int all_k = true;
    for (int e : matches) {
      all_k &= (e == k);
    }

    if (all_k) {
      if (pos_in_t == start_pos - (ornament.size() - k)) {
        return false;
      }
    }
  }

  return potentially_problem;
}

bool KmerIndex::is_problem_transcript(std::string t, std::string ornament, int orn_tid, int k) {
  if (t.size() < k || ornament.size() < k) {
    return false;
  }
  std::string pad(ornament.size() - k, '$');
  std::string padded_t = pad + t + pad;
  bool is_problem = false;
  for (int i = 0; i < padded_t.size() - ornament.size(); i++) {
    is_problem |= get_array_matches(i, padded_t, ornament, orn_tid, k); 
  }
  return is_problem;
}



std::vector<size_t> KmerIndex::get_substring_positions(std::string t, std::string ornament, int orn_tid, int k) {
  if (t.size() < k || ornament.size() < k) {
    return {};
  }
  size_t orn_position = 0;
  size_t pos_in_t = get_ornament_position(orn_tid);
  const char* transcript_seq = t.c_str();
  const char* ornament_seq = ornament.c_str();
  std::vector<size_t> repeat_positions;
  /* need helper to get position of ornament in haploid transcript */
  for (int j = 0; j < ornament.size() - k; j++) {
    /* ornament[j:j+k] is a k-mer to consider */
    bool timeit = false;
    for (int i = 0; i < t.size() - k; i++) {
      /* search for this k-mer within the transcript */
      /* disregard if its at the exact same position as implied by the ornament */
      bool appears_here = false;
      if (timeit) {
        /* t[i]; */
        /* if t[i + k - 1]; */
        continue;
      }
      if (::strncmp(transcript_seq + i, ornament_seq + j, k) == 0) {
        appears_here = true;
      }

      timeit = true;
      if (appears_here && i != j + pos_in_t) {
        repeat_positions.push_back(i);
        /* return repeat_positions; */
      }
    }
  }
  return repeat_positions;
}

void KmerIndex::BuildTranscripts(const ProgramOptions& opt) {
  // read input
  std::unordered_set<std::string> unique_names;
  int k = opt.k;
  for (auto& fasta : opt.transfasta) {
    std::cerr << "[build] loading fasta file " << fasta
              << std::endl;
  }
  std::cerr << "[build] k-mer length: " << k << std::endl;


  std::vector<std::string> seqs;

  // read fasta file  
  gzFile fp = 0;
  int aw_fuck = 0;
  kseq_t *seq;
  int l = 0;
  std::mt19937 gen(42);
  int countNonNucl = 0;
  int countUNuc = 0;
  int polyAcount = 0;
  // keep track of the last non-variant transcript
  int last_full_tr = -1;
  time_t start, end;
  time(&start);
  for (auto& fasta : opt.transfasta) {
    fp = gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);

    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      seqs.emplace_back(seq->seq.s);
      std::string& str = *seqs.rbegin();
      auto n = str.size();
      for (auto i = 0; i < n; i++) {
        char c = str[i];
        c = ::toupper(c);
        if (c=='U') {
          str[i] = 'T';
          countUNuc++;
        } else if (c !='A' && c != 'C' && c != 'G' && c != 'T') {
          str[i] = Dna(gen()); // replace with pseudorandom string
          countNonNucl++;
        }
      }
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);

      if (str.size() >= 10 && str.substr(str.size()-10,10) == "AAAAAAAAAA") {
        // clip off polyA tail
        polyAcount++;
        int j;
        for (j = str.size()-1; j >= 0 && str[j] == 'A'; j--) {}
        str = str.substr(0,j+1);
      }


      target_lens_.push_back(seq->seq.l);
      std::string name(seq->name.s);
      size_t p = name.find(' ');
      if (p != std::string::npos) {
        name = name.substr(0,p);
      }

      if (unique_names.find(name) != unique_names.end()) {
        if (!opt.make_unique) {
          std::cerr << "Error: repeated name in FASTA file " << fasta << "\n" << name << "\n\n" << "Run with --make-unique to replace repeated names with unique names" << std::endl;
          exit(1);
        } else {
          for (int i = 1; ; i++) { // potential bug if you have more than 2^32 repeated names
            std::string new_name = name + "_" + std::to_string(i);
            if (unique_names.find(new_name) == unique_names.end()) {
              name = new_name;
              break;
            }
          }
        }
      }
      unique_names.insert(name);
      target_names_.push_back(name);

      // Assumes that variants are listed after the original ref, and before any new original ref transcripts
      if (name.find("_refsnp_") != std::string::npos) {
          std::string tname = name.substr(0, name.find("_refsnp_")); // TODO: need to make this safe against underscores :(
          std::string variant = name.substr(name.find("_refsnp_"), name.size());
          std::vector<std::string>::iterator it
                = std::find(target_names_.begin() + last_full_tr, target_names_.end(), tname);
          int index = std::distance(target_names_.begin(), it);
          // Update this shade to point to the original transcript
          shadeToColorTranscriptMap[target_names_.size() - 1] = index;
      } else{
          last_full_tr = target_names_.size() - 1;
      }
    }
    gzclose(fp);
    std::cerr << "Found " << shadeToColorTranscriptMap.size() << " shades" << std::endl;
    fp=0;
  }
  time(&end);
  std::cerr << "Transcripts took " << (end - start) << " seconds" << std::endl;

  if (polyAcount > 0) {
    std::cerr << "[build] warning: clipped off poly-A tail (longer than 10)" << std::endl << "        from " << polyAcount << " target sequences" << std::endl;
  }

  
  if (countNonNucl > 0) {
    std::cerr << "[build] warning: replaced " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl << "        with pseudorandom nucleotides" << std::endl;
  }
  if (countUNuc > 0) {
    std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
  }
  
  num_trans = seqs.size();
  
  // for each target, create it's own equivalence class
  for (int i = 0; i < seqs.size(); i++ ) {
    std::vector<int> single(1,i);
    ecmap.push_back(single);
    ecmapinv.insert({single,i});
  }
  
  BuildDeBruijnGraph(opt, seqs);
  BuildEquivalenceClasses(opt, seqs);
  //BuildEdges(opt);

}

void KmerIndex::BuildDeBruijnGraph(const ProgramOptions& opt, const std::vector<std::string>& seqs) {
  

  std::cerr << "[build] counting k-mers ... "; std::cerr.flush();
  // gather all k-mers
  time_t starter, ender;
  time(&starter);
  for (int i = 0; i < seqs.size(); i++) {
    const char *s = seqs[i].c_str();
    KmerIterator kit(s),kit_end;
    for (; kit != kit_end; ++kit) {
      kmap.insert({kit->first.rep(), KmerEntry()}); // don't care about repeats
    }
  }
  time(&ender);
  std::cerr << "done." << std::endl;
  std::cerr << "This process took " << difftime(ender, starter) << " seconds." << std::endl;
  std::cerr << "[build] building target de Bruijn graph ... " << std::endl; std::cerr.flush();
  // find out how much we can skip ahead for each k-mer.
  for (auto& kv : kmap) {
    if (kv.second.contig == -1) {
      // ok we haven't processed the k-mer yet
      std::vector<Kmer> flist, blist;

      // iterate in forward direction
      Kmer km = kv.first;
      Kmer end = km;
      Kmer last = end;
      Kmer twin = km.twin();
      bool selfLoop = false;
      flist.push_back(km);

      while (fwStep(end,end)) {
        if (end == km) {
          // selfloop
          selfLoop = true;
          break;
        } else if (end == twin) {
          selfLoop = (flist.size() > 1); // hairpins are not loops
          // mobius loop
          break;
        } else if (end == last.twin()) {
          // hairpin
          break;
        }
        flist.push_back(end);
        last = end;
      }

      Kmer front = twin;
      Kmer first = front;

      if (!selfLoop) {
        while (fwStep(front,front)) {
          if (front == twin) {
            // selfloop
            selfLoop = true;
            break;
          } else if (front == km) {
            // mobius loop
            selfLoop = true;
            break;
          } else if (front == first.twin()) {
            // hairpin
            break;
          }
          blist.push_back(front);
          first = front;
        }
      }

      std::vector<Kmer> klist;
      for (auto it = blist.rbegin(); it != blist.rend(); ++it) {
        klist.push_back(it->twin());
      }
      for (auto x : flist) {
        klist.push_back(x);
      }


      Contig contig;
      contig.id = dbGraph.contigs.size();
      contig.length = klist.size();
      contig.seq = klist[0].toString();
      contig.seq.reserve(contig.length + k-1);

      bool forward;
      for (int i = 0; i < klist.size(); i++) {
        Kmer x = klist[i];
        Kmer xr = x.rep();
        forward = (x==xr);
        auto it = kmap.find(xr);
        assert(it->second.contig==-1);
        it->second = KmerEntry(contig.id, contig.length, i, forward);
        if (i > 0) {
          contig.seq.push_back(x.toString()[k-1]);
        }
      }
      
      dbGraph.contigs.push_back(contig);
      dbGraph.ecs.push_back(-1);
    }
  }
  std::cerr << " done " << std::endl;

}

bool contigHasTrid(std::vector<ContigToTranscript>& trans, int trid) {
  bool ret = false;
  for (auto& ctt : trans) {
    if (ctt.trid == trid) {
      ret = true;
    }
  }
  return ret;
}

void KmerIndex::BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs) {
  std::cerr << "[build] creating equivalence classes ... "; std::cerr.flush();

  std::vector<std::vector<TRInfo>> trinfos(dbGraph.contigs.size());
  for (int i = 0; i < seqs.size(); i++) {
    int seqlen = seqs[i].size() - k + 1; // number of k-mers
    const char *s = seqs[i].c_str();
    KmerIterator kit(s), kit_end;
    for (; kit != kit_end; ++kit) {
      Kmer x = kit->first;
      Kmer xr = x.rep();
      auto search = kmap.find(xr);
      bool forward = (x==xr);
      KmerEntry val = search->second;
      std::vector<TRInfo>& trinfo = trinfos[val.contig];
      Contig& contig = dbGraph.contigs[val.contig];

      

      TRInfo tr;
      tr.trid = i;
      int jump = kit->second;
      if (forward == val.isFw()) {
        tr.sense = true;
        tr.start = val.getPos();
        if (contig.length - tr.start > seqlen - kit->second) {
          // tartget stops
          tr.stop = tr.start + seqlen - kit->second;
          jump = seqlen;
        } else {
          tr.stop = contig.length;
          jump = kit->second + (tr.stop - tr.start)-1;
        }
      } else {
        tr.sense = false;
        tr.stop = val.getPos()+1;
        int stpos = tr.stop - (seqlen - kit->second);
        if (stpos > 0) {
          tr.start = stpos;
          jump = seqlen;
        } else {
          tr.start = 0;
          jump = kit->second + (tr.stop - tr.start) - 1;
        }
      }

      // // debugging -->
      //std::cout << "covering seq" << std::endl << seqs[i].substr(kit->second, jump-kit->second +k) << std::endl;
      //std::cout << "id = " << tr.trid << ", (" << tr.start << ", " << tr.stop << ")" << std::endl;
      //std::cout << "contig seq" << std::endl;
      // if (forward == val.isFw()) {
      //   //std::cout << contig.seq << std::endl;
      //   assert(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start) == seqs[i].substr(kit->second, jump-kit->second +k) );
      // } else {
      //   //std::cout << revcomp(contig.seq) << std::endl;
      //   assert(revcomp(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start)) == seqs[i].substr(kit->second, jump-kit->second +k));
      // }
      // if (jump == seqlen) {
      //   //std::cout << std::string(k-1+(tr.stop-tr.start)-1,'-') << "^" << std::endl;
      // }

      // // <-- debugging
      
      trinfo.push_back(tr);
      kit.jumpTo(jump);
    }
  }

  
  FixSplitContigs(opt, trinfos);

  
  int perftr = 0;
  for (int i = 0; i < trinfos.size(); i++) {
    bool all = true;

    int contigLen = dbGraph.contigs[i].length;
    for (auto x : trinfos[i]) {
      if (x.start!=0 || x.stop !=contigLen) {
        all = false;
      }
    }
    if (all) {
      perftr++;
    } 
  }

  assert(dbGraph.contigs.size() == trinfos.size());
  // for each contig
  for (int i = 0; i < trinfos.size(); i++) {
    std::vector<int> u;
    // push back transcript id for all transcripts in this contig
    for (auto x : trinfos[i]) {
      u.push_back(x.trid);

      if (shadeToColorTranscriptMap.find(x.trid) != shadeToColorTranscriptMap.end()) {
          u.push_back(shadeToColorTranscriptMap[x.trid]);
      }
    }

    sort(u.begin(), u.end());
    if (!isUnique(u)) {
      std::vector<int> v = unique(u);
      swap(u,v);
    }

    assert(!u.empty());

    auto search = ecmapinv.find(u);
    int ec = -1;
    if (search != ecmapinv.end()) {
      // insert contig -> ec info
      ec = search->second;
    } else {
      ec = ecmapinv.size();
      ecmapinv.insert({u,ec});
      ecmap.push_back(u);
    }
    dbGraph.ecs[i] = ec;
    assert(ec != -1);
    
    // record the transc
    Contig& contig = dbGraph.contigs[i];
    contig.ec = ec;
    // correct ec of all k-mers in contig
  }

  // map transcripts to contigs
  //std::cout << std::endl;
  for (int i = 0; i < seqs.size(); i++) {
    int seqlen = seqs[i].size() - k + 1; // number of k-mers
    // debugging
    std::string stmp;
    const char *s = seqs[i].c_str();
    ////std::cout << "sequence number " << i << std::endl;
    //std::cout << ">" << target_names_[i] << std::endl;
    //std::cout << seqs[i] << std::endl;
    KmerIterator kit(s), kit_end;
    for (; kit != kit_end; ++kit) {
      Kmer x = kit->first;
      //std::cout << "position = " << kit->second << ", mapping " << x.toString() << std::endl;
      Kmer xr = x.rep();
      auto search = kmap.find(xr);
      bool forward = (x==xr);
      KmerEntry val = search->second;
      Contig& contig = dbGraph.contigs[val.contig];

      ContigToTranscript info;
      info.trid = i;
      info.pos = kit->second;
      info.sense = (forward == val.isFw());
      int jump = kit->second + contig.length-1;
      //std::cout << "mapped to contig " << val.contig << ", len = " << contig.length <<  ", pos = " << val.getPos() << ", sense = " << info.sense << std::endl;
      contig.transcripts.push_back(info);
      // debugging
      if (info.sense) {

        if (info.pos == 0) {
          stmp.append(contig.seq);
          //std::cout << contig.seq << std::endl;
        } else {
          stmp.append(contig.seq.substr(k-1));
          //std::cout << contig.seq.substr(k-1) << std::endl;
        }
      } else {
        std::string r = revcomp(contig.seq);
        if (info.pos == 0) {
          stmp.append(r);
          //std::cout << r << std::endl;
        } else {
          stmp.append(r.substr(k-1));
          //std::cout << r.substr(k-1) << std::endl;
        }
      }
      //std::cout << stmp << std::endl;
      //std::cout << "covering seq" << std::endl << seqs[i].substr(kit->second, jump-kit->second +k) << std::endl;
      //std::cout << "id = " << tr.trid << ", (" << tr.start << ", " << tr.stop << ")" << std::endl;
      //std::cout << "contig seq" << std::endl;
      // if (forward == val.isFw()) {
      //   //std::cout << contig.seq << std::endl;
      //   assert(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start) == seqs[i].substr(kit->second, jump-kit->second +k) );
      // } else {
      //   //std::cout << revcomp(contig.seq) << std::endl;
      //   assert(revcomp(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start)) == seqs[i].substr(kit->second, jump-kit->second +k));
      // }
      // if (jump == seqlen) {
      //   //std::cout << std::string(k-1+(tr.stop-tr.start)-1,'-') << "^" << std::endl;
      // }

      // // <-- debugging
      
      kit.jumpTo(jump);
    }
    if (seqlen > 0 && seqs[i] != stmp) {
      /*std::cout << ">" << target_names_[i] << std::endl
                << seqs[i] << std::endl
                << stmp << std::endl;*/
      assert(false);
    }
  }

  // double check the contigs
  // for every contig
  /* std::cerr << std::endl; */
  for (auto &c : dbGraph.contigs) {
    // for every transcript that is in this contig
    int old_size = c.transcripts.size();
    /* std::cerr << c.id << ": "; */
    std::string r;
    for (int i = 0; i < old_size; i++) {
      auto info = c.transcripts[i];
      if (info.sense) {
        r = c.seq;
      } else {
        r = revcomp(c.seq);
      }
      assert(r == seqs[info.trid].substr(info.pos,r.size()));
      // if this contig is in an alternative shade, add the trid corresponding to the original ref
      if (shadeToColorTranscriptMap.find(info.trid) != shadeToColorTranscriptMap.end()) { // &&
        /* if (!contigHasTrid(c.transcripts, shadeToColorTranscriptMap.at(info.trid))) { */
        /* std::cerr << target_names_[info.trid] << std::endl; */
        int color = shadeToColorTranscriptMap.at(info.trid);

        ContigToTranscript altContig;
        altContig.trid = color; // shadeToColorTranscriptMap.at(info.trid);
        // position of alternative is always 0, since it's its own contig
        // Not actually the case! can have a contig split an alt contig down the middle
        altContig.sense = info.sense;
        /* if (std::find_if( */
        /*       c.transcripts.begin(), */
        /*       c.transcripts.end(), */
        /*       [&ac = altContig] */
        /*       (const ContigToTranscript& ct) -> bool { return ac.trid == ct.trid; }) == c.transcripts.end()) { */
        std::vector<int> refsnp_vals;
        std::vector<int> altsnp_vals;
        std::vector<int> alt2snp_vals;
        int start_pos;
        convert_ornament_to_snps(
            target_names_[info.trid],
            refsnp_vals,
            altsnp_vals,
            alt2snp_vals,
            &start_pos);
        altContig.pos = start_pos;
        if (altsnp_vals.size() != 0 || alt2snp_vals.size() != 0) {
          c.transcripts.push_back(altContig);
        }
        /* } */
      }
      /* } */
    }
    /* std::cerr << std::endl << r << std::endl << std::endl; */
  }

  std::cerr << " done" << std::endl;
  std::cerr << "[build] target de Bruijn graph has " << dbGraph.contigs.size() << " contigs and contains "  << kmap.size() << " k-mers " << std::endl;
}

std::vector<std::string> KmerIndex::getTranscriptNamesForKmer(const char* kmer) {
    KmerIterator kit(kmer), kit_end;
    std::vector<std::string> transcript_names;
    if (kit == kit_end) {
        return transcript_names;
    }

    Kmer x = kit->first;
    Kmer xr = x.rep();
    bool forward = (x==xr);

    auto search = kmap.find(xr);
    if (search == kmap.end()) {
        return transcript_names;
    }

    KmerEntry val = search->second;
    Contig& contig = dbGraph.contigs[val.contig];

    // By definition, all the kmers in a contig have the exact same set of transcripts
    std::vector<int> trids;
    for (auto t : contig.transcripts) {
        trids.push_back(t.trid);
    }

    for (int trid : trids) {
        transcript_names.push_back(target_names_[trid]);
    }

    return transcript_names;
}

void KmerIndex::FixSplitContigs(const ProgramOptions& opt, std::vector<std::vector<TRInfo>>& trinfos) {

  int perftr = 0;
  int orig_size = trinfos.size();

  for (int i = 0; i < orig_size; i++) {
    bool all = true;

    int contigLen = dbGraph.contigs[i].length;
    for (auto x : trinfos[i]) {
      if (x.start!=0 || x.stop !=contigLen) {
        all = false;
      }
      assert(x.start < x.stop);
    }

    if (all) {
      perftr++;
    } else {
      // break up equivalence classes
      // sort by start/stop
      std::vector<int> brpoints;
      for (auto& x : trinfos[i]) {
        brpoints.push_back(x.start);
        brpoints.push_back(x.stop);
      }
      sort(brpoints.begin(), brpoints.end());
      assert(brpoints[0] == 0);
      assert(brpoints[brpoints.size()-1]==contigLen);

      // find unique points
      if (!isUnique(brpoints)) {
        std::vector<int> u = unique(brpoints);
        swap(u,brpoints);
      }

      assert(!brpoints.empty());
      
      // copy sequence
      std::string seq = dbGraph.contigs[i].seq;
      // copy old trinfo
      std::vector<TRInfo> oldtrinfo = trinfos[i];
      
      for (int j = 1; j < brpoints.size(); j++) {
        assert(brpoints[j-1] < brpoints[j]);
        Contig newc;
        newc.seq = seq.substr(brpoints[j-1], brpoints[j]-brpoints[j-1]+k-1);
        newc.length = brpoints[j]-brpoints[j-1];

        if (j>1) {
          newc.id = dbGraph.contigs.size();
          dbGraph.contigs.push_back(newc);
          dbGraph.ecs.push_back(-1);
        } else {
          newc.id = i;
          dbGraph.contigs[i] = newc;
        }

        // repair k-mer mapping
        KmerIterator kit(newc.seq.c_str()), kit_end;
        for (; kit != kit_end; ++kit) {
          Kmer x = kit->first;
          Kmer xr = x.rep();
          auto search = kmap.find(xr);
          assert(search != kmap.end());
          bool forward = (x==xr);
          search->second = KmerEntry(newc.id, newc.length,  kit->second, forward);
        }

        // repair tr-info
        std::vector<TRInfo> newtrinfo;
        for (auto x : oldtrinfo) {
          if (!(x.stop <= brpoints[j-1] || x.start >= brpoints[j])) {
            TRInfo trinfo;
            trinfo.sense = x.sense;
            trinfo.trid = x.trid;
            trinfo.start = 0;
            trinfo.stop = newc.length;
            newtrinfo.push_back(trinfo);
          }
        }
        if (j>1) {
          trinfos.push_back(newtrinfo);
        } else {
          trinfos[i] = newtrinfo;
        }
      }
    }
  }
}




void KmerIndex::write(const std::string& index_out, bool writeKmerTable) {
  std::cerr << index_out << std::endl;
  std::ofstream out;
  out.open(index_out, std::ios::out | std::ios::binary);

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }

  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. write k
  out.write((char *)&k, sizeof(k));

  // 3. write number of targets
  out.write((char *)&num_trans, sizeof(num_trans));

  // 4. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  size_t kmap_size = kmap.size();

  if (writeKmerTable) {
    // 5. write number of k-mers in map
    out.write((char *)&kmap_size, sizeof(kmap_size));

    // 6. write kmer->ec values
    for (auto& kv : kmap) {
      out.write((char *)&kv.first, sizeof(kv.first));
      out.write((char *)&kv.second, sizeof(kv.second));
    }
  } else {
    // 5. write fake k-mer size
    kmap_size = 0;
    out.write((char *)&kmap_size, sizeof(kmap_size));

    // 6. write none of the kmer->ec values
  }
  // 7. write number of equivalence classes
  size_t tmp_size;
  tmp_size = ecmap.size();
  out.write((char *)&tmp_size, sizeof(tmp_size));

  // 8. write out each equiv class
  //  for (auto& kv : ecmap) {
  for (int ec = 0; ec < ecmap.size(); ec++) {
    out.write((char *)&ec, sizeof(ec));
    auto& v = ecmap[ec];
    // 8.1 write out the size of equiv class
    tmp_size = v.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));
    // 8.2 write each member
    for (auto& val: v) {
      out.write((char *)&val, sizeof(val));
    }
  }

  // 9. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 9.1 write out how many bytes
    // XXX: Note: this doesn't actually encore the max targ id size.
    // might cause problems in the future
    // tmp_size = tid.size();
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 9.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }

  // 10. write out contigs
  if (writeKmerTable) {
    assert(dbGraph.contigs.size() == dbGraph.ecs.size());
    tmp_size = dbGraph.contigs.size();
    out.write((char*)&tmp_size, sizeof(tmp_size));
    for (auto& contig : dbGraph.contigs) {
      out.write((char*)&contig.id, sizeof(contig.id));
      out.write((char*)&contig.length, sizeof(contig.length));
      tmp_size = strlen(contig.seq.c_str());
      out.write((char*)&tmp_size, sizeof(tmp_size));
      out.write(contig.seq.c_str(), tmp_size);

      // 10.1 write out transcript info
      tmp_size = contig.transcripts.size();
      out.write((char*)&tmp_size, sizeof(tmp_size));
      for (auto& info : contig.transcripts) {
        out.write((char*)&info.trid, sizeof(info.trid));
        out.write((char*)&info.pos, sizeof(info.pos));
        out.write((char*)&info.sense, sizeof(info.sense));
      }
    }
    
    // 11. write out ecs info
    for (auto ec : dbGraph.ecs) {
      out.write((char*)&ec, sizeof(ec));
    }


  } else {
    // write empty dBG
    tmp_size = 0;
    out.write((char*)&tmp_size, sizeof(tmp_size));
  }

  // 12. write out variant information
  out.flush();
  out.close();
}

bool KmerIndex::fwStep(Kmer km, Kmer& end) const {
  int j = -1;
  int fw_count = 0;
  for (int i = 0; i < 4; i++) {
    Kmer fw_rep = end.forwardBase(Dna(i)).rep();
    auto search = kmap.find(fw_rep);
    if (search != kmap.end()) {
      j = i;
      ++fw_count;
      if (fw_count > 1) {
        return false;
      }
    }
  }

  if (fw_count != 1) {
    return false;
  }

  Kmer fw = end.forwardBase(Dna(j));

  int bw_count = 0;
  for (int i = 0; i < 4; i++) {
    Kmer bw_rep = fw.backwardBase(Dna(i)).rep();
    if (kmap.find(bw_rep) != kmap.end()) {
      ++bw_count;
      if (bw_count > 1) {
        return false;
      }
    }
  }

  if (bw_count != 1) {
    return false;
  } else {
    if (fw != km) {
      end = fw;
      return true;
    } else {
      return false;
    }
  }

}

void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable) {
  std::string& index_in = opt.index;
  std::ifstream in;


  in.open(index_in, std::ios::in | std::ios::binary);

  if (!in.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index input file could not be opened!";
    exit(1);
  }

  // 1. read version
  size_t header_version = 0;
  in.read((char *)&header_version, sizeof(header_version));

  if (header_version != INDEX_VERSION) {
    std::cerr << "Error: incompatible indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate";
    exit(1);
  }

  // 2. read k
  in.read((char *)&k, sizeof(k));
  if (Kmer::k == 0) {
    //std::cerr << "[index] no k has been set, setting k = " << k << std::endl;
    Kmer::set_k(k);
    opt.k = k;
  } else if (Kmer::k == k) {
    //std::cerr << "[index] Kmer::k has been set and matches" << k << std::endl;
    opt.k = k;
  } else {
    std::cerr << "Error: Kmer::k was already set to = " << Kmer::k << std::endl
              << "       conflicts with value of k  = " << k << std::endl;
    exit(1);
  }

  // 3. read in number of targets
  in.read((char *)&num_trans, sizeof(num_trans));

  // 4. read in length of targets
  target_lens_.clear();
  target_left_lens_.clear();
  target_right_lens_.clear();
  target_lens_.reserve(num_trans);
  target_left_lens_.reserve(num_trans);
  target_right_lens_.reserve(num_trans);

  for (int i = 0; i < num_trans; i++) {
    int tlen;
    in.read((char *)&tlen, sizeof(tlen));
    target_lens_.push_back(tlen);
    target_left_lens_.push_back(tlen);
    target_right_lens_.push_back(tlen);
  }

  // 5. read number of k-mers
  size_t kmap_size;
  in.read((char *)&kmap_size, sizeof(kmap_size));

  std::cerr << "[index] k-mer length: " << k << std::endl;
  std::cerr << "[index] number of targets: " << pretty_num(num_trans)
    << std::endl;
  std::cerr << "[index] number of k-mers: " << pretty_num(kmap_size)
    << std::endl;

  kmap.clear();
  if (loadKmerTable) {
    kmap.reserve(kmap_size,true);
  }

  // 6. read kmer->ec values
  Kmer tmp_kmer;
  KmerEntry tmp_val;
  for (size_t i = 0; i < kmap_size; ++i) {
    in.read((char *)&tmp_kmer, sizeof(tmp_kmer));
    in.read((char *)&tmp_val, sizeof(tmp_val));

    if (loadKmerTable) {
      kmap.insert({tmp_kmer, tmp_val});
    }
  }

  // 7. read number of equivalence classes
  size_t ecmap_size;
  in.read((char *)&ecmap_size, sizeof(ecmap_size));

  std::cerr << "[index] number of equivalence classes: "
    << pretty_num(ecmap_size) << std::endl;
  ecmap.resize(ecmap_size);
  int tmp_id;
  int tmp_ecval;
  size_t vec_size;
  // 8. read each equiv class
  for (size_t ec = 0; ec < ecmap_size; ++ec) {
    in.read((char *)&tmp_id, sizeof(tmp_id));

    // 8.1 read size of equiv class
    in.read((char *)&vec_size, sizeof(vec_size));

    // 8.2 read each member
    std::vector<int> tmp_vec;
    tmp_vec.reserve(vec_size);
    for (size_t j = 0; j < vec_size; ++j ) {
      in.read((char *)&tmp_ecval, sizeof(tmp_ecval));
      tmp_vec.push_back(tmp_ecval);
    }
    //ecmap.insert({tmp_id, tmp_vec});
    ecmap[tmp_id] = tmp_vec;
    ecmapinv.insert({tmp_vec, tmp_id});
  }

  // 9. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  size_t tmp_size;
  size_t bufsz = 1024;
  char *buffer = new char[bufsz];
  int last_full_tr = 0;
  for (auto i = 0; i < num_trans; ++i) {
    // 9.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size +1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }
    
    // clear the buffer 
    memset(buffer,0,bufsz);
    // 9.2 read in the character string
    in.read(buffer, tmp_size);

    /* std::string tmp_targ_id( buffer ); */
    std::string name = std::string(buffer);
    target_names_.push_back(name);
//    std::cerr << std::endl;
      // Assumes that variants are listed after the original ref, and before any new original ref transcripts
    if (name.find("_refsnp_") != std::string::npos) {
      std::string tname = name.substr(0, name.find("_refsnp_"));
      std::string variant = name.substr(name.find("_refsnp_"), name.size());
      std::vector<std::string>::iterator it
              = std::find(target_names_.begin() + last_full_tr, target_names_.end(), tname);
      int index = std::distance(target_names_.begin(), it);
      // Update this shade to point to the original transcript
      shadeToColorTranscriptMap[target_names_.size() - 1] = index;
    } else{
      last_full_tr = target_names_.size() - 1;
    }
  }

  
  // 10. read contigs
  size_t contig_size;
  in.read((char *)&contig_size, sizeof(contig_size));
  dbGraph.contigs.clear();
  dbGraph.contigs.reserve(contig_size);
  for (auto i = 0; i < contig_size; i++) {
    Contig c;
    in.read((char *)&c.id, sizeof(c.id));
    in.read((char *)&c.length, sizeof(c.length));
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size + 1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    memset(buffer,0,bufsz);
    in.read(buffer, tmp_size);
    c.seq = std::string(buffer); // copy
    
    // 10.1 read transcript info
    in.read((char*)&tmp_size, sizeof(tmp_size));
    c.transcripts.clear();
    c.transcripts.reserve(tmp_size);

    for (auto j = 0; j < tmp_size; j++) {
      ContigToTranscript info;
      in.read((char*)&info.trid, sizeof(info.trid));
      /* If this trid is a shade, and if this shade would not be in the underlying genotype, */
      /* then drop remove the base transcript from this contig */
      in.read((char*)&info.pos, sizeof(info.pos));
      in.read((char*)&info.sense, sizeof(info.sense));
      c.transcripts.push_back(info);
    }

    dbGraph.contigs.push_back(c);
  }

  // 11. read ecs info
  dbGraph.ecs.clear();
  dbGraph.ecs.reserve(contig_size);
  int tmp_ec;
  for (auto i = 0; i < contig_size; i++) {
    in.read((char *)&tmp_ec, sizeof(tmp_ec));
    dbGraph.ecs.push_back(tmp_ec);
  }


  // delete the buffer
  delete[] buffer;
  buffer=nullptr;
  
  in.close();

  std::cerr << opt.phased << std::endl;
  time_t start, end;
  time(&start);
  initVariantParamEntry(opt.vcf_file, opt.sample_name, opt.phased);
  time(&end);
  std::cerr << "Initializing variant parameter entries took " << difftime(end, start) << " seconds." << std::endl;
  time(&start);
  filterVariants(opt, opt.phased);
  time(&end);
  std::cerr << "Filtering and pruning DBG took " << difftime(end, start) << " seconds." << std::endl;
}

void KmerIndex::convert_string_to_arr(std::string variant_info, std::vector<int>& arr) const {
  std::stringstream s_stream(variant_info);
  std::string curr;
  while(std::getline(s_stream, curr, ',')) {
    arr.push_back(std::stoi(curr));
  }
}

int KmerIndex::get_ornament_position(int shade) const {
  std::string ornament_name = target_names_[shade];
  std::string tname = ornament_name.substr(0, ornament_name.find("_refsnp_"));
  std::string variant = ornament_name.substr(ornament_name.find("_refsnp_") + 1, ornament_name.size());
  auto second_split = variant.find_last_of("|");
  return std::stoi(variant.substr(second_split + 1, variant.size()));
}

void KmerIndex::convert_ornament_to_snps(
    std::string ornament_name,
    std::vector<int>& refsnp_vals,
    std::vector<int>& altsnp_vals,
    std::vector<int>& altsnp2_vals,
    int* start_pos) const {
  assert(ornament_name.find("snp") != std::string::npos);
  std::string tname = ornament_name.substr(0, ornament_name.find("_refsnp_"));
  std::string variant = ornament_name.substr(ornament_name.find("_refsnp_") + 1, ornament_name.size());
  auto split = variant.find("|");
  auto second_split = variant.find_last_of("|");
  std::string refsnp = variant.substr(0, split);
  std::string altsnp = variant.substr(split + 1, second_split - split - 1);
  *start_pos = std::stoi(variant.substr(second_split + 1, variant.size()));
  refsnp.pop_back();
  altsnp.pop_back();

  convert_string_to_arr(refsnp.substr(refsnp.find("(") + 1, refsnp.size()), refsnp_vals);

  int alt0_start = altsnp.find("(") + 1;
  int alt0_end = altsnp.find(")");
  std::string alts = altsnp.substr(alt0_start, alt0_end - alt0_start);
  convert_string_to_arr(alts, altsnp_vals);

  int alt2pos = altsnp.find("altsnp1_");
  if (alt2pos != std::string::npos) {
    int alt2end = second_split - 1;
    std::string alt2s = altsnp.substr(alt2pos + 9, alt2end - (alt2pos + 1));
    convert_string_to_arr(alt2s, altsnp2_vals);
  }
}

void KmerIndex::filterVariants(const ProgramOptions& opt, bool filter_phased) {
  std::unordered_map<int, int> spurious_ornaments;
  std::unordered_set<int> spurious_shades;
  std::unordered_set<int> spurious_colors;
  std::vector<int> refsnpVals;
  std::vector<int> altsnpVals;
  std::vector<int> alt2snpsVals;
  int start_pos;
  int num_pruned = 0;
  for (auto& kv : shadeToColorTranscriptMap) {
    std::string name = target_names_[kv.first];
    refsnpVals.clear();
    altsnpVals.clear();
    alt2snpsVals.clear();
    convert_ornament_to_snps(name, refsnpVals, altsnpVals, alt2snpsVals, &start_pos);

    /* Verify that every one of the above snps is consistent with this individual's homozygous genotypes */
    /* For shades which are not consistent, go to every contig that has this segment and remove the base trid from it */
    int color_id = kv.second;
    
    bool retain = true;
    bool color_retain = false;
    bool bc_retain = true;
    int num_het = 0;

    if (transcriptToParam.find(color_id) == transcriptToParam.end()) {
      /* If this transcript is not in the param map, means that individual was all homo */
      auto homo_alts = homo_alt_snps[color_id];
      for (int rsnp : refsnpVals) {
        if (homo_alts.find(rsnp) != homo_alts.end()) {
          retain = false;
        }
      }
      
      for (int asnp : altsnpVals) {
        if (homo_alts.find(asnp) == homo_alts.end()) {
          retain = false;
        }
      }

      bc_retain = false;
    } else {
      /* If there is a parameter for this position, then we are het and so we are good */
      auto& trans_snp_map = transcriptToParam.at(color_id);
      auto homo_alts = homo_alt_snps[color_id];
      /* This transcript is in the map, meaning we have some het snp in this transcript */
      /* Go through all rsnps, if we are homo alt at any of those pos then discard shade */
      /* Go through all asnps, if we are homo ref at any of those pos then discard shade */
      for (int rsnp : refsnpVals) {
        // This ornament contains a reference shade. If the sample is homozygous reference
        // or heterozygous, then it's fine to keep this ornament. If it's homozygous alternative,
        // then we need to remove the shade.
        if (homo_alts.find(rsnp) != homo_alts.end()) {
          retain = false;
        }
      }

      // This ornament contains an alternative shade. If the sample is homozygous reference, then we
      // need to remove the shade. A sample is homozygous reference if: 1. it is not homozygous alt,
      // and 2. it is not heterozygous
      for (int asnp : altsnpVals) {
        if (trans_snp_map.find(asnp) == trans_snp_map.end() &&
            homo_alts.find(asnp) == homo_alts.end()) {
          retain = false;
        }
      }

      /* If we are running in phased mode, perform an additional check to make sure that this ornament adheres to the phasing */
      /* char phase_direction = 0; */
      /* for (int rsnp : refsnpVals) { */
      /*   auto s = trans_snp_map.find(rsnp); */
      /*   if (s != trans_snp_map.end()) { */
      /*     num_het++; */
      /*     if (filter_phased) { */
      /*       int param_index = s->second; */
      /*       auto& vpe_list = effectiveASEParams.at(param_index); */
      /*       for (auto& vpe : vpe_list) { */
      /*         /1* For the current heterozygous snp entry... *1/ */
      /*         if (vpe.trid == kv.second && vpe.pos == rsnp) { */
      /*           int curr_phase = vpe.orientation == 'r' ? 1 : -1; */
      /*           if (phase_direction == 0) { */
      /*             phase_direction = curr_phase; */
      /*           } else { */
      /*             if (phase_direction != curr_phase) retain = false; */
      /*           } */
      /*           break; */
      /*         } else { */
      /*           continue; */
      /*         } */
      /*       } */
      /*     } */
      /*   } */
      /* } */

      /* /1* std::cerr  << name << ", " << retain << std::endl; *1/ */
      /* for (int asnp : altsnpVals) { */
      /*   auto s = trans_snp_map.find(asnp); */
      /*   if (s != trans_snp_map.end()) { */
      /*     num_het++; */
      /*     if (filter_phased) { */
      /*       int param_index = s->second; */
      /*       auto& vpe_list = effectiveASEParams.at(param_index); */
      /*       for (auto& vpe : vpe_list) { */
      /*         /1* For the current heterozygous snp entry... *1/ */
      /*         if (vpe.trid == kv.second && vpe.pos == asnp && vpe.is_first) { */
      /*           int curr_phase = vpe.orientation == 'a' ? 1 : -1; */
      /*           if (phase_direction == 0) { */
      /*             phase_direction = curr_phase; */
      /*           } else { */
      /*             if (phase_direction != curr_phase) retain = false; */
      /*           } */
      /*           break; */
      /*         } else { */
      /*           continue; */
      /*         } */
      /*       } */
      /*     } */
      /*   } */
      /* } */

      /* for (int asnp : alt2snpsVals) { */
      /*   auto s = trans_snp_map.find(asnp); */
      /*   if (s != trans_snp_map.end()) { */
      /*     num_het++; */
      /*     if (filter_phased) { */
      /*       int param_index = s->second; */
      /*       auto& vpe_list = effectiveASEParams.at(param_index); */
      /*       for (auto& vpe : vpe_list) { */
      /*         /1* For the current heterozygous snp entry... *1/ */
      /*         if (vpe.trid == kv.second && vpe.pos == asnp && !vpe.is_first) { */
      /*           int curr_phase = vpe.orientation == 'a' ? 1 : -1; */
      /*           if (phase_direction == 0) { */
      /*             phase_direction = curr_phase; */
      /*           } else { */
      /*             if (phase_direction != curr_phase) retain = false; */
      /*           } */
      /*           break; */
      /*         } else { */
      /*           continue; */
      /*         } */
      /*       } */
      /*     } */
      /*   } */
      /* } */

      /* Lastly, if all snps in the ornament are homozygous for this individual, then remove even if consistent */
    }

    // 

    /* If retain is equal to false, but this ornament appears somewhere in what would have been */
    /* reference if this was just one sample, then remove the ornament rather than the color */
    /* if (bc_weights.find(kv.first) != bc_weights.end()) { */
    /*   if (!retain || num_het == 0) { */
    /*     color_retain = true; */
    /*     /1* Also, delete this ornament out of bc_weights, since we will no longer ever pick it up *1/ */
    /*     bc_weights.erase(kv.first); */
    /*   } */
    /* } */

    // We want to remove the ornament if retain is false.
    // This involves removing the contig if it only contains that ornament's color,
    // or only removing the shade (and its color) if another color exists there

    /* if (bc_weights.find(kv.first) != bc_weights.end()) { */
    /*   color_retain = true; */
    /* } else { */
    /*   if (!retain) color_retain = false; */
    /*   /1* if (!retain && num_het == 0) { *1/ */
    /*   /1*   color_retain = false; *1/ */
    /*   /1* } else if (!retain && num_het > 0) { *1/ */
    /*   /1*   /2* Over multiple SNPs within k bp, one het but genotype or phase of this ornament doesn't match *2/ *1/ */

    /*   /1* } *1/ */
    /* } */

    /* Else if retain */ 

    /* Currently doing this separately, but maybe we should remove all ornaments that wouldn't be here otherwise, even if the genotypes aren't consistent? */
    /* TODO: look into this */
    /* if (!bc_retain || num_het == 0) { */
    /*   bc_weights.erase(kv.first); */
    /* } */

    /* Actually maybe this is fine? */
    /* This is not fine, need to retain ornament over homozygous site? */
    /* This site would have been imputed into the genome for one sample. */
    /* So, delete shade if genotype matches, and delete color and shade both otherwise */
    /* if (!retain || num_het == 0) { */
    if (!retain) {
      num_pruned++;
      spurious_shades.insert(kv.first);
    }
  }

  for (auto &c : dbGraph.contigs) {
    // Remove the shade if it is spurious
    auto shade_it = std::find_if(
        c.transcripts.begin(),
        c.transcripts.end(),
        [&spurious_shades=spurious_shades](const ContigToTranscript& info) {
          return spurious_shades.find(info.trid) != spurious_shades.end();
        });
    if (shade_it != c.transcripts.end()) {
      int ind = shade_it - c.transcripts.begin();
      int color = shadeToColorTranscriptMap.at(c.transcripts[ind].trid);
      c.transcripts.erase(shade_it);

      // Remove ONE instance of the color if the shade was spurious
      auto color_it = std::find_if(
          c.transcripts.begin(),
          c.transcripts.end(),
          [&color=color](const ContigToTranscript& info) {
            return info.trid == color;
          });
      if (color_it != c.transcripts.end()) {
        c.transcripts.erase(color_it);
      }
    }

    auto& ec = ecmap[dbGraph.ecs[c.id]];
    auto shade_it_ec = std::find_if(
        ec.begin(),
        ec.end(),
        [&spurious_shades=spurious_shades](const int trid) {
          return spurious_shades.find(trid) != spurious_shades.end();
        });
    if (shade_it_ec != ec.end()) {
      int ind_ec = shade_it_ec - ec.begin();
      int color_ec = shadeToColorTranscriptMap.at(ec[ind_ec]);
      ec.erase(shade_it_ec);

      // Remove ONE instance of the color if the shade was spurious
      auto color_it_ec = std::find(ec.begin(), ec.end(), color_ec);
      if (color_it_ec != ec.end()) {
        ec.erase(color_it_ec);
      }
    }
  }
  std::cerr << "number of ornaments pruned: " << num_pruned << std::endl;
}

std::unordered_map<std::string, double> load_transcript_weights(std::string file_name, char delim) {
  std::ifstream is(file_name);
  std::string line;
  std::unordered_map<std::string, double> tw;
  while(std::getline(is, line)) {
    std::stringstream ss(line);
    std::string token;
    std::string tname;
    double tweight;

    // Read the respective parameters from the current line. Currently only reading first 10
    // TODO: Change effective params to be one indexed
    int i = 0;
    while (std::getline(ss, token, delim)) {
      switch (i) {
        case 0: // first item in line
          tname = token;
          break;
        case 1: // second item in line
          tweight = std::stod(token);
          break;
      }
      i += 1;
    }

    tw[tname] = tweight;
  }

  return tw;
}

// Setup variant param entry structure
void KmerIndex::initVariantParamEntry(std::string vcf_file, std::string sample, bool is_phased) {
  std::cerr << std::endl << "Reading vcf file for sample: " << vcf_file << std::endl;
  std::ifstream is(vcf_file);
  std::string line;
  bool skip = false;
  std::vector<VariantParameterEntry> currentParameterRun;
  bool prevWasPhased = false;
  std::string prevTranscriptName; // .empty checks if this has been set yet
  int prevTpos;
  effectiveASEParams[0] = {};
  std::unordered_map<std::string, int> nameToID;
  int ind = -1;
  /* Modify this to also include ornaments? So that we can set up bias correction */
  for (int i = 0; i < target_names_.size(); i++) {
    if (shadeToColorTranscriptMap.find(i) == shadeToColorTranscriptMap.end()) {
      nameToID[target_names_[i]] = i;
    }
  }

  while(std::getline(is, line)) {
    if (line.find("CHROM") != std::string::npos) {
      skip = true;
      // Get the index corresponding to the current sample
      std::string token;
      std::stringstream ss(line);
      int i = 0;
      while (std::getline(ss, token, '\t')) {
        trim(token);
        if (sample == token) {
          ind = i; 
          break;
        }
        i += 1;
      }

      if (ind == -1) {
        std::cerr << "ERROR: Specified sample was not found in specified VCF file." << std::endl;
        exit(1);
      } else {
        std::cerr << "Sample index found at position " << ind << std::endl;
      }
      continue;
    }

    if (!skip) { // skip all header lines
      continue;
    }

    std::stringstream ss(line);
    std::string token;
    std::string tname;
    int tpos;
    std::string ref;
    std::string alt;
    std::string genotype;

    // Read the respective parameters from the current line. Currently only reading first 10
    int i = 0;
    while (std::getline(ss, token, '\t')) {
      switch (i)
      {
        case 0:
          tname = token;
          break;
        case 1:
          tpos = std::stoi(token) - 1;
          break;
        case 3:
          ref = token;
          break;
        case 4:
          alt = token;
          break;
      }

      if (i == ind) {
        genotype = token;
        break;
      }
      i += 1;
    }

    bool currPhased;
    int genotype_l;
    int genotype_r;

    // skip sites that are male X, Y, or MT since there is no phasing to be done
    if (genotype.length() == 1) {
      if (nameToID.find(tname) != nameToID.end()) {
        int trid = nameToID[tname];
        if (genotype == "1") {
          homo_alt_snps[trid].insert(tpos);
        } 
      }
      continue;
    }

    int delim = genotype.find('|');
    if (delim != std::string::npos) {
      currPhased = is_phased;
      // TODO: Make this safe against triallelic sites
      genotype_l = std::stoi(genotype.substr(0, delim));
      genotype_r = std::stoi(genotype.substr(delim+1, genotype.size()));
    } else {
      delim = genotype.find('/');
      assert(delim !=  std::string::npos);
      currPhased = false;
      genotype_l = std::stoi(genotype.substr(0, delim));
      genotype_r = std::stoi(genotype.substr(delim+1, genotype.size()));
    }

    if (!prevTranscriptName.empty()) {
      if (prevTranscriptName != tname) {
//        std::cerr << "Switch " << prevTranscriptName << " to " << tname << std::endl;
        if (currentParameterRun.size() > 0) {
          if (nameToID.find(prevTranscriptName) != nameToID.end()) {
            effectiveASEParams[effectiveASEParams.size()] = currentParameterRun;
          }
          currentParameterRun.clear();
        }
        prevWasPhased = false;
      }
    }

    // We want to skip the homozygous sites for this individual
    if (genotype_l == genotype_r) {
      prevTranscriptName = tname;
      prevTpos = tpos;
      if (nameToID.find(tname) != nameToID.end()) {
        int trid = nameToID[tname];
        if (genotype_l == 1) {
          homo_alt_snps[trid].insert(tpos);
        }
      }
      continue;
    }

    // If we are here we know the current site is heterozygous
    // Want just an array of vectors, first corresponds to ref, next correspond to the alts
    if (prevWasPhased) {
      // If both are phased and we are here, they are part of the same transcript. Group the parameters
      if (currPhased) {
        // If multiallelic, treat the first occurring one as ref
        if (nameToID.find(tname) != nameToID.end()) {
          int trid = nameToID.at(tname);
          VariantParameterEntry vpe(
            trid,
            tpos,
            (genotype_l == 0 ? 'r' : 'a'),
            ref,
            alt
          );

          /* Check and see if this is a multiallelic site. Assumes sorted input file */
          if (currentParameterRun.size() > 0) {
            auto& prev_vpe = currentParameterRun[currentParameterRun.size() - 1];
            if (prev_vpe.pos == vpe.pos) {
              vpe.is_first = false;
            }
          }
          if (vpe.is_indel()) {
            if (genotype_l == 1) target_left_lens_[trid] += vpe.delta_length();
            if (genotype_r == 1) target_right_lens_[trid] += vpe.delta_length();
          }
          currentParameterRun.push_back(vpe);
        }
        prevWasPhased = true; // not needed
      } else {
        if (nameToID.find(tname) != nameToID.end()) {
          effectiveASEParams[effectiveASEParams.size()] = currentParameterRun;
        }
        currentParameterRun.clear();

        if (nameToID.find(tname) != nameToID.end()) {
          int trid = nameToID.at(tname);
          VariantParameterEntry vpe(
            trid,
            tpos,
            (genotype_l == 0 ? 'r' : 'a'),
            ref,
            alt
          );
          if (vpe.is_indel()) {
            if (genotype_l == 1) target_left_lens_[trid] += vpe.delta_length();
            if (genotype_r == 1) target_right_lens_[trid] += vpe.delta_length();
          }
          currentParameterRun.push_back(vpe);
          effectiveASEParams[effectiveASEParams.size()] = currentParameterRun;
        }
        currentParameterRun.clear();
        prevWasPhased = false;
      }
    } else {
      if (currPhased) {
        if (nameToID.find(tname) != nameToID.end()) {
          int trid = nameToID.at(tname);
          VariantParameterEntry vpe(
            trid,
            tpos,
            (genotype_l == 0 ? 'r' : 'a'),
            ref,
            alt
          );
          if (vpe.is_indel()) {
            if (genotype_l == 1) target_left_lens_[trid] += vpe.delta_length();
            if (genotype_r == 1) target_right_lens_[trid] += vpe.delta_length();
          }
          currentParameterRun.push_back(vpe);
        }
        prevWasPhased = true;
      } else {
        if (nameToID.find(tname) != nameToID.end()) {
          int trid = nameToID.at(tname);
          VariantParameterEntry vpe(
            trid,
            tpos,
            'r', // for unphased parameters, treat the parameter we are learning as REF! 
            ref,
            alt
          );

          if (effectiveASEParams[effectiveASEParams.size() - 1].size() > 0) {
            auto& prev_vpe = effectiveASEParams[effectiveASEParams.size() - 1][0];
            if (prev_vpe.trid == vpe.trid && prev_vpe.pos == vpe.pos) {
              vpe.is_first = false;
            }
          }
          // What to do here for multiallelic site?
          if (vpe.is_indel()) {
            if (genotype_l == 1) target_left_lens_[trid] += vpe.delta_length();
            if (genotype_r == 1) target_right_lens_[trid] += vpe.delta_length();
          }
          currentParameterRun.push_back(vpe);
          effectiveASEParams[effectiveASEParams.size()] = currentParameterRun;
        }
        currentParameterRun.clear();
        prevWasPhased = false;
      }
    }

    prevTranscriptName = tname;
    prevTpos = tpos;
  }

  // write the last parameter to the list of params
  if (currentParameterRun.size() > 0) {
    if (nameToID.find(prevTranscriptName) != nameToID.end()) {
      effectiveASEParams[effectiveASEParams.size()] = currentParameterRun;
    }
  }

  for (auto& it : effectiveASEParams) {
    for (VariantParameterEntry& vpe : it.second) {
      // trid corresponds to the color, pos corresponds to the shade
      // subtract one above because we use 0-indexing in the ornaments file
      if (vpe.trid != -1) {
        transcriptToParam[vpe.trid][vpe.pos] = it.first;
      }
    }
  }

  std::cerr << "[index] Effective parameters array is setup and has " << effectiveASEParams.size() - 1 << " entries" << std::endl;
  std::cerr << "[index] Transcript to parameter map is setup and has " << transcriptToParam.size() << " entries" << std::endl;
}

int KmerIndex::mapPair(const char *s1, int l1, const char *s2, int l2, int ec) const {
  bool d1 = true;
  bool d2 = true;
  int p1 = -1;
  int p2 = -1;
  int c1 = -1;
  int c2 = -1;


  KmerIterator kit1(s1), kit_end;

  bool found1 = false;
  for (; kit1 != kit_end; ++kit1) {
    Kmer x = kit1->first;
    Kmer xr = x.rep();
    auto search = kmap.find(xr);
    bool forward = (x==xr);

    if (search != kmap.end()) {
      found1 = true;
      KmerEntry val = search->second;
      c1 = val.contig;
      if (forward == val.isFw()) {
        p1 = val.getPos() - kit1->second;
        d1 = true;
      } else {
        p1 = val.getPos() + k + kit1->second;
        d1 = false;
      }
      break;
    }
  }

  if (!found1) {
    return -1;
  }

  

  KmerIterator kit2(s2);
  bool found2 = false;

  for (; kit2 != kit_end; ++kit2) {
    Kmer x = kit2->first;
    Kmer xr = x.rep();
    auto search = kmap.find(xr);
    bool forward = (x==xr);

    if (search != kmap.end()) {
      found2 = true;
      KmerEntry val = search->second;
      c2 = val.contig;
      if (forward== val.isFw()) {
        p2 = val.getPos() - kit2->second;
        d2 = true;
      } else {
        p2 = val.getPos() + k + kit2->second;
        d2 = false;
      }
      break;
    }
  }

  if (!found2) {
    return -1;
  }

  if (c1 != c2) {
    return -1;
  }

  if ((d1 && d2) || (!d1 && !d2)) {
    return -1;
  }

  if (p1>p2) {
    return p1-p2;
  } else {
    return p2-p1;
  }

}

// use:  match(s,l,v), s is string (read), l is length of read
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s
void KmerIndex::match(const char *s, int l, std::vector<std::pair<KmerEntry, int>>& v) const {
  KmerIterator kit(s), kit_end;
  bool backOff = false;
  int nextPos = 0; // nextPosition to check
//  std::cerr << std::endl;
  for (int i = 0;  kit != kit_end; ++i,++kit) {
    // need to check it
    auto search = kmap.find(kit->first.rep());
    int pos = kit->second;

    if (search != kmap.end()) {

      KmerEntry val = search->second;

      v.push_back({val, kit->second});

      // see if we can skip ahead
      // bring thisback later
      bool forward = (kit->first == search->first);
      int dist = val.getDist(forward);


      //const int lastbp = 10;
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) {
          // if we can jump beyond the read, check the end
          nextPos = l-k;
        }

        // check next position
        KmerIterator kit2(kit);
        kit2.jumpTo(nextPos);
        if (kit2 != kit_end) {
          Kmer rep2 = kit2->first.rep();
          auto search2 = kmap.find(rep2);
          bool found2 = false;
          int  found2pos = pos+dist;
          if (search2 == kmap.end()) {
            found2=true;
            found2pos = pos;
          } else if (val.contig == search2->second.contig) {
            found2=true;
            found2pos = pos+dist;
          }
          if (found2) {
            // great, a match (or nothing) see if we can move the k-mer forward
            if (found2pos >= l-k) {
              if (val.contig == search2->second.contig) {
                v.push_back({search2->second, l-k});
              } else {
                v.push_back({val, l-k});
              }
              break; //
            } else {
              if (val.contig == search2->second.contig) {
                v.push_back({search2->second, found2pos});
              } else {
                v.push_back({val, found2pos});
              }
              kit = kit2; // move iterator to this new position
            }
          } else {
            // this is weird, let's try the middle k-mer
            bool foundMiddle = false;
            if (dist > 4) {
              int middlePos = (pos + nextPos)/2;
              int middleContig = -1;
              int found3pos = pos+dist;
              KmerIterator kit3(kit);
              kit3.jumpTo(middlePos);
              KmerEntry val3;
              if (kit3 != kit_end) {
                Kmer rep3 = kit3->first.rep();
                auto search3 = kmap.find(rep3);
                if (search3 != kmap.end()) {
                  middleContig = search3->second.contig;
                  if (middleContig == val.contig) {
                    foundMiddle = true;
                    found3pos = middlePos;
                  } else if (middleContig == search2->second.contig) {
                    foundMiddle = true;
                    found3pos = pos+dist;
                  }
                }


                if (foundMiddle) {
                  v.push_back({search3->second, found3pos});
                  if (nextPos >= l-k) {
                    break;
                  } else {
                    kit = kit2; 
                  }
                }
              }
            }


            if (!foundMiddle) {
              ++kit;
              backOff = true;
              goto donejumping; // sue me Dijkstra!
            }
          }
        } else {
          // the sequence is messed up at this point, let's just take the match
          break;
        }
      }
    }

donejumping:

    if (backOff) {
      // backup plan, let's play it safe and search incrementally for the rest, until nextStop
      for (int j = 0; kit != kit_end; ++kit,++j) {
        if (j==skip) {
          j=0;
        }
        if (j==0) {
          // need to check it
          Kmer rep = kit->first.rep();
          auto search = kmap.find(rep);
          if (search != kmap.end()) {
            // if k-mer found
            v.push_back({search->second, kit->second}); // add equivalence class, and position
          }
        }

        if (kit->second >= nextPos) {
          backOff = false;
          break; // break out of backoff for loop
        }
      }
    }
  }
}

std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, int p) const {
  auto it = kmap.find(km.rep());
  if (it != kmap.end()) {
    KmerEntry val = it->second;
    return findPosition(tr, km, val, p);
  } else {
    return {-1,true};
  }
}

//use:  (pos,sense) = index.findPosition(tr,km,val,p)
//pre:  index.kmap[km] == val,
//      km is the p-th k-mer of a read
//      val.contig maps to tr
//post: km is found in position pos (1-based) on the sense/!sense strand of tr
std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, KmerEntry val, int p) const {
  bool fw = (km == km.rep());
  bool csense = (fw == val.isFw());

  int trpos = -1;
  bool trsense = true;
  if (val.contig < 0) {
    return {-1, true};
  }

  const Contig &c = dbGraph.contigs[val.contig];
  for (auto x : c.transcripts) {
    if (x.trid == tr) {
      trpos = x.pos;
      trsense = x.sense;
      break;
    }
  }

  if (trpos == -1) {
    return {-1,true};
  }



  if (trsense) {
    if (csense) {
      return {trpos + val.getPos() - p + 1, csense}; // 1-based, case I
    } else {
      return {trpos + val.getPos() + k + p, csense}; // 1-based, case III
    }
  } else {
    if (csense) {
      return {trpos + (c.length - val.getPos() -1) + k + p, !csense};  // 1-based, case IV
    } else {
      return {trpos + (c.length - val.getPos())  - p, !csense}; // 1-based, case II
    }
  }
}

std::pair<int, bool> KmerIndex::findFirstOrnamentAgnosticPosition(int tr, Kmer km, KmerEntry& val, int p) const {
  auto x = findPosition(tr, km, val, p);
  const Contig &c = dbGraph.contigs[val.contig];
  int start_pos = -1;
  int end_pos = -1;
  std::vector<int> refsnps;
  std::vector<int> altsnps;
  std::vector<int> alt2snps;
  for (auto t : c.transcripts) {
    auto se = shadeToColorTranscriptMap.find(t.trid);
    if (se != shadeToColorTranscriptMap.end() && se->second == tr) {
      convert_ornament_to_snps(target_names_[t.trid], refsnps, altsnps, alt2snps, &start_pos);
      end_pos = start_pos + target_lens_[t.trid];
    }
  }
  if (start_pos == -1 || (refsnps.size() > 0 && altsnps.size() == 0)) {
    return {x.first, x.second};
  } else {
    return {x.first + start_pos, x.second};
  }

  return {-1, true};
}

// use:  res = intersect(ec,v)
// pre:  ec is in ecmap, v is a vector of valid targets
//       v is sorted in increasing order
// post: res contains the intersection  of ecmap[ec] and v sorted increasing
//       res is empty if ec is not in ecma
std::vector<int> KmerIndex::intersect(int ec, const std::vector<int>& v) const {
  std::vector<int> res;
  if (ec < ecmap.size()) {
    auto& u = ecmap[ec];
    res.reserve(v.size());

    auto a = u.begin();
    auto b = v.begin();

    while (a != u.end() && b != v.end()) {
      if (*a < *b) {
        ++a;
      } else if (*b < *a) {
        ++b;
      } else {
        // match
        res.push_back(*a);
        ++a;
        ++b;
      }
    }
  }
  return res;
}


void KmerIndex::loadTranscriptSequences() const {
  if (target_seqs_loaded) {
    return;
  }
  
  std::vector<std::vector<std::pair<int, ContigToTranscript>>> trans_contigs(num_trans);
  for (auto &c : dbGraph.contigs) {
    for (auto &ct : c.transcripts) {
      trans_contigs[ct.trid].push_back({c.id, ct});
    }
  }

  auto &target_seqs = const_cast<std::vector<std::string>&>(target_seqs_);
  
  for (int i = 0; i < trans_contigs.size(); i++) {
    auto &v = trans_contigs[i];
    std::sort(v.begin(), v.end(), [](std::pair<int,ContigToTranscript> a, std::pair<int,ContigToTranscript> b) {
        return a.second.pos < b.second.pos;
      });

    std::string seq;
    seq.reserve(target_lens_[i]);
    for (auto &pct : v) {
      auto ct = pct.second;
      int start = (ct.pos==0) ? 0 : k-1;
      const auto& contig = dbGraph.contigs[pct.first];
      if (ct.sense) {
        seq.append(contig.seq.substr(start));
      } else {
        seq.append(revcomp(contig.seq).substr(start));
      }
    }
    target_seqs.push_back(seq);
  }

  bool &t = const_cast<bool&>(target_seqs_loaded);
  t = true;//target_seqs_loaded = true;
  return;
}

void KmerIndex::clear() {
  kmap.clear_table();
  ecmap.resize(0);
  dbGraph.ecs.resize(0);
  dbGraph.contigs.resize(0);
  {
    std::unordered_map<std::vector<int>, int, SortedVectorHasher> empty;
    std::swap(ecmapinv, empty);
  }
  
  target_lens_.resize(0);
  target_names_.resize(0);
  target_seqs_.resize(0);
}

void KmerIndex::writePseudoBamHeader(std::ostream &o) const {
  // write out header
  o << "@HD\tVN:1.0\n";
  for (int i = 0; i < num_trans; i++) {
    o << "@SQ\tSN:" << target_names_[i] << "\tLN:" << target_lens_[i] << "\n";
  }
  o << "@PG\tID:kallisto\tPN:kallisto\tVN:"<< KALLISTO_VERSION << "\n";
  o.flush();
}
