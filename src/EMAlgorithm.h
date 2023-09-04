#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "KmerIndex.h"
#include "MinCollector.h"
#include "weights.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <vector>

// smallest weight we expect is ~10^-4
// on most machines, TOLERANCE should be 2.22045e-15
const double TOLERANCE = std::numeric_limits<double>::denorm_min();
const double eps = 1e-08;

struct SortedIntPairHasher {
    size_t operator()(const std::pair<int, int>& p) const {
      uintmax_t hash = std::hash<int>{}(p.first);
      uintmax_t hash2 = std::hash<int>{}(p.second);
      hash ^= hash2 + 0x9e3779b9 + (hash<<6) + (hash>>2);
      return hash;
    }
};

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
  os<<"{ ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
  os<<"}";
  return os;
}

struct EMAlgorithm {
  EMAlgorithm(const std::vector<int>& counts,
              const std::unordered_map<int, int>& aseCounts,
              const KmerIndex& index,
              const MinCollector& tc,
              const std::vector<double>& all_means,
              const ProgramOptions& opt) :
    index_(index),
    tc_(tc),
    num_trans_(index.target_names_.size()),
    ecmap_(index.ecmap),
    aseCounts_(aseCounts),
    counts_(counts),
    target_names_(index.target_names_),
    post_bias_(4096,1.0),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    rho_(num_trans_, 0.0),
    rho_set_(false),
    all_fl_means(all_means),
    opt(opt),
    indexTransform(index_.target_names_.size(), -1)
  {
    assert(all_fl_means.size() == index_.target_lens_.size());
    eff_lens_ = calc_eff_lens(index_.target_lens_, all_fl_means);
    eff_lens_left_ = eff_lens_;
    eff_lens_right_ = eff_lens_;
    eff_nv_weights_ = eff_lens_;
    for (int i = 0; i < eff_lens_.size(); i++) {
      eff_lens_left_[i] -= (index_.target_lens_[i] - index_.target_left_lens_[i]);
      eff_lens_right_[i] -= (index_.target_lens_[i] - index_.target_right_lens_[i]);

      eff_nv_weights_[i] = (2 / (eff_lens_left_[i] + eff_lens_right_[i]));
    }

    for (int i = 0; i < index_.target_names_.size(); i++) {
      if (index.target_names_[i].find("_refsnp_") == std::string::npos) {
        indexTransform[i] = thetas.size();
        reverseIndexTransform.push_back(i);
        thetas.push_back(1);
        transcriptToEParam_.push_back({});
      }
    }

    
    w_e_fast_ = get_nonvariant_weight_map(ecmap_, eff_lens_left_, eff_lens_right_);
    std::cerr << "Finished making weights vector for nonvariant equivalence classes" << std::endl;

    std::vector<std::vector<int>> w_v_lambda_fast(index_.vcmap.size());
    std::vector<std::vector<double>> split_complex_multimap(index_.vcmap.size());
    std::unordered_set<int> interesting_snps;
    std::unordered_set<int> nonvariant_complex_transcripts;

    // want map from transcript id to the effective parameter for it. assuming that each transcript has one effective param for now.
    for (int vc = 0; vc < index_.vcmap.size(); vc++) {
      if (aseCounts_.at(vc) == 0) {
        std::cerr << "By construction, should never get here" << std::endl;
        continue;
      }

      auto &pair = index_.vcmap[vc];
      std::vector<int> trids = pair.first;
      vc_trids_.push_back({});
      /* vc_trids_unique_.push_back(trids); */
      vc_counts_.push_back(aseCounts_.at(vc));
      std::vector<int> effectiveParams = pair.second;
      
      w_v_lambda_fast[vc] = {};
      split_complex_multimap[vc] = {};
      std::unordered_set<int> remaining_ids(std::begin(trids), std::end(trids));
      for (int i = 0; i < trids.size(); i++) {
        for (int param : effectiveParams) { // if p negative, evidence for other side
          int ind_p = (param > 0) ? param : -1 * param;
          if (index_.effectiveASEParams.at(ind_p)[0].trid == trids[i]) { // the parameter in this vc belongs to this transcript
            w_v_lambda_fast[vc].push_back(param);
            split_complex_multimap[vc].push_back(1.0);
            vc_trids_[vc_trids_.size() - 1].push_back(trids[i]);
            remaining_ids.erase(trids[i]);

            auto& eparams = transcriptToEParam_[indexTransform[trids[i]]];
            if (std::find(eparams.begin(), eparams.end(), ind_p) == eparams.end()) {
              eparams.push_back(ind_p);
            }
          }
        }
      }

      for (int rtrid : remaining_ids) {
        w_v_lambda_fast[vc].push_back(0);
        nonvariant_complex_transcripts.insert(rtrid);
        split_complex_multimap[vc].push_back(2.0);
        /* split_complex_multimap[vc].push_back(1.0); */
        vc_trids_[vc_trids_.size() - 1].push_back(rtrid);
      }

      // basically want to check and see if the last entry of w_v_lambda contains a 0 but also 
      auto &lam = w_v_lambda_fast[vc];
      bool is_interesting = false;
      for (int param : lam) {
        if (param == 0) {
          is_interesting = true;
        }
      }

      if (is_interesting) {
        for (int param : lam) {
          if (param != 0) {
            interesting_snps.insert((param > 0) ? param : -1 * param);
          }
        }
      }

    }

    /* // complex transcripts are transcripts can be both variant and nonvariant */
		/* std::ofstream snp_info(opt.output + "/complex_snps.txt"); */
		/* std::ofstream complex_info(opt.output + "/complex_transcripts.txt"); */
    /* for (auto param : interesting_snps) { */
    /*   for (auto& vpe : index_.effectiveASEParams.at(param)) { */
    /*     int trid = vpe.trid; */
    /*     int pos = vpe.pos; */
    /*     snp_info << index_.target_names_[trid] << "_" << pos << "\t"; */
    /*     complex_info << index_.target_names_[trid] << "\n"; */
    /*   } */
    /*   snp_info << std::endl; */
    /* } */

    /* for (int trid : nonvariant_complex_transcripts) { */
    /*   complex_info << index_.target_names_[trid] << "\n"; */
    /* } */

    
    std::cerr << "Finished making weights vector and lambda vector for variant equivalence classes" << std::endl;
    w_v_lambda_fast_ = w_v_lambda_fast;
    split_complex_multimap_ = split_complex_multimap;
    w_v_fast_ = get_variant_map(vc_trids_, w_v_lambda_fast_, eff_lens_left_, eff_lens_right_);
    assert(target_names_.size() == eff_lens_.size());

    std::cerr << "Number of variant classes: " << index_.vcmap.size() << std::endl;

		/* std::ofstream tid_info(opt.output + "/tid2.dat"); */
    /* for (int i = 0; i < num_trans_; i++) { */
      /* if (index_.target_names_[i].find("_refsnp_") == std::string::npos) { */
        /* tid_info << i << "\t" << index_.target_names_[i] << std::endl; */
      /* } */
    /* } */

    /* std::ofstream ec_info(opt.output + "/ec2.dat"); */
    /* std::vector<int> emp; */
    /* for (int i = 0; i < ecmap_.size(); i++) { */
      /* /1* if (counts_[i] == 0) continue; *1/ */
      /* ec_info << counts_[i] << "\t" << ecmap_[i] << "\t" << emp << std::endl; */
      /* /1* std::cerr << counts_[i] << "\t" << ecmap_[i] << "\t" << emp << std::endl; *1/ */
    /* } */

    /* for (int i = 0; i < vc_trids_.size(); i++) { */
      /* /1* if (vc_counts_[i] == 0) continue; *1/ */
      /* ec_info << vc_counts_[i] << "\t" << vc_trids_[i] << "\t" << w_v_lambda_fast_[i] << std::endl; */
      /* /1* std::cerr << vc_counts_[i] << "\t" << vc_trids_[i] << "\t" << w_v_lambda_fast_[i] << std::endl; *1/ */
    /* } */
  }

  ~EMAlgorithm() {}

  std::vector<std::vector<double>> get_nonvariant_weight_map(
      const EcMap& ecmap,
      const std::vector<double>& eff_lens_left_,
      const std::vector<double>& eff_lens_right_) {
    WeightMap w_e_fast(ecmap.size());
    // first T ec's are just the transcripts themselves. so, only iterate through nontrivial equivalence classes
    for (int ec = num_trans_; ec < ecmap.size(); ec++) {
      std::vector<int> trids = ecmap[ec];
      for (int i = 0; i < trids.size(); i++) {
        // TODO: waste to store the first index_.target_names_.size() elements of this
        w_e_fast[ec].push_back(2 / (eff_lens_left_[trids[i]] + eff_lens_right_[trids[i]])); 
      }
    }
    return w_e_fast;
  }

  std::vector<std::vector<double>> get_variant_map(
      const std::vector<std::vector<int>> vc_trids,
      const std::vector<std::vector<int>> w_v_lambda_fast_,
      const std::vector<double>& eff_lens_left_,
      const std::vector<double>& eff_lens_right_) {
    std::vector<std::vector<double>> w_v_fast(vc_trids.size());
    for (int vc = 0; vc < vc_trids.size(); vc++) {
      for (int i = 0; i < vc_trids[vc].size(); i++) {
        if (w_v_lambda_fast_[vc][i] > 0) {
          w_v_fast[vc].push_back(1 / eff_lens_left_[vc_trids[vc][i]]); 
        } else if (w_v_lambda_fast_[vc][i] < 0) {
          w_v_fast[vc].push_back(1 / eff_lens_right_[vc_trids[vc][i]]);
        } else {
          /* Need to change this to be linear combination of the two */
         /* this is where we push the weight for this count. if param > 0, left. < 0, right. */
          w_v_fast[vc].push_back(2 / (eff_lens_left_[vc_trids[vc][i]] + eff_lens_right_[vc_trids[vc][i]])); 
        }
      }
    }
    return w_v_fast;
  }

  void printVC(int vc) {
    auto& pair = index_.vcmap[vc];
    std::vector<int> trids = pair.first;
    std::vector<int> effectiveParams = pair.second;
    std::cerr << "VC " << vc << ": " << aseCounts_.at(vc) << std::endl;
    for (int t : trids) {
        std::cerr << index_.target_names_[t] << std::endl;
        for (int param : effectiveParams) {
          int ind_p = (param > 0) ? param : -1 * param;
          auto& vpe = index_.effectiveASEParams.at(ind_p)[0];
          if (vpe.trid == t) {
            if (param > 0) {
              std::cerr << "\t" << index_.target_names_[t] << " (Eff Len " << eff_lens_[t] << ")" << ": L, " << vpe.pos << std::endl;
            } else {
              std::cerr << "\t" << index_.target_names_[t] << " (Eff len " << eff_lens_[t] << ")" << ": R, " << vpe.pos << std::endl;
            }
          }

          std::cerr << std::endl;
        }
    }
  }

  double convergence(int conv) {
    if (conv == 1) {
      return 5 * 1e-2;
    } else if (conv == 2) {
      return 1e-2;
    } else if (conv == 3) {
      return 1e-3;
    } else if (conv == 4) {
      return 1e-4;
    } else if (conv == 5) {
      return 1e-5;
    } else {
      return 1e-2;
    }
  }

  void run(size_t n_iter = 10000, size_t min_rounds=50, bool verbose = true, bool recomputeEffLen = true) {
    std::cerr << "Number of variant equivalence classes: " << index_.vcmap.size() << std::endl;
    const double alpha_limit = 1e-7;
    const double alpha_change_limit = 1e-2;
    const double alpha_change = convergence(opt.conv);
    std::cerr << alpha_change << ", " << opt.smooth << std::endl;
    bool finalRound = false;

    // TESTING FOR DIFFERENT VALUES
    bool phased_em = opt.phased;
    assert(thetas.size() == target_names_.size() - index_.shadeToColorTranscriptMap.size());
    // relevant allele specific parameters
    std::vector<double> ps(index_.effectiveASEParams.size(), .5);
    for (int j = 0; j < ps.size(); j++) {
      ps[j] = .5; //((double) std::rand() / (RAND_MAX));
      if (phased_em) {
        if (j == 0) continue;
        auto& vpe = index_.effectiveASEParams.at(j);
        double left_len = eff_lens_left_[vpe[0].trid];
        double right_len = eff_lens_right_[vpe[0].trid];
        ps[j] = (left_len / (left_len + right_len));
      }
      pdenoms.push_back(.5);
    }

    // calculate weights for each
    int num_iterations;
    std::vector<double> theta_sum(thetas.size(), 0.0);
    std::vector<double> newPs(ps.size(), 0.0);
    std::vector<double> denoms(ps.size(), 0.0);
    double eta_t = 1;
    double eta_p = 1e-2;

    time_t start, end;
    time(&start);

    /* Implement kallisto EM update equations to see if it's a pseudoalignment issue or if it's something else */
    std::vector<double> diploid_thetas(2 * thetas.size(), 1.0);
    std::vector<double> diploid_thetas_sum(2 * thetas.size(), 0.0);
    for (num_iterations = 0; num_iterations < 1000000; num_iterations++) {
      if (recomputeEffLen && (num_iterations == min_rounds || num_iterations == min_rounds + 500)) {
        std::vector<double> alphas(index_.target_names_.size());
        for (int i = 0; i < num_trans_; i++) {
          if (index_.target_names_[i].find("_refsnp_") == std::string::npos) {
            alphas[i] = thetas[indexTransform[i]];
          }
        }

        eff_lens_ = update_eff_lens(all_fl_means, tc_, index_, alphas, eff_lens_, post_bias_, opt);
        /* Update the weight maps here with the new eff lens */
        w_e_fast_ = get_nonvariant_weight_map(ecmap_, eff_lens_left_, eff_lens_right_);
        w_v_fast_ = get_variant_map(vc_trids_, w_v_lambda_fast_, eff_lens_left_, eff_lens_right_);
        /* if (!phased_em) { */
        /*   w_v_unique_ = get_variant_map(vc_trids_unique_, eff_lens_left_, eff_lens_right_); */
        /* } */
      }

      for (int i = 0; i < thetas.size(); i++) {
        int t = reverseIndexTransform[i]; // t is kallisto index for transcript
        theta_sum[i] = counts_[t];
      }

      for (int ec = index_.target_names_.size(); ec < ecmap_.size(); ec++) {
        if (counts_[ec] == 0) {
          continue;
        }

        double normalize = 0;
        auto& trids = ecmap_[ec];
        std::vector<double>& weights_for_ec = w_e_fast_[ec];
        for (int i = 0; i < trids.size(); i++) {
          double w = 0.0;
          if (!phased_em) {
            w = 1 / (eff_lens_left_[trids[i]] + eff_lens_right_[trids[i]]);
          } else {
            w = eff_nv_weights_[trids[i]];
          }
          weights_for_ec[i] = w;
          normalize += thetas[indexTransform[trids[i]]] * weights_for_ec[i];
        }

        for (int i = 0; i < trids.size(); i++) {
          double contribution = (counts_[ec] * thetas[indexTransform[trids[i]]] * weights_for_ec[i] / normalize);
          theta_sum[indexTransform[trids[i]]] += contribution;
        }
      }

      /* VARIANT EQUIVALENCE CLASSES */
      if (phased_em) {
        for (int vc = 0; vc < vc_trids_.size(); vc++) {
          if (vc_counts_[vc] == 0) {
            continue;
          }

          double normalize = 0;
          auto& trids = vc_trids_[vc];
          std::vector<int>& lambda_vec = w_v_lambda_fast_[vc];
          std::vector<double>& weights_for_vc = w_v_fast_[vc];
          for (int i = 0; i < trids.size(); i++) {
            if (lambda_vec[i] == 0) {
              normalize += thetas[indexTransform[trids[i]]] * eff_nv_weights_[trids[i]];
            } else if (lambda_vec[i] > 0) {
              normalize += thetas[indexTransform[trids[i]]] * ps[lambda_vec[i]] * weights_for_vc[i];
            } else {
              normalize += thetas[indexTransform[trids[i]]] * (1 - ps[-1 * lambda_vec[i]]) * weights_for_vc[i];
            }
          }

          for (int i = 0; i < trids.size(); i++) {
            double contribution;
            if (lambda_vec[i] == 0) {
              contribution = 1 * thetas[indexTransform[trids[i]]] * eff_nv_weights_[trids[i]] / normalize;
            } else if (lambda_vec[i] > 0) {
              contribution = ps[lambda_vec[i]] * thetas[indexTransform[trids[i]]] * weights_for_vc[i] / normalize;
              newPs[lambda_vec[i]] += vc_counts_[vc] * contribution;
              denoms[lambda_vec[i]] += vc_counts_[vc] * contribution;
            } else {
              contribution = (1 - ps[-1 * lambda_vec[i]]) *  thetas[indexTransform[trids[i]]] * weights_for_vc[i] / normalize;
              denoms[-1 * lambda_vec[i]] += vc_counts_[vc] * contribution; 
            }
            theta_sum[indexTransform[trids[i]]] += vc_counts_[vc] * contribution;
          }
        }
      } else {
        for (int vc = 0; vc < vc_trids_.size(); vc++) {
          if (vc_counts_[vc] == 0) {
            continue;
          }

          double normalize = 0;
          auto& trids = vc_trids_[vc];
          auto& trid_weights = split_complex_multimap_[vc];
          for (int i = 0; i < trids.size(); i++) {
            double w = (trid_weights[i] / (eff_lens_left_[trids[i]] + eff_lens_right_[trids[i]]));
            normalize += thetas[indexTransform[trids[i]]] * w;
          }

          for (int i = 0; i < trids.size(); i++) {
            double w = (trid_weights[i] / (eff_lens_left_[trids[i]] + eff_lens_right_[trids[i]]));
            double contribution = thetas[indexTransform[trids[i]]] * w / normalize;
            theta_sum[indexTransform[trids[i]]] += vc_counts_[vc] * contribution;
          }
        }
      }

      // update parameters
      bool stopEM = false;
      int chcount = 0;
      /* double renorm = 0.0; */
      for (int i = 0; i < thetas.size(); i++) {
        int true_trid = reverseIndexTransform[i];
        if (phased_em) {
          for (int ind_p : transcriptToEParam_[i]) {
            if (denoms[ind_p] > 0) {
              double bias = (finalRound) ? 0 : 1e-2;
              bias = (opt.smooth ? bias : 0);
              ps[ind_p] = (bias + newPs[ind_p]) / (2*bias + denoms[ind_p]);
              if (eff_lens_left_[true_trid] != eff_lens_right_[true_trid]) {
                eff_nv_weights_[true_trid] = 
                  (ps[ind_p] / eff_lens_left_[true_trid]) + ((1 - ps[ind_p]) / eff_lens_right_[true_trid]);
              }
              pdenoms[ind_p] = denoms[ind_p];
              newPs[ind_p] = 0.0;
              denoms[ind_p] = 0.0;
            } else {
              ps[ind_p] = .5;
              pdenoms[ind_p] = 0;
            }
          }
        }

        double newVal = theta_sum[i];
        if (newVal > alpha_change_limit && (std::fabs(newVal - thetas[i]) / (thetas[i]) > alpha_change)) {
          chcount++;
        }
        thetas[i] = newVal;
        theta_sum[i] = 0.0;
      }

      if (chcount == 0 && num_iterations > min_rounds) {
        stopEM = true;
      }


      if (finalRound) {
        break;
      }

      if (stopEM) {
        finalRound = true;
        for (int i = 0; i < thetas.size(); i++) {
          if (thetas[i] < 2*alpha_limit/10.0) {
            thetas[i] = 0.0;
          }
        }
      }
    }

    time(&end);
    std::cerr << "EM actually took " << difftime(end, start) << " seconds." << std::endl;
    std::cerr << "EM Algorithm took " << num_iterations << " rounds" << std::endl;


    double total {0.0};
    for (auto i = 0; i < thetas.size(); i++) {
      double t = thetas[i] * eff_nv_weights_[reverseIndexTransform[i]];
      tpms.push_back(t);
      total += t;
    }

    for (auto& t : tpms) {
      t /= total;
      t *= 1000000;
    }

    std::ofstream out_data(opt.output + "/thetas.txt");
    std::ofstream out_tpms(opt.output + "/tpms.txt");
    out_data << "{" << std::endl;
    out_tpms << "{" << std::endl;
    for (int i = 0; i < index_.target_names_.size(); i++) {
      if (index_.target_names_[i].find("_refsnp_") == std::string::npos) {
        out_data << "  " << index_.target_names_[i] << ": " << ((thetas[indexTransform[i]])) << "," << std::endl;
        out_tpms << "  " << index_.target_names_[i] << ": " << ((tpms[indexTransform[i]])) << "," << std::endl;
      }
    }
    out_data << "}" << std::endl;
    out_tpms << "}" << std::endl;

    if (phased_em) {
      std::ofstream out_data_p(opt.output + "/ps.txt");
      std::ofstream out_data_pdenoms(opt.output + "/pdenoms.txt");
      out_data_p.precision(28);
      out_data_p << "{" << std::endl;
      // 0 doesn't exist as a key
      for (auto& kv : index_.effectiveASEParams) {
        if (kv.second.size() == 0) {
          continue;
        }
        auto& vpe = kv.second[0];

        if (vpe.trid < index_.target_names_.size() && ps[kv.first] != .5) {
          out_data_p << "  " 
                     << index_.target_names_[vpe.trid] << "_L_" 
                     << vpe.pos << ": " << ps[kv.first] << "," << std::endl;
        }

        if (vpe.trid < index_.target_names_.size()) {
          out_data_pdenoms << "  " 
                           << index_.target_names_[vpe.trid] << "_L_" 
                           << vpe.pos << ": " << pdenoms[kv.first] << "," << std::endl;
        }
      }
      out_data_p << "}" << std::endl;
    } else {
      std::vector<double> p_refs(index_.effectiveASEParams.size(), 0.0);
      std::vector<double> p_alts(index_.effectiveASEParams.size(), 0.0);
      /* Extract expected read counts over loci... */
      for (int vc = 0; vc < vc_trids_.size(); vc++) {
        if (vc_counts_[vc] == 0) {
          continue;
        }
        
        /* Define a map of nonvariant posterior weight for each transcript in this equivalence class */
        std::unordered_map<int, double> trid_to_nv_posterior;
        double normalize = 0;
        auto& trids_unique = vc_trids_[vc];
        auto& trids_weights = split_complex_multimap_[vc];
        /* auto& weights_for_vc_unique = w_v_unique_[vc]; */
        for (int i = 0; i < trids_unique.size(); i++) {
          double w = (trids_weights[i] / (eff_lens_left_[trids_unique[i]] + eff_lens_right_[trids_unique[i]]));
          normalize += thetas[indexTransform[trids_unique[i]]] * w;
        }

        for (int i = 0; i < trids_unique.size(); i++) {
          double w = (trids_weights[i] / (eff_lens_left_[trids_unique[i]] + eff_lens_right_[trids_unique[i]]));
          trid_to_nv_posterior[trids_unique[i]] =
            thetas[indexTransform[trids_unique[i]]] * w / normalize;
        }

        auto& trids = vc_trids_[vc];
        std::vector<int>& lambda_vec = w_v_lambda_fast_[vc];

        for (int i = 0; i < trids.size(); i++) {
          /* In a given equivalence class, we only witness either ref or alt for a parameter */
          double v_weight = 1.0;
          if (tc_.effParamsBias.find(lambda_vec[i]) != tc_.effParamsBias.end()) {
            v_weight = tc_.effParamsBias.at(lambda_vec[i]);
          }

          if (lambda_vec[i] > 0) {
            p_refs[lambda_vec[i]] += vc_counts_[vc] * trid_to_nv_posterior[trids[i]] * v_weight;
          } else {
            p_alts[-1 * lambda_vec[i]] += vc_counts_[vc] * trid_to_nv_posterior[trids[i]] * v_weight;
          }
        }
      }

      std::ofstream out_data_read_counts(opt.output + "/allele_counts.txt");

      /* Output format should be chr, pos, ref, alt, alt2, gt (unphased, so NA), ref, alt, other */
      // 0 doesn't exist as a key
      for (auto& kv : index_.effectiveASEParams) {
        if (kv.second.size() == 0) {
          continue;
        }
        auto& vpe = kv.second[0];
        if (vpe.trid < index_.target_names_.size()) {
          out_data_read_counts << index_.target_names_[vpe.trid] << "\t"
                               << vpe.pos + 1 << "\t"
                               << vpe.ref << "\t"
                               << vpe.alt << "\t"
                               << "NA" << "\t"
                               << p_refs[kv.first] << "\t"
                               << p_alts[kv.first] << "\t"
                               << 0 << "\n";
        }
      }
    }
  }

  void compute_rho() {
    if (rho_set_) {
      // rho has already been set, let's clear it
      std::fill(rho_.begin(), rho_.end(), 0.0);
    }

    double total {0.0};
    for (auto i = 0; i < alpha_.size(); ++i) {
      if (eff_lens_[i] < TOLERANCE) {
        std::cerr << "Should actually never really get here... tid: "  << i <<
            std::endl;
        continue;
      }
      rho_[i] = thetas[indexTransform[i]]; // this is fucked
      total += rho_[i];
    }

    for (auto& r : rho_) {
      r /= total;
    }

    rho_set_ = true;
  }

  // DEPRECATED:
  void write(const std::string& out_fname) const {
    std::ofstream out;
    out.open(out_fname, std::ios::out);

    if (!out.is_open()) {
      std::cerr << "Error opening '" << out_fname << "'" <<
          std::endl;
      exit(1);
    }

    out.precision(15);

    out <<
        "target_id" << "\t" <<
        "kallisto_id" << "\t" <<
        "rho" << "\t" <<
        "tpm" << "\t" <<
        "est_counts" <<
        std::endl;

    const double MILLION = 1e6;

    for (auto i = 0; i < rho_.size(); ++i) {
      out <<
          target_names_[i] << "\t" <<
          i << "\t" <<
          rho_[i] << "\t" <<
          rho_[i] * MILLION << "\t" <<
          alpha_[i] <<
          std::endl;
    }

    out.flush();
    out.close();
  }

  void set_start(const EMAlgorithm& em_start) {
    assert(em_start.alpha_before_zeroes_.size() == alpha_.size());
    double big = 1.0;
    double sum_counts = std::accumulate(counts_.begin(), counts_.end(), 0.0);
    double sum_big = 0.0;
    int count_big = 0;
    for (auto x : em_start.alpha_before_zeroes_) {
      if (x >= big) {
        sum_big += x;
        count_big++;
      }
    }
    int n = alpha_.size();
    for (auto i = 0; i < n; i++) {
      if (em_start.alpha_before_zeroes_[i] >= big) {
        alpha_[i] = em_start.alpha_before_zeroes_[i];
      } else {
        alpha_[i] = sum_counts/(n - count_big);
      }
    }

    //std::cout << sum_big << " " << count_big << " " << n << std::endl;

    std::copy(em_start.alpha_before_zeroes_.begin(), em_start.alpha_before_zeroes_.end(),
        alpha_.begin());
  }


  int num_trans_;
  const KmerIndex& index_;
  const MinCollector& tc_;
  const EcMap& ecmap_;
  const std::vector<int>& counts_;
  const std::unordered_map<int, int>& aseCounts_;
  const std::vector<std::string>& target_names_;
  const std::vector<double>& all_fl_means;
  std::vector<int> vc_counts_;
  std::vector<double> eff_lens_;
  std::vector<double> eff_nv_weights_;
  std::vector<double> eff_lens_left_;
  std::vector<double> eff_lens_right_;
  std::vector<double> post_bias_;
  WeightMap weight_map_;
  std::vector<double> alpha_;
  std::vector<double> alpha_before_zeroes_;
  std::vector<double> thetas;
  std::vector<double> tpms;
  std::vector<double> pdenoms;
  std::vector<std::vector<double>> w_e_fast_;
  std::vector<std::vector<double>> w_v_fast_;
  std::vector<std::vector<int>> w_v_lambda_fast_;
  std::vector<std::vector<double>> split_complex_multimap_;
  std::vector<std::vector<int>> vc_trids_;
  std::vector<std::vector<int>> transcriptToEParam_;
  std::vector<int> indexTransform; // transcript index to theta index
  std::vector<int> reverseIndexTransform; // theta index to transcript index
  std::vector<double> rho_;
  bool rho_set_;
  const ProgramOptions& opt;
};


#endif // KALLISTO_EMALGORITHM_H
