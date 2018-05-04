#ifndef READ_PROBABILITY_CALCULATOR_H__
#define READ_PROBABILITY_CALCULATOR_H__

#include <map>
#include "path_aligner.h"

#include "config.pb.h"

struct SingleProbabilityChange {
  vector<Path> added_paths;
  vector<Path> removed_paths;

  vector<SingleReadAlignment> added_alignments;
  vector<SingleReadAlignment> removed_alignments;

  int new_paths_length;

  vector<Path> new_paths;
};

struct PairedProbabilityChange {
  vector<Path> added_paths;
  vector<Path> removed_paths;

  vector<PairedReadAlignment> added_alignments;
  vector<PairedReadAlignment> removed_alignments;

  int new_paths_length;

  vector<Path> new_paths;
};

struct HICProbabilityChange {
  vector<Path> added_paths;
  vector<Path> removed_paths;
  vector<Path> still_paths;

  int new_paths_length;

  vector<Path> new_paths;
};



struct ProbabilityChanges {
  // single object for every probability change
  vector<SingleProbabilityChange> single_read_changes;
  vector<PairedProbabilityChange> paired_read_changes;
  vector<HICProbabilityChange> hic_read_changes;

  long long getLength() {
    if (!single_read_changes.empty()) return single_read_changes[0].new_paths_length;
    if (!paired_read_changes.empty()) return paired_read_changes[0].new_paths_length;
    if (!hic_read_changes.empty()) return hic_read_changes[0].new_paths_length;
    return 0;
  }
};

class SingleReadProbabilityCalculator {
 public:
  SingleReadProbabilityCalculator(
      SingleShortReadSet<>* read_set, //double mismatch_prob,
      double insert_prob, double del_prob, double subst_prob,
      double min_prob_start, double min_prob_per_base,
      double penalty_constant, int penalty_step) :
        read_set_(read_set), path_aligner_(read_set),
        //mismatch_prob_(mismatch_prob),
        insert_prob_(insert_prob), del_prob_(del_prob), subst_prob_(subst_prob),
        min_prob_start_(min_prob_start), min_prob_per_base_(min_prob_per_base),
        penalty_constant_(penalty_constant), penalty_step_(penalty_step),
        old_paths_length_(1) {
    read_probs_.resize(read_set_->size());
    total_log_prob_ = InitTotalLogProb();
    total_reads_length_ = 0;
    for (size_t i = 0; i < read_set_->size(); i++) total_reads_length_ += (*read_set_)[i].size();
  }

  // Call this first
  double GetPathsProbability(
      const vector<Path>& paths, SingleProbabilityChange& prob_change);

  // Call this after you are happy with current result (i.e. you got better
  // probability)
  void CommitProbabilityChange(const SingleProbabilityChange &prob_change);

  int GetPathsLength(const vector<Path>& paths) const;
  // Evals change with filled added and removed paths
  void EvalProbabilityChange(SingleProbabilityChange& prob_change);

  void RemovePathFromCache(const Path& p) {
    path_aligner_.RemovePathFromCache(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    path_aligner_.RemovePathsFromCache(paths);
  }

  long long GetTotalReadsLength() {
    return total_reads_length_;
  }

  int GetUnalignedReadsCount();

  int GetReadsCount() {
    return (int)read_set_->size();
  }

  string GetProbHist() {
    stringstream ss;
    map<int, int> log_bins;
    for (int i = 0; i < read_probs_.size(); i++) {
      const double prob = read_probs_[i];
      const int bin = (int)ceil(log(prob));
      log_bins[bin] = 1 + log_bins[bin];
    }

    const int threshold_bin = (int)ceil(log(GetMinLogProbability((*read_set_)[0].size())));
    for (int i = 0; i >= -1000; i--) {
      if (i > threshold_bin) {
        ss << log_bins[i] << "\t";
      }
      else {
        ss << "**" << log_bins[i] << "\t";
      }
      if (abs(i) % 10 == 0 ) ss << "||";
    }
    return ss.str();
  }

 private:
  double InitTotalLogProb();

  double GetMinLogProbability(int read_length) const;

  // max(min_prob, prob)
  double GetRealReadProbability(double prob, int read_id) const;



  // Get total probability from change and cached data
  double EvalTotalProbabilityFromChange(const SingleProbabilityChange& prob_change, bool write=false);



  double GetAlignmentProbSimple(int dist, int read_length) const;
  double GetAlignmentProb(int matches, int inserts, int dels, int substs) const;

  SingleShortReadSet<>* read_set_;
  SingleShortReadPathAligner path_aligner_;
  //double mismatch_prob_;
  double insert_prob_;
  double del_prob_;
  double subst_prob_;
  double min_prob_start_;
  double min_prob_per_base_;
  double penalty_constant_;
  int penalty_step_;

  long long total_reads_length_;

  vector<double> read_probs_;
  double total_log_prob_;
  int old_paths_length_;
  vector<Path> old_paths_;
};

class PairedReadProbabilityCalculator {
 public:
  PairedReadProbabilityCalculator(
      ShortPairedReadSet<>* read_set,
      //double mismatch_prob,
      double insert_prob,
      double del_prob,
      double subst_prob,
      double min_prob_start,
      double min_prob_per_base,
      double penalty_constant,
      int penalty_step,
      double mean_distance,
      double std_distance,
      bool use_as_advice,
      int uncovered_threshold,
      double uncovered_penalty,
      int uncovered_start_ignore,
      int lower_bound_advice_count
  ): read_set_(read_set), path_aligner_(read_set), //mismatch_prob_(mismatch_prob),
     insert_prob_(insert_prob), del_prob_(del_prob), subst_prob_(subst_prob),
     min_prob_start_(min_prob_start), min_prob_per_base_(min_prob_per_base),
     penalty_constant_(penalty_constant), penalty_step_(penalty_step),
     old_paths_length_(1), mean_distance_(mean_distance),
     std_distance_(std_distance), use_as_advice_(use_as_advice),
     uncovered_threshold_(uncovered_threshold), uncovered_penalty_(uncovered_penalty),
     uncovered_start_ignore_(uncovered_start_ignore),
     lower_bound_advice_count_(lower_bound_advice_count){
    read_probs_.resize(read_set_->size());
    total_log_prob_ = InitTotalLogProb();

    read_als_.resize(read_set_->size(), 0);
    read_part_als_1.resize(read_set_->size(), 0);
    read_part_als_2.resize(read_set_->size(), 0);

    total_reads_length_ = 0;
    for (size_t i = 0; i < read_set_->size(); i++) total_reads_length_ += read_set_->reads_1_[i].size() + read_set_->reads_2_[i].size();
    uncovered_bases_count_ = 0;

  }
  // Call this first
  double GetPathsProbability(
      const vector<Path>& paths, PairedProbabilityChange& prob_change);

  // Call this after you are happy with current result (i.e. you got better
  // probability)
  void CommitProbabilityChange(const PairedProbabilityChange& prob_change);


  bool use_as_advice_;
  PairedReadPathAligner path_aligner_;
  double mean_distance_;
  double std_distance_;

  int uncovered_threshold_;
  double uncovered_penalty_;
  int uncovered_bases_count_;
  int uncovered_start_ignore_;
  int lower_bound_advice_count_;

  int GetPathsLength(const vector<Path>& paths) const;
  void RemovePathFromCache(const Path& p) {
    path_aligner_.RemovePathFromCache(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    path_aligner_.RemovePathsFromCache(paths);
  }

  long long GetTotalReadsLength() {
    return total_reads_length_;
  }
  int GetUnalignedReadsCount();

  double GetAvgReadProb() {
    double sum = 0;
    int count = 0;
    for(int i = 0; i < read_probs_.size(); i++) {
      if (read_als_[i] > 0) {
        sum += read_probs_[i];
        count += 1;
      }
    }
    return sum / (double)count;
  }

  string GetProbHist() {
    stringstream ss;
    map<int, int> log_bins;
    for (int i = 0; i < read_probs_.size(); i++) {
      if (read_als_[i] == 0) continue;
      const double prob = read_probs_[i];
      if (prob <= 0) log_bins[1] = 1 + log_bins[1];
      else {
        const int bin = (int)ceil(log(prob));
        log_bins[bin] = 1 + log_bins[bin];
      }
    }

    const int threshold_bin = (int)ceil(GetMinLogProbability((*read_set_)[0].first.size() +(*read_set_)[0].second.size()));

    for (int i = 1; i >= -1000; i--) {
      if (i == 1) {
        ss << "[" << log_bins[1] << "]\t";
      }
      else if (i > threshold_bin) {
        ss << log_bins[i] << "\t";
      }
      else {
        ss << "**" << log_bins[i] << "\t";
      }
      if (abs(i) % 10 == 0 ) ss << "||";
    }
    return ss.str();
  }

  int GetReadsCount() {
    return (int)read_set_->size();
  }

  int GetCompletelyUnalignedReadsCount();

  int GetUncoveredBasesCount(const Path& p);

  double GetAlignmentProb(const PairedReadAlignment& al) const;
 private:

  // Evals change with filled added and removed paths
  void EvalProbabilityChange(PairedProbabilityChange& prob_change, bool debug_output=true);
  // Get total probability from change and cached data
  double EvalTotalProbabilityFromChange(const PairedProbabilityChange& prob_change, bool write=false) ;

  double InitTotalLogProb();

  double GetMinLogProbability(int read_length) const;

  // max(min_prob, prob)
  double GetRealReadProbability(double prob, int read_id) const;

  ShortPairedReadSet<>* read_set_;
  //double mismatch_prob_;
  double insert_prob_;
  double del_prob_;
  double subst_prob_;
  double min_prob_start_;
  double min_prob_per_base_;
  double penalty_constant_;
  int penalty_step_;
  int old_paths_length_;
  double total_log_prob_;
  //double mean_distance_;
  //double std_distance_;

  long long total_reads_length_;

  vector<double> read_probs_;
  vector<Path> old_paths_;

  vector<int> read_als_;
  vector<int> read_part_als_1, read_part_als_2;

};

class HICReadProbabilityCalculator {
 public:
  HICReadProbabilityCalculator(
      HICReadSet<>* read_set,
      double insert_prob,
      double del_prob,
      double subst_prob,
      double translocation_prob,
      int binsize,
      double min_prob_start,
      double min_prob_per_base,
      double penalty_constant,
      int penalty_step,
      bool use_as_advice,
      int lower_bound_advice_count
  ): read_set_(read_set), path_aligner_(HICReadPathAligner(read_set, binsize)), //mismatch_prob_(mismatch_prob),
     insert_prob_(insert_prob), del_prob_(del_prob), subst_prob_(subst_prob),
     translocation_prob_(translocation_prob),
     min_prob_start_(min_prob_start), min_prob_per_base_(min_prob_per_base),
     penalty_constant_(penalty_constant), penalty_step_(penalty_step),
     old_paths_length_(1), binsize_(binsize), use_as_advice_(use_as_advice),
     lower_bound_advice_count_(lower_bound_advice_count){
    //read_probs_.resize(read_set_->size());
    read_cis_phis_.resize(read_set_->size());
    read_trans_psis_.resize(read_set_->size());
    total_log_prob_ = InitTotalLogProb();

    total_reads_length_ = 0;
    for (size_t i = 0; i < read_set_->size(); i++) total_reads_length_ += read_set_->reads_1_[i].size() + read_set_->reads_2_[i].size();

    read_cis_als_.resize(read_set_->size(), 0);
    read_trans_als_.resize(read_set_->size(), 0);

    cis_constant_ = 1;
    trans_constant_ = 1;
  }
  // Call this first
  double GetPathsProbability(
      const vector<Path>& paths, HICProbabilityChange& prob_change);

  // Call this after you are happy with current result (i.e. you got better
  // probability)
  void CommitProbabilityChange(const HICProbabilityChange& prob_change);


  bool use_as_advice_;
  HICReadPathAligner path_aligner_;

  // Evals change with filled added and removed paths
  void EvalProbabilityChange(HICProbabilityChange& prob_change, bool debug_output=true);
  // Get total probability from change and cached data
  double EvalTotalProbabilityFromChange(const HICProbabilityChange &prob_change,
                                        bool write=false);
  int GetPathsLength(const vector<Path>& paths) const;
  void RemovePathFromCache(const Path& p) {
    path_aligner_.RemovePathFromCache(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    path_aligner_.RemovePathsFromCache(paths);
  }

  long long GetTotalReadsLength() {
    return total_reads_length_;
  }

  int GetReadsCount() {
    return (int)read_set_->size();
  }

  int GetUnalignedReadsCount();
  int GetCompletelyUnalignedReadsCount();

  double GetAvgReadProb() {
    double sum = 0;
    int count = 0;
    for(int i = 0; i < read_cis_als_.size(); i++) {
      if (read_cis_als_[i] + read_trans_als_[i] > 0) {
        sum += cis_constant_ * read_cis_phis_[i] + trans_constant_ * read_trans_als_[i];
        count += 1;
      }
    }
    return sum / (double)count;
  }

  string GetProbHist() {
    stringstream ss;
    map<int, int> log_bins;
    for (int i = 0; i < read_cis_als_.size(); i++) {
      if (read_cis_als_[i] + read_trans_als_[i] == 0) continue;
      const double prob = cis_constant_ * read_cis_phis_[i] + trans_constant_ * read_trans_als_[i];
      if (prob <= 0) log_bins[1] = 1 + log_bins[1];
      else {
        const int bin = (int)ceil(log(prob));
        log_bins[bin] = 1 + log_bins[bin];
      }
    }

    const int threshold_bin = (int)ceil(GetMinLogProbability((*read_set_)[0].first.size() +(*read_set_)[0].second.size()));
    for (int i = 1; i >= -1000; i--) {
      if (i == 1) {
        ss << "[" << log_bins[1] << "]\t";
      }
      else if (i > threshold_bin) {
        ss << log_bins[i] << "\t";
      }
      else {
        ss << "**" << log_bins[i] << "\t";
      }
      if (abs(i) % 10 == 0 ) ss << "||";
    }
    return ss.str();
  }

  int lower_bound_advice_count_;

 private:

  double InitTotalLogProb();

  double GetMinLogProbability(int read_length) const;

  double GetCisScore(SingleReadAlignment al1, int a1_len, SingleReadAlignment al2, int a2_len, double lambda);

  double GetTransScore (SingleReadAlignment al1, SingleReadAlignment al2);

  // max(min_prob, prob)
  double GetRealReadProbability(double prob, int read_id) const;

  //double GetAlignmentProb(const HICReadAlignment& al) const;

  HICReadSet<>* read_set_;
  //double mismatch_prob_;
  double insert_prob_;
  double del_prob_;
  double subst_prob_;
  double translocation_prob_;
  double min_prob_start_;
  double min_prob_per_base_;
  double penalty_constant_;
  int penalty_step_;
  int old_paths_length_;
  double total_log_prob_;
  int binsize_;

  double cis_constant_;
  double trans_constant_;

  long long total_reads_length_;

  //vector<double> read_probs_;
  vector<double> read_cis_phis_;
  vector<double> read_trans_psis_;
  vector<Path> old_paths_;

  vector<int> read_cis_als_;
  vector<int> read_trans_als_;

};

class GlobalProbabilityCalculator {
 public:
  GlobalProbabilityCalculator(const Config &config);

  // Call this first
  double GetPathsProbability(
      const vector<Path>& paths, ProbabilityChanges& prob_changes);

  // Call this after you are happy with current result (i.e. you got better
  // probability)
  void CommitProbabilityChanges(const ProbabilityChanges &prob_changes);

  double GetAprioriPathsLogProbability(const vector<Path>& paths);

  void RemovePathFromCache(const Path& p) {
    for (auto &calc: single_read_calculators_) {
      calc.first.RemovePathFromCache(p);
    }
    for (auto &calc: paired_read_calculators_) {
      calc.first.RemovePathFromCache(p);
    }
    for (auto &calc: hic_read_calculators_) {
      calc.first.RemovePathFromCache(p);
    }
  }

  void RemovePathsFromCache(const vector<Path>& paths) {
    for (auto &calc: single_read_calculators_) {
      calc.first.RemovePathsFromCache(paths);
    }
    for (auto &calc: paired_read_calculators_) {
      calc.first.RemovePathsFromCache(paths);
    }
    for (auto &calc: hic_read_calculators_) {
      calc.first.RemovePathsFromCache(paths);
    }
  }

  string GetUnalignedReadsDebug();
  string GetUnalignedReadsLog();

  string GetProbHists() {
    stringstream ss;
    for (auto &calc: single_read_calculators_) {
      ss << "S: " << calc.first.GetProbHist() << endl;
    }
    for (auto &calc: paired_read_calculators_) {
      ss << "P: " << calc.first.GetProbHist() << endl;
    }
    for (auto &calc: hic_read_calculators_) {
      ss << "H: " << calc.first.GetProbHist() << endl;
    }
    return ss.str();
  }

  int GetUncoveredBasesCount() {
    int res = 0;
    for (auto &pc: paired_read_calculators_) {
      res += pc.first.uncovered_bases_count_;
    }
    return res;
  }

  vector<SingleShortReadSet<>*> single_short_read_sets_;
  vector<ShortPairedReadSet<>*> paired_read_sets_;
  vector<HICReadSet<>*> hic_read_sets_;
  // (prob calculator, weight)
  vector<pair<SingleReadProbabilityCalculator, double>> single_read_calculators_;
  vector<pair<PairedReadProbabilityCalculator, double>> paired_read_calculators_;
  vector<pair<HICReadProbabilityCalculator, double>> hic_read_calculators_;

 private:


};

#endif
