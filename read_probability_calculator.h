#ifndef READ_PROBABILITY_CALCULATOR_H__
#define READ_PROBABILITY_CALCULATOR_H__

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
  // @TODO create an interface for ProbabilityChange?
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
      bool use_as_advice
  ): read_set_(read_set), path_aligner_(read_set), //mismatch_prob_(mismatch_prob),
     insert_prob_(insert_prob), del_prob_(del_prob), subst_prob_(subst_prob),
     min_prob_start_(min_prob_start), min_prob_per_base_(min_prob_per_base),
     penalty_constant_(penalty_constant), penalty_step_(penalty_step),
     old_paths_length_(1), mean_distance_(mean_distance),
     std_distance_(std_distance), use_as_advice_(use_as_advice) {
    read_probs_.resize(read_set_->size());
    total_log_prob_ = InitTotalLogProb();

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

  int GetPathsLength(const vector<Path>& paths) const;
  void RemovePathFromCache(const Path& p) {
    path_aligner_.RemovePathFromCache(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    path_aligner_.RemovePathsFromCache(paths);
  }
 private:

  // Evals change with filled added and removed paths
  void EvalProbabilityChange(PairedProbabilityChange& prob_change, bool debug_output=true);
  // Get total probability from change and cached data
  double EvalTotalProbabilityFromChange(const PairedProbabilityChange& prob_change, bool write=false) ;

  double InitTotalLogProb();

  double GetMinLogProbability(int read_length) const;

  // max(min_prob, prob)
  double GetRealReadProbability(double prob, int read_id) const;

  double GetAlignmentProb(const PairedReadAlignment& al) const;

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

  vector<double> read_probs_;
  vector<Path> old_paths_;

};

class HICReadProbabilityCalculator {
  // @TODO implement HIC stuff
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
      bool use_as_advice
  ): read_set_(read_set), path_aligner_(read_set), //mismatch_prob_(mismatch_prob),
     insert_prob_(insert_prob), del_prob_(del_prob), subst_prob_(subst_prob),
     translocation_prob_(translocation_prob),
     min_prob_start_(min_prob_start), min_prob_per_base_(min_prob_per_base),
     penalty_constant_(penalty_constant), penalty_step_(penalty_step),
     old_paths_length_(1), binsize_(binsize), use_as_advice_(use_as_advice) {
    //read_probs_.resize(read_set_->size());
    read_cis_phis_.resize(read_set_->size());
    read_trans_psis_.resize(read_set_->size());
    total_log_prob_ = InitTotalLogProb();

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

  //vector<double> read_probs_;
  vector<double> read_cis_phis_;
  vector<double> read_trans_psis_;
  vector<Path> old_paths_;

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

  vector<SingleShortReadSet<>*> single_short_read_sets_;
  vector<ShortPairedReadSet<>*> paired_read_sets_;
  vector<HICReadSet<>*> hic_read_sets_;
  // (prob calculator, weight)
  // @TODO do you speak polymorphism? (vytvorit template na read calculatory) (i3)
  vector<pair<SingleReadProbabilityCalculator, double>> single_read_calculators_;
  vector<pair<PairedReadProbabilityCalculator, double>> paired_read_calculators_;
  vector<pair<HICReadProbabilityCalculator, double>> hic_read_calculators_;

 private:


};

#endif
