#include "read_probability_calculator.h"
#include <algorithm>
#include <cmath>
#include <cassert>

double SingleReadProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, SingleProbabilityChange& prob_change) {
  prob_change.added_paths.clear();
  prob_change.removed_paths.clear();
  prob_change.added_alignments.clear();
  prob_change.removed_alignments.clear();
  ComparePathSets(old_paths_, paths, prob_change.added_paths, prob_change.removed_paths);

  prob_change.new_paths_length = GetPathsLength(paths);
  prob_change.new_paths = paths;

  EvalProbabilityChange(prob_change);

  return EvalTotalProbabilityFromChange(prob_change);
}

void SingleReadProbabilityCalculator::EvalProbabilityChange(
    SingleProbabilityChange& prob_change) {
  for (size_t i = 0; i < prob_change.added_paths.size(); i++) {
    auto &p = prob_change.added_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.added_alignments.insert(prob_change.added_alignments.end(), als.begin(), als.end()); 
    printf("\rdone %d/%d evals", (int) i+1, (int) prob_change.added_paths.size()); 
    fflush(stdout);
  }
  for (size_t i = 0; i < prob_change.removed_paths.size(); i++) {
    auto &p = prob_change.removed_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.removed_alignments.insert(prob_change.removed_alignments.end(), als.begin(), als.end()); 
    printf("\rdone %d/%d evals", (int) i+1, (int) prob_change.removed_paths.size()); 
    fflush(stdout);
  }
  printf("\n");
}

double SingleReadProbabilityCalculator::EvalTotalProbabilityFromChange(
    const SingleProbabilityChange& prob_change, bool write) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  // (read_id, prob_change)
  vector<pair<int, double>> changes;
  for (auto &a: prob_change.added_alignments) {
    //changes.push_back(make_pair(a.read_id, GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
    //changes.emplace_back(a.read_id, GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size()));
    changes.emplace_back(a.read_id, GetAlignmentProb(a.matches, a.inserts, a.dels, a.substs));
  }
  for (auto &a: prob_change.removed_alignments) {
    //changes.push_back(make_pair(a.read_id, -GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
    //changes.emplace_back(a.read_id, -GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size()));
    changes.emplace_back(a.read_id, -GetAlignmentProb(a.matches, a.inserts, a.dels, a.substs));
  }
  sort(changes.begin(), changes.end());
  int last_read_id = -47;
  double accumulated_prob = 0;
  for (auto &ch: changes) {
    if (ch.first != last_read_id && last_read_id != -47) {
      //new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
      new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id);
      //new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
      new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id);
      if (write) {
        read_probs_[last_read_id] += accumulated_prob;
      }
      accumulated_prob = 0;
    }
    accumulated_prob += ch.second;
    last_read_id = ch.first;
  }
  if (last_read_id != -47) {
    //new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
    new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id);
    //new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
    new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id);
    if (write) {
      read_probs_[last_read_id] += accumulated_prob;
    }
  }
  if (write) total_log_prob_ = new_prob;
  return new_prob;
}

inline double SimpleReadProbModel (int dist, int read_length, double mismatch_prob) {
  // @DONE corrected epsilon for e/3 (error from the GAML article)
  return pow(mismatch_prob/3, dist) * pow(1 - mismatch_prob, read_length - dist);
}

inline double BetterReadProbModel (int matches, int inserts, int dels, int substs, double insert_prob, double del_prob, double subst_prob) {
  return pow(insert_prob, inserts) * pow(del_prob, dels) * pow(subst_prob, substs) * pow(1 - insert_prob - del_prob - subst_prob, matches);
}

double SingleReadProbabilityCalculator::GetAlignmentProbSimple(
    int dist, int read_length) const {
  return 0;
  //return SimpleReadProbModel(dist, read_length, mismatch_prob_);
}

double SingleReadProbabilityCalculator::GetAlignmentProb(
    int matches, int inserts, int dels, int substs) const {
  return BetterReadProbModel(matches, inserts, dels, substs, insert_prob_, del_prob_, subst_prob_);
}

void SingleReadProbabilityCalculator::CommitProbabilityChange(
    const SingleProbabilityChange &prob_change) {
  EvalTotalProbabilityFromChange(prob_change, true);
  old_paths_ = prob_change.new_paths;
  old_paths_length_ = prob_change.new_paths_length; 
}

double SingleReadProbabilityCalculator::InitTotalLogProb() {
  double ret = 0;
  for (size_t i = 0; i < read_set_->size(); i++) {
    read_probs_[i] = 0;
    //ret += GetMinLogProbability((*read_set_)[i].size()) / read_set_->size();
    ret += GetMinLogProbability((*read_set_)[i].size());
  }
  return ret;
}

double SingleReadProbabilityCalculator::GetMinLogProbability(int read_length) const {
  return min_prob_start_ + read_length * min_prob_per_base_; 
}

double SingleReadProbabilityCalculator::GetRealReadProbability(double prob, int read_id) const {
  return max(log(max(0.0, prob)), GetMinLogProbability((*read_set_)[read_id].size()));
}

int SingleReadProbabilityCalculator::GetPathsLength(const vector<Path>& paths) const {
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}

double PairedReadProbabilityCalculator::InitTotalLogProb() {
  double ret = 0;
  for (size_t i = 0; i < read_set_->size(); i++) {
    read_probs_[i] = 0;
    //ret += GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size()) / read_set_->size();
    ret += GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size());
  }
  return ret;
}
void PairedReadProbabilityCalculator::CommitProbabilityChange(const PairedProbabilityChange &prob_change) {
  EvalTotalProbabilityFromChange(prob_change, true);
  old_paths_ = prob_change.new_paths;
  old_paths_length_ = prob_change.new_paths_length;

}
double PairedReadProbabilityCalculator::GetPathsProbability(const vector<Path> &paths,
                                                            PairedProbabilityChange &prob_change) {
  prob_change.added_paths.clear();
  prob_change.removed_paths.clear();
  prob_change.added_alignments.clear();
  prob_change.removed_alignments.clear();
  ComparePathSets(old_paths_, paths, prob_change.added_paths, prob_change.removed_paths);

  prob_change.new_paths_length = GetPathsLength(paths);
  prob_change.new_paths = paths;

  EvalProbabilityChange(prob_change);

  return EvalTotalProbabilityFromChange(prob_change);
}
void PairedReadProbabilityCalculator::EvalProbabilityChange(PairedProbabilityChange &prob_change, bool debug_output) {
  for (size_t i = 0; i < prob_change.added_paths.size(); i++) {
    auto &p = prob_change.added_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.added_alignments.insert(prob_change.added_alignments.end(), als.begin(), als.end());
    if (debug_output) printf("\n(added_paths) done %d/%d evals; %d alignments", (int) i+1, (int) prob_change.added_paths.size(), (int) als.size());
    if (debug_output) fflush(stdout);
  }
  if (debug_output) printf("\n");
  for (size_t i = 0; i < prob_change.removed_paths.size(); i++) {
    auto &p = prob_change.removed_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.removed_alignments.insert(prob_change.removed_alignments.end(), als.begin(), als.end());
    if (debug_output) printf("\n(removed_paths) done %d/%d evals;  %d alignments", (int) i+1, (int) prob_change.removed_paths.size(), (int) als.size());
    if (debug_output) fflush(stdout);
  }
  if (debug_output) printf("\n");
}

double PairedReadProbabilityCalculator::GetMinLogProbability(int read_length) const {
  // see "Reads that have no good alignment to A" in the  GAML article
  return min_prob_start_ + read_length * min_prob_per_base_;
}

double PairedReadProbabilityCalculator::GetRealReadProbability(double prob, int read_id) const {
  return max(log(max(0.0, prob)), GetMinLogProbability((*read_set_)[read_id].first.size() + (*read_set_)[read_id].second.size()));
}

double PairedReadProbabilityCalculator::EvalTotalProbabilityFromChange(const PairedProbabilityChange &prob_change,
                                                                       bool write) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  // (read_id, prob_change)
  vector< pair<int, double> > changes;
  for (auto &a: prob_change.added_alignments) {
    //changes.push_back(make_pair(a.read_id, GetAlignmentProb(a)));
    changes.emplace_back(a.read_id, GetAlignmentProb(a));
  }

  for (auto &a: prob_change.removed_alignments) {
    //changes.push_back(make_pair(a.read_id, -GetAlignmentProb(a)));
    changes.emplace_back(a.read_id, -GetAlignmentProb(a));
  }

  sort(changes.begin(), changes.end());

  const int BLANK = -47;
  int last_read_id = BLANK;
  double accumulated_prob = 0;
  for (auto &ch: changes) {
    if (ch.first != last_read_id && last_read_id != -47) {
      //new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
      new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id);
      //new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
      new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id);
      if (write) {
        read_probs_[last_read_id] += accumulated_prob;
      }
      accumulated_prob = 0;
    }
    accumulated_prob += ch.second;
    last_read_id = ch.first;
  }
  if (last_read_id != -47) {
    //new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
    new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id);
    //new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
    new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id);
    if (write) {
      read_probs_[last_read_id] += accumulated_prob;
    }
  }
  if (write) total_log_prob_ = new_prob;
  return new_prob;
}

int PairedReadProbabilityCalculator::GetPathsLength(const vector<Path> &paths) const {
  // same as for SingleReadProbabilityCalculator
  // @TODO refactor later
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}

double pdf_normal(const double x, const double mean, const double std) {
  static const double inv_sqrt_2pi = 0.3989422804014327;
  const double a = (mean - x)/std;
  return inv_sqrt_2pi / std * exp(-0.5 * a * a);
}

bool is_same_orientation (const string& o1, const string& o2) {
  if (o1 == o2) {
    return true;
  }
  else if ((o1 == "FF" || o1 == "RR") && (o2 == "FF" || o2 == "RR") ) {
    return true;
  }
  else {
    return false;
  }
}

double PairedReadProbabilityCalculator::GetAlignmentProb(const PairedReadAlignment& al) const {
  if (! is_same_orientation(al.orientation, read_set_->orientation_)) {return 0;}
  //const double pp1 = SimpleReadProbModel(al.al1.dist, (*read_set_)[al.read_id].first.size(), mismatch_prob_);
  const double pp1 = BetterReadProbModel(al.al1.matches, al.al1.inserts, al.al1.dels, al.al1.substs, insert_prob_, del_prob_, subst_prob_);
  //const double pp2 = SimpleReadProbModel(al.al2.dist, (*read_set_)[al.read_id].second.size(), mismatch_prob_);
  const double pp2 = BetterReadProbModel(al.al2.matches, al.al2.inserts, al.al2.dels, al.al2.substs, insert_prob_, del_prob_, subst_prob_);

  return pp1 * pp2 * pdf_normal(al.insert_length, mean_distance_, std_distance_);
}

GlobalProbabilityCalculator::GlobalProbabilityCalculator(const Config& config) {
  for (auto &single_short_reads: config.single_short_reads()) {
    SingleShortReadSet<>* rs = new SingleShortReadSet<>();
    rs->LoadReadSet(single_short_reads.filename());
    single_short_read_sets_.push_back(rs);
    single_read_calculators_.push_back(make_pair(SingleReadProbabilityCalculator(
          rs,
          //single_short_reads.mismatch_prob(),
          single_short_reads.insert_prob(),
          single_short_reads.del_prob(),
          single_short_reads.subst_prob(),
          single_short_reads.min_prob_start(),
          single_short_reads.min_prob_per_base(),
          single_short_reads.penalty_constant(),
          single_short_reads.penalty_step()), single_short_reads.weight()));
  }


  for (auto &paired_reads: config.paired_reads()) {
    ShortPairedReadSet<>* rs = new ShortPairedReadSet<>();
    rs->LoadReadSet(paired_reads.filename1(), paired_reads.filename2(), paired_reads.orientation());
    paired_read_sets_.push_back(rs);
    paired_read_calculators_.push_back(make_pair(
        PairedReadProbabilityCalculator(
            rs,
            //paired_reads.mismatch_prob(),
            paired_reads.insert_prob(),
            paired_reads.del_prob(),
            paired_reads.subst_prob(),
            paired_reads.min_prob_start(),
            paired_reads.min_prob_per_base(),
            paired_reads.penalty_constant(),
            paired_reads.penalty_step(),
            paired_reads.mean_distance(),
            paired_reads.std_distance(),
            paired_reads.use_as_advice()
        ),
        paired_reads.weight()
    ));
  }

  for (auto &hic_reads: config.hic_reads()) {
    HICReadSet<>* rs = new HICReadSet<>();
    rs->LoadReadSet(hic_reads.filename1(), hic_reads.filename2());
    hic_read_sets_.push_back(rs);
    hic_read_calculators_.push_back(make_pair(
        HICReadProbabilityCalculator(
            rs,
            hic_reads.insert_prob(),
            hic_reads.del_prob(),
            hic_reads.subst_prob(),
            hic_reads.translocation_prob(),
            hic_reads.min_prob_start(),
            hic_reads.min_prob_per_base(),
            hic_reads.penalty_constant(),
            hic_reads.penalty_step(),
            hic_reads.binsize(),
            hic_reads.use_as_advice()
        ),
        hic_reads.weight()
    ));
  }
}

double GlobalProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, ProbabilityChanges& prob_changes) {
  double total_prob = 0;
  // @TODO add weighting of datasets based on their lengths
  prob_changes.single_read_changes.clear();
  for (auto &single_read_calculator: single_read_calculators_) {
    SingleProbabilityChange ch;
    double prob = single_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * single_read_calculator.second;
    prob_changes.single_read_changes.push_back(ch);
  }

  prob_changes.paired_read_changes.clear();
  for (auto &paired_read_calculator: paired_read_calculators_) {
    PairedProbabilityChange ch;
    double prob = paired_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * paired_read_calculator.second;
    prob_changes.paired_read_changes.push_back(ch);
  }

  prob_changes.hic_read_changes.clear();
  for (auto &hic_read_calculator: hic_read_calculators_) {
    HICProbabilityChange ch;
    double prob = hic_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * hic_read_calculator.second;
    prob_changes.hic_read_changes.push_back(ch);
  }

  total_prob += GetAprioriPathsLogProbability(paths);
  return total_prob;
}

void GlobalProbabilityCalculator::CommitProbabilityChanges(
    const ProbabilityChanges &prob_changes) {
  assert(prob_changes.single_read_changes.size() == single_read_calculators_.size());
  for (size_t i = 0; i < single_read_calculators_.size(); i++) {
    single_read_calculators_[i].first.CommitProbabilityChange(prob_changes.single_read_changes[i]);
  }

  assert(prob_changes.paired_read_changes.size() == paired_read_calculators_.size());
  for (size_t i = 0; i < paired_read_calculators_.size(); i++) {
    paired_read_calculators_[i].first.CommitProbabilityChange(prob_changes.paired_read_changes[i]);
  }

  assert(prob_changes.hic_read_changes.size() == hic_read_calculators_.size());
  for (size_t i = 0; i < hic_read_calculators_.size(); i++) {
    hic_read_calculators_[i].first.CommitProbabilityChange(prob_changes.hic_read_changes[i]);
  }
}
double GlobalProbabilityCalculator::GetAprioriPathsLogProbability(const vector<Path> &paths) {
  int prob = 0;
  for (auto &p: paths) {
    const auto length = (long long)p.ToString(true).size();
    prob += length;
  }
  return -prob * log(4);
}



double HICReadProbabilityCalculator::InitTotalLogProb() {
  double ret = 0;
  for (size_t i = 0; i < read_set_->size(); i++) {
    //read_probs_[i] = 0;
    //ret += GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size()) / read_set_->size();
    ret += GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size());
    read_cis_phis_[i] = 0;
    read_trans_psis_[i] = 0;
  }
  return ret;
}
double HICReadProbabilityCalculator::GetMinLogProbability(int read_length) const {
  // see "Reads that have no good alignment to A" in the  GAML article
  return min_prob_start_ + read_length * min_prob_per_base_;
}
double HICReadProbabilityCalculator::GetPathsProbability(const vector<Path> &paths, HICProbabilityChange &prob_change) {
  prob_change.added_paths.clear();
  prob_change.removed_paths.clear();
  prob_change.still_paths.clear();

  ComparePathSetsWithStill(old_paths_, paths, prob_change.added_paths, prob_change.removed_paths, prob_change.still_paths);

  prob_change.new_paths_length = GetPathsLength(paths);
  prob_change.new_paths = paths;

  // skip optimisation step
  //EvalProbabilityChange(prob_change);

  return EvalTotalProbabilityFromChange(prob_change, false);
}
int HICReadProbabilityCalculator::GetPathsLength(const vector<Path> &paths) const {
  // same as for SingleReadProbabilityCalculator
  // @TODO refactor later
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}
void HICReadProbabilityCalculator::CommitProbabilityChange(const HICProbabilityChange &prob_change) {
  EvalTotalProbabilityFromChange(prob_change, true);
  old_paths_ = prob_change.new_paths;
  old_paths_length_ = prob_change.new_paths_length;
}
void HICReadProbabilityCalculator::EvalProbabilityChange(HICProbabilityChange &prob_change, bool debug_output) {
  // @TODO implement
}

double exp_pdf(const int x, const double lambda) {
  if (x < 0) return 0;
  else return lambda * exp(-lambda * x);
}

double HICReadProbabilityCalculator::GetCisScore(const SingleReadAlignment al1, const int a1_len, const SingleReadAlignment al2, const int a2_len, const double lambda) {
  string orientation;
  int insert_length;
  tie(orientation, insert_length) = eval_orientation(al1, a1_len, al2, a2_len);
  return BetterReadProbModel(al1.matches, al1.inserts, al1.dels, al1.substs, insert_prob_, del_prob_, subst_prob_) *
      BetterReadProbModel(al2.matches, al2.inserts, al2.dels, al2.substs, insert_prob_, del_prob_, subst_prob_) *
      exp_pdf(insert_length/binsize_, lambda);
}

double HICReadProbabilityCalculator::GetTransScore (const SingleReadAlignment al1, const SingleReadAlignment al2) {
  return BetterReadProbModel(al1.matches, al1.inserts, al1.dels, al1.substs, insert_prob_, del_prob_, subst_prob_) *
      BetterReadProbModel(al2.matches, al2.inserts, al2.dels, al2.substs, insert_prob_, del_prob_, subst_prob_);
}

double HICReadProbabilityCalculator::EvalTotalProbabilityFromChange(const HICProbabilityChange &prob_change,
                                                                    bool write) {
  // cis part
  // revoke old lambda from aligner
  // remove old phis using old lambda

  vector<double> new_cis_phis(read_cis_phis_);
  vector<SingleReadAlignment> lefties, righties;
  //vector<SingleReadAlignment> als_left, als_right;
  int counter = 0;
  for (auto &p: prob_change.removed_paths) {
    const double old_lambda = path_aligner_.eval_lambda(p);
    const auto &als_left = path_aligner_.GetPartAlignmentsForPath(p, 0);
    const auto &als_right = path_aligner_.GetPartAlignmentsForPath(p, 1);

    // assume alignments are sorted
    auto it_al1 = als_left.begin();
    auto it_al2 = als_right.begin();

    for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
      while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
      while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
      lefties.clear();
      righties.clear();
      while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
        lefties.push_back(*it_al1);
        it_al1++;
      }
      while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
        righties.push_back(*it_al2);
        it_al2++;
      }

      if (!lefties.empty() && !righties.empty()) {
        for (auto al1: lefties) {
          for (auto al2: righties) {
            new_cis_phis[read_id] -= GetCisScore(al1, read_set_->reads_1_[read_id].size(), al2, read_set_->reads_2_[read_id].size(), old_lambda);
          }
        }
      }
    }
    counter += 1;
    cout << "(removed cis path) done " << counter << " / " << prob_change.removed_paths.size() << endl;
  }

  // eval new lambda by aligner
  // add new phis using new lambda
  counter = 0;
  for (auto &p: prob_change.added_paths) {
    const double new_lambda = path_aligner_.eval_lambda(p);
    const auto &als_left = path_aligner_.GetPartAlignmentsForPath(p, 0);
    const auto &als_right = path_aligner_.GetPartAlignmentsForPath(p, 1);

    // assume alignments are sorted
    auto it_al1 = als_left.begin();
    auto it_al2 = als_right.begin();
    //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
    for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
      while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
      while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
      lefties.clear();
      righties.clear();
      while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
        lefties.push_back(*it_al1);
        it_al1++;
      }
      while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
        righties.push_back(*it_al2);
        it_al2++;
      }

      if (!lefties.empty() && !righties.empty()) {
        for (auto al1: lefties) {
          for (auto al2: righties) {
            new_cis_phis[read_id] += GetCisScore(al1, read_set_->reads_1_[read_id].size(), al2, read_set_->reads_2_[read_id].size(), new_lambda);
          }
        }
      }
    }
    counter += 1;
    if (counter % 100 == 0) cout << "(added cis paths) done: " << counter << "/" << prob_change.added_paths.size() << endl;
  }

  // eval new cis weighting constant
  double new_cis_constant_ = (1.00 - translocation_prob_) / prob_change.new_paths_length;

  // trans part
  vector<double> new_trans_psis(read_trans_psis_);
  // remove old stuff
  // \Psi(r, R)
  counter = 0;
  for (int i = 0; i < (int)prob_change.removed_paths.size(); i++) {
    for (int j = 0; j < (int)prob_change.removed_paths.size(); j++) {
      if (i == j) continue; // only different paths
      auto &pl = prob_change.removed_paths[i];
      auto &pr = prob_change.removed_paths[j];

      const auto &als_left = path_aligner_.GetPartAlignmentsForPath(pl, 0);
      const auto &als_right = path_aligner_.GetPartAlignmentsForPath(pr, 1);

      // assume alignments are sorted
      auto it_al1 = als_left.begin();
      auto it_al2 = als_right.begin();
      //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
      for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
        while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
        while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
        lefties.clear();
        righties.clear();
        while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
          lefties.push_back(*it_al1);
          it_al1++;
        }
        while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
          righties.push_back(*it_al2);
          it_al2++;
        }

        if (!lefties.empty() && !righties.empty()) {
          for (auto al1: lefties) {
            for (auto al2: righties) {
              new_trans_psis[read_id] -= GetTransScore(al1, al2);
            }
          }
        }
      }
      counter += 1;
      if (counter % 100 == 0) cout << "(removed trans paths) done: " << counter << " / " << prob_change.removed_paths.size() * (prob_change.removed_paths.size() -1) << endl;
    }
  }
  // psi(s_r, s_k)
  for (int i = 0; i < (int)prob_change.removed_paths.size(); i++) {
    for (int j = 0; j < (int)prob_change.still_paths.size(); j++) {
      auto &pl = prob_change.removed_paths[i];
      auto &pr = prob_change.still_paths[j];

      const auto &als_left = path_aligner_.GetPartAlignmentsForPath(pl, 0);
      const auto &als_right = path_aligner_.GetPartAlignmentsForPath(pr, 1);

      // assume alignments are sorted
      auto it_al1 = als_left.begin();
      auto it_al2 = als_right.begin();
      //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
      for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
        while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
        while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
        lefties.clear();
        righties.clear();
        while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
          lefties.push_back(*it_al1);
          it_al1++;
        }
        while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
          righties.push_back(*it_al2);
          it_al2++;
        }

        if (!lefties.empty() && !righties.empty()) {
          for (auto al1: lefties) {
            for (auto al2: righties) {
              new_trans_psis[read_id] -= GetTransScore(al1, al2);
            }
          }
        }
      }
      counter += 1;
      if (counter % 100 == 0) cout << "(removed trans R-K paths) done: " << counter << " / " << prob_change.removed_paths.size() * (prob_change.still_paths.size()) << endl;
    }
  }
  // psi(s_k, s_s)
  counter = 0;
  for (int i = 0; i < (int)prob_change.still_paths.size(); i++) {
    for (int j = 0; j < (int)prob_change.removed_paths.size(); j++) {
      auto &pl = prob_change.still_paths[i];
      auto &pr = prob_change.removed_paths[j];

      const auto &als_left = path_aligner_.GetPartAlignmentsForPath(pl, 0);
      const auto &als_right = path_aligner_.GetPartAlignmentsForPath(pr, 1);

      // assume alignments are sorted
      auto it_al1 = als_left.begin();
      auto it_al2 = als_right.begin();
      //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
      for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
        while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
        while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
        lefties.clear();
        righties.clear();
        while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
          lefties.push_back(*it_al1);
          it_al1++;
        }
        while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
          righties.push_back(*it_al2);
          it_al2++;
        }

        if (!lefties.empty() && !righties.empty()) {
          for (auto al1: lefties) {
            for (auto al2: righties) {
              new_trans_psis[read_id] -= GetTransScore(al1, al2);
            }
          }
        }
      }
      counter += 1;
      if (counter % 100 == 0) cout << "(removed trans K-R paths) done: " << counter << " / " << prob_change.removed_paths.size() * (prob_change.still_paths.size()) << endl;
    }
  }

  // add new stuff
  // Psi(A)
  counter = 0;
  for (int i = 0; i < (int)prob_change.added_paths.size(); i++) {
    for (int j = 0; j < (int)prob_change.added_paths.size(); j++) {
      if (i == j) continue; // only different paths
      auto &pl = prob_change.added_paths[i];
      auto &pr = prob_change.added_paths[j];

      const auto &als_left = path_aligner_.GetPartAlignmentsForPath(pl, 0);
      const auto &als_right = path_aligner_.GetPartAlignmentsForPath(pr, 1);

      // assume alignments are sorted
      auto it_al1 = als_left.begin();
      auto it_al2 = als_right.begin();
      //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
      for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
        while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
        while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
        lefties.clear();
        righties.clear();
        while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
          lefties.push_back(*it_al1);
          it_al1++;
        }
        while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
          righties.push_back(*it_al2);
          it_al2++;
        }

        if (!lefties.empty() && !righties.empty()) {
          for (auto al1: lefties) {
            for (auto al2: righties) {
              new_trans_psis[read_id] += GetTransScore(al1, al2);
            }
          }
        }
      }
      counter += 1;
      if (counter % 100 == 0) cout << "(added trans paths) done: " << counter << " / " << prob_change.added_paths.size() * (prob_change.added_paths.size() -1) << endl;
    }
  }

  // psi(s_a, s_k)
  counter = 0;
  for (int i = 0; i < (int)prob_change.added_paths.size(); i++) {
    for (int j = 0; j < (int)prob_change.still_paths.size(); j++) {
      auto &pl = prob_change.added_paths[i];
      auto &pr = prob_change.still_paths[j];

      const auto &als_left = path_aligner_.GetPartAlignmentsForPath(pl, 0);
      const auto &als_right = path_aligner_.GetPartAlignmentsForPath(pr, 1);

      // assume alignments are sorted
      auto it_al1 = als_left.begin();
      auto it_al2 = als_right.begin();
      //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
      for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
        while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
        while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
        lefties.clear();
        righties.clear();
        while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
          lefties.push_back(*it_al1);
          it_al1++;
        }
        while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
          righties.push_back(*it_al2);
          it_al2++;
        }

        if (!lefties.empty() && !righties.empty()) {
          for (auto al1: lefties) {
            for (auto al2: righties) {
              new_trans_psis[read_id] += GetTransScore(al1, al2);
            }
          }
        }
      }
      counter += 1;
      if (counter % 100 == 0) cout << "(added trans A-K paths) done: " << counter << " / " << prob_change.added_paths.size() * (prob_change.still_paths.size()) << endl;
    }
  }
  // psi(s_a, s_s)
  counter = 0;
  for (int i = 0; i < (int)prob_change.still_paths.size(); i++) {
    for (int j = 0; j < (int)prob_change.added_paths.size(); j++) {
      auto &pl = prob_change.still_paths[i];
      auto &pr = prob_change.added_paths[j];

      const auto &als_left = path_aligner_.GetPartAlignmentsForPath(pl, 0);
      const auto &als_right = path_aligner_.GetPartAlignmentsForPath(pr, 1);

      // assume alignments are sorted
      auto it_al1 = als_left.begin();
      auto it_al2 = als_right.begin();
      //vector<SingleReadAlignment> lefties, righties; // we'll just use old ones
      for (int read_id = 0; read_id < (int)read_set_->size(); read_id++) {
        while (it_al1 != als_left.end() && it_al1->read_id < read_id) it_al1++;
        while (it_al2 != als_right.end() && it_al2->read_id < read_id) it_al2++;
        lefties.clear();
        righties.clear();
        while (it_al1 != als_left.end() && it_al1->read_id == read_id) {
          lefties.push_back(*it_al1);
          it_al1++;
        }
        while (it_al2 != als_right.end() && it_al2->read_id == read_id) {
          righties.push_back(*it_al2);
          it_al2++;
        }

        if (!lefties.empty() && !righties.empty()) {
          for (auto al1: lefties) {
            for (auto al2: righties) {
              new_trans_psis[read_id] += GetTransScore(al1, al2);
            }
          }
        }
      }
      counter += 1;
      if (counter % 100 == 0) cout << "(added trans K-A paths) done: " << counter << " / " << prob_change.added_paths.size() * (prob_change.still_paths.size()) << endl;
    }
  }

  // eval new trans_constant
  double new_trans_constant = 0;
  vector<double> multis;
  multis.reserve(prob_change.new_paths.size() * (prob_change.new_paths.size() - 1) / 2);
  for (int i = 0; i < (int)prob_change.new_paths.size() - 1; i++) {
    for (int j = i+1; j < (int)prob_change.new_paths.size(); j++) {
      multis.push_back(1.00 * prob_change.new_paths[i].ToString(true).size() * prob_change.new_paths[j].ToString(true).size());
    }
  }
  sort(multis.begin(), multis.end());
  for (auto x: multis) new_trans_constant += x;
  new_trans_constant = translocation_prob_ / new_trans_constant;

  double new_total_prob = 0;

  for (int i = 0; i < (int)read_set_->size(); i++) {
    const double read_prob = new_cis_constant_ * new_cis_phis[i] + new_trans_constant * new_trans_psis[i];
    new_total_prob += GetRealReadProbability(read_prob, i);
  }

  if (write) {
    read_cis_phis_ = new_cis_phis;
    read_trans_psis_ = new_trans_psis;
    total_log_prob_ = new_total_prob;
  }

  return new_total_prob;
}
double HICReadProbabilityCalculator::GetRealReadProbability(double prob, int read_id) const {
  return max(log(max(0.0, prob)), GetMinLogProbability((*read_set_)[read_id].first.size() + (*read_set_)[read_id].second.size()));
}
