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



double pdf_normal(const double x, const double mean, const double stde) {
  const double inv_sqrt_2pi = 0.3989422804014327;
  const double a = (mean - x)/stde;
  return inv_sqrt_2pi / stde * exp(-0.5 * a * a);
}

inline double SimpleReadProbModel (int dist, int read_length, double mismatch_prob) {
  // @DONE corrected epsilon for e/3 (error from the GAML article)
  return pow(mismatch_prob/3, dist) * pow(1 - mismatch_prob, read_length - dist);
}

inline double BetterReadProbModel (int matches, int inserts, int dels, int substs, double insert_prob, double del_prob, double subst_prob) {
  const double res = pow(insert_prob, inserts) * pow(del_prob, dels) * pow(subst_prob, substs) * pow(1 - insert_prob - del_prob - subst_prob, matches);
  //if (rand() % 1000 == 0) cerr << "PROBMODEL: " << matches << "\t" << inserts << "\t" << dels << "\t" << substs << "\t" << res << endl;
  return res;
  //return exp(inserts * log(insert_prob) + dels * log(del_prob) + substs*log(subst_prob) + matches * (1 - insert_prob - del_prob - subst_prob));
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
  path_aligner_.RemovePathsFromCache(prob_change.removed_paths);
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
    ret += p.GetLength();
  }
  return ret;
}
int SingleReadProbabilityCalculator::GetUnalignedReadsCount() {
  int res = 0;
  for (int i = 0; i < read_set_->size(); i++) {
    if (log(max(0.0, read_probs_[i])) <= GetMinLogProbability((int)(*read_set_)[i].size())) res += 1;
  }
  return res;
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

  path_aligner_.RemovePathsFromCache(prob_change.removed_paths);
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

  return EvalTotalProbabilityFromChange(prob_change, false);
}
void PairedReadProbabilityCalculator::EvalProbabilityChange(PairedProbabilityChange &prob_change, bool debug_output) {
  for (size_t i = 0; i < prob_change.added_paths.size(); i++) {
    auto &p = prob_change.added_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.added_alignments.insert(prob_change.added_alignments.end(), als.begin(), als.end());
    if (debug_output) printf("\n(added_paths) done %d/%d evals; %d alignments", (int) i+1, (int) prob_change.added_paths.size(), (int) als.size());
    if (debug_output) fflush(stdout);
    // @DEBUG
    for (auto &al: als) {
      if (GetAlignmentProb(al) == 0) {
        for (auto &a: {al.al1, al.al2}) {
          cerr << "(" << a.matches << "," << a.substs << "," << a.inserts << "," << a.dels << "|"<< BetterReadProbModel(a.matches,a.inserts, a.dels, a.substs, insert_prob_, del_prob_, subst_prob_)<< ") ";
        }
        cerr << "orient: " << al.orientation << " ";
        cerr << " len:" << al.insert_length << "|" << pdf_normal(al.insert_length, mean_distance_, std_distance_) << endl;

      }
    }
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
  //new_prob += log(old_paths_length_);
  //new_prob -= log(prob_change.new_paths_length);

  // because we are not averageing anymore
  new_prob += log(2 * old_paths_length_) * read_set_->size();
  new_prob -= log(2 * prob_change.new_paths_length) * read_set_->size();

  map<string, int> orient_stats;

  // (read_id, prob_change)
  vector< pair<int, double> > changes;
  for (auto &a: prob_change.added_alignments) {
    //changes.push_back(make_pair(a.read_id, GetAlignmentProb(a)));
    changes.emplace_back(a.read_id, GetAlignmentProb(a));
    if (write) read_als_[a.read_id]++;
    orient_stats[a.orientation] = 1 + orient_stats[a.orientation];
  }

  for (auto &a: prob_change.removed_alignments) {
    //changes.push_back(make_pair(a.read_id, -GetAlignmentProb(a)));
    changes.emplace_back(a.read_id, -GetAlignmentProb(a));
    if (write) read_als_[a.read_id]--;
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
  cerr  << "ORIENTATION STATS GLOBAL: ";
  for (auto &o: { "FR", "RF","FF", "RR", ""}) {
    cerr << o << ": " << orient_stats[o] << " \t";
  }
  cerr << endl;
  return new_prob;
}

int PairedReadProbabilityCalculator::GetPathsLength(const vector<Path> &paths) const {
  // same as for SingleReadProbabilityCalculator
  // @TODO refactor later
  int ret = 0;
  for (auto &p: paths) {
    ret += p.GetLength();
  }
  return ret;
}


bool is_same_orientation (const string& o1, const string& o2) {
  if (o1 == o2) {
    return true;
  }
  else if ((o1 == "FF" || o1 == "RR") && (o2 == "FF" || o2 == "RR") ) {
    return true;
  }
  else if ((o1 == "RF" && o2 == "FR") || (o1 == "FR" && o2 == "RF")) {
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

  //if(rand() % 10000 == 0) cerr << "PAIRED: " << pp1 << "\t" << pp2 << "\t" <<  pdf_normal(al.insert_length, mean_distance_, std_distance_) << endl;

  return pp1 * pp2 * pdf_normal(al.insert_length, mean_distance_, std_distance_);
}
int PairedReadProbabilityCalculator::GetUnalignedReadsCount() {
  int res = 0;
  for (int i = 0; i < read_set_->size(); i++) {
    if (log(max(0.0, read_probs_[i])) <= GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size())) res += 1;
    //if (read_probs_[i] == 0) res++;
    //if (read_als_[i] == 0) res++;
  }
  return res;
}
int PairedReadProbabilityCalculator::GetCompletelyUnalignedReadsCount() {
  int res = 0;
  for (int i = 0; i < read_set_->size(); i++) {
    if (read_als_[i] == 0) res++;
  }
  return res;
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
            hic_reads.binsize(),
            hic_reads.min_prob_start(),
            hic_reads.min_prob_per_base(),
            hic_reads.penalty_constant(),
            hic_reads.penalty_step(),
            hic_reads.use_as_advice()
        ),
        hic_reads.weight()
    ));
  }
}

double GlobalProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, ProbabilityChanges& prob_changes) {
  double total_prob = 0;
  long long total_reads_length = 0;
  // @TODO add weighting of datasets based on their lengths
  prob_changes.single_read_changes.clear();
  for (auto &single_read_calculator: single_read_calculators_) {
    SingleProbabilityChange ch;
    double prob = single_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * single_read_calculator.second;
    prob_changes.single_read_changes.push_back(ch);

    total_reads_length += single_read_calculator.first.GetTotalReadsLength();
  }

  prob_changes.paired_read_changes.clear();
  for (auto &paired_read_calculator: paired_read_calculators_) {
    PairedProbabilityChange ch;
    double prob = paired_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * paired_read_calculator.second;
    prob_changes.paired_read_changes.push_back(ch);

    total_reads_length += paired_read_calculator.first.GetTotalReadsLength();
  }

  prob_changes.hic_read_changes.clear();
  for (auto &hic_read_calculator: hic_read_calculators_) {
    HICProbabilityChange ch;
    double prob = hic_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * hic_read_calculator.second;
    prob_changes.hic_read_changes.push_back(ch);

    total_reads_length += hic_read_calculator.first.GetTotalReadsLength();
  }

  total_prob += GetAprioriPathsLogProbability(paths);
  total_prob /= total_reads_length;
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
    const auto length = (long long)p.GetLength();
    prob += length;
  }
  return -prob * log(4);
}
string GlobalProbabilityCalculator::GetUnalignedReadsDebug() {
  stringstream res;
  for (auto pc: single_read_calculators_) {
    //res << "S:" << pc.first.GetUnalignedReadsCount() << "(" << pc.first.GetCompletelyUnalignedReadsCount() << ")/" << pc.first.GetReadsCount()<< "|AVG:" << pc.first.GetAvgReadProb() << ", ";
  }
  for (auto pc: paired_read_calculators_) {
    res << "P:" << pc.first.GetUnalignedReadsCount() << "(" << pc.first.GetCompletelyUnalignedReadsCount() << ")/" << pc.first.GetReadsCount() << "|AVG:" << pc.first.GetAvgReadProb() << ", ";
  }
  for (auto pc: hic_read_calculators_) {
    res << "H:" << pc.first.GetUnalignedReadsCount() << "(" << pc.first.GetCompletelyUnalignedReadsCount() << ")/" << pc.first.GetReadsCount() << "|AVG:" << pc.first.GetAvgReadProb() << ", ";
  }

  return res.str();
}

string GlobalProbabilityCalculator::GetUnalignedReadsLog() {
  stringstream res;
  for (auto pc: single_read_calculators_) {
    //res << pc.first.GetUnalignedReadsCount() << ", ";
  }
  for (auto pc: paired_read_calculators_) {
    res << pc.first.GetUnalignedReadsCount()  << ", " << pc.first.GetCompletelyUnalignedReadsCount() << ", ";
  }
  for (auto pc: hic_read_calculators_) {
    res << pc.first.GetUnalignedReadsCount()  << ", " << pc.first.GetCompletelyUnalignedReadsCount() << ", ";
  }

  return res.str();
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
    ret += p.GetLength();
  }
  return ret;
}
void HICReadProbabilityCalculator::CommitProbabilityChange(const HICProbabilityChange &prob_change) {
  EvalTotalProbabilityFromChange(prob_change, true);
  old_paths_ = prob_change.new_paths;
  old_paths_length_ = prob_change.new_paths_length;
  path_aligner_.RemovePathsFromCache(prob_change.removed_paths);
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
  double new_total_prob = 0;
  int counter = 0;
  // -------- CIS CONNECTIONS ------------------
  vector<double> new_cis_phis(read_cis_phis_);
  vector<tuple<vector<Path>, int>> params_cis;
  params_cis.emplace_back(prob_change.added_paths, +1);
  params_cis.emplace_back(prob_change.removed_paths, -1);
  for (const auto &param: params_cis) {
    counter = 0;
    const vector<Path>& paths = get<0>(param);
    const double sygn = get<1>(param);
    vector<SingleReadAlignment> als1, als2;
    for (const auto &p: paths) {
      counter += 1;
      const double lambda = path_aligner_.eval_lambda(p);
      if (lambda == 0) continue;

      als1 = path_aligner_.GetPartAlignmentsForPath(p, 0);
      if (als1.empty()) continue;
      als2 = path_aligner_.GetPartAlignmentsForPath(p, 1);
      if (als2.empty()) continue;
      auto it1 = als1.begin();
      auto it2 = als2.begin();


      int subcounter = 0;
      while (it1 != als1.end() && it2 != als2.end()) {
        while (it2 != als2.end() && it2->read_id < it1->read_id) it2++;
        if (it2 == als2.end()) {
          break;
        }
        else if (it2->read_id == it1->read_id) {
          int righties_count = 0;
          const int read_id = it1->read_id;
          while (it2 != als2.end() && it2->read_id == read_id) {
            righties_count += 1;
            it2++;
          }
          while (it1 != als1.end() && it1->read_id == read_id) {
            it2 -= righties_count;
            for (int i = 0; i < righties_count; i++) {
              new_cis_phis[read_id] += sygn * GetCisScore(*it1, (int)read_set_->reads_1_[read_id].size(), *it2, (int)read_set_->reads_2_[read_id].size(), lambda);
              if (write) read_cis_als_[read_id] += sygn;
              subcounter++;
              it2++;
            }
            it1++;
          }
        }
        else { //  (it2->read_id > read, (i.e. we didn't hit any good reads)
          while (it1 != als1.end() && it1->read_id < it2->read_id) it1++;
        }
      }

      cout << "cis " << (sygn == 1 ? "added ":"removed ") << counter << " " << subcounter << " aligments; lambda: "<< lambda <<  endl;
    }
  }


  //---------- TRANS CONNECTIONS ---------------
  vector<double> new_trans_psis(read_trans_psis_);

  //       (((((( INNER TRANS )))))))))))
  vector<tuple<vector<Path>, int>> params_trans_inner;
  params_trans_inner.emplace_back(prob_change.added_paths, +1);
  params_trans_inner.emplace_back(prob_change.removed_paths, -1);

  for (auto &param: params_trans_inner) {
    const auto &paths = get<0>(param);
    const int sygn = get<1>(param);
    vector<SingleReadAlignment> als1, als2;
    counter = 0;
    for (int i = 0; i < (int)paths.size(); i++) {
      counter++;
      als1 = path_aligner_.GetPartAlignmentsForPath(paths[i], 0);
      if (als1.empty()) {
        counter += (int)paths.size() - 1;
        continue;
      }
      for (int j = 0; j < (int)paths.size(); j++) {
        if (i == j) continue; // only for different paths
        als2 = path_aligner_.GetPartAlignmentsForPath(paths[j], 1);
        if (als2.empty()) {
          counter++;
          continue;
        }

        auto it1 = als1.begin();
        auto it2 = als2.begin();

        while (it1 != als1.end() && it2 != als2.end()) {
          while (it2 != als2.end() && it2->read_id < it1->read_id) it2++;
          if (it2 == als2.end()) {
            break;
          } else if (it2->read_id == it1->read_id) {
            int righties_count = 0;
            const int read_id = it1->read_id;
            while (it2 != als2.end() && it2->read_id == read_id) {
              righties_count += 1;
              it2++;
            }
            while (it1 != als1.end() && it1->read_id == read_id) {
              it2 -= righties_count;
              for (int k = 0; k < righties_count; k++) {
                new_trans_psis[read_id] += sygn * GetTransScore(*it1, *it2);
                if (write) read_trans_als_[read_id] += sygn;
                it2++;
              }
              it1++;
            }
          } else { //  (it2->read_id > read, (i.e. we didn't hit any good reads)
            while (it1 != als1.end() && it1->read_id < it2->read_id) it1++;
          }
        }
        if (counter % 1000 == 0) cout << "\rtrans inner " << (sygn == 1 ? "added ":"removed ") << counter << " / " << paths.size() * ((int)paths.size() - 1);
      }
    }
    cout << "\rtrans inner " << (sygn == 1 ? "added ":"removed ") << counter << " / " << paths.size() * ((int)paths.size() - 1) << endl;
  }

  //       ((((((( OUTER TRANS ))))))))))
  vector<tuple<vector<Path>, vector<Path>, int>> params_trans_outer;
  params_trans_outer.emplace_back(prob_change.still_paths, prob_change.removed_paths, -1);
  params_trans_outer.emplace_back( prob_change.removed_paths, prob_change.still_paths, -1);
  params_trans_outer.emplace_back(prob_change.still_paths, prob_change.added_paths, +1);
  params_trans_outer.emplace_back( prob_change.added_paths, prob_change.still_paths, +1);

  for (const auto &param: params_trans_outer) {
    const vector<Path> &paths1 = get<0>(param);
    const vector<Path> &paths2 = get<1>(param);
    if (paths1.empty() || paths2.empty()) continue;
    const int sygn = get<2>(param);

    counter = 0;
    vector<SingleReadAlignment> als1, als2;

    for (const auto &p1: paths1) {
      als1 = path_aligner_.GetPartAlignmentsForPath(p1, 0);
      if (als1.empty()) {
        counter += paths2.size();
        continue;
      }
      for (const auto &p2: paths2) {
        als2 = path_aligner_.GetPartAlignmentsForPath(p2, 1);
        if (als2.empty()) {
          counter++;
          continue;
        }
        auto it1 = als1.begin();
        auto it2 = als2.begin();

        while (it1 != als1.end() && it2 != als2.end()) {
          while (it2 != als2.end() && it2->read_id < it1->read_id) it2++;
          if (it2 == als2.end()) {
            break;
          } else if (it2->read_id == it1->read_id) {
            int righties_count = 0;
            const int read_id = it1->read_id;
            while (it2 != als2.end() && it2->read_id == read_id) {
              righties_count += 1;
              it2++;
            }
            while (it1 != als1.end() && it1->read_id == read_id) {
              it2 -= righties_count;
              for (int k = 0; k < righties_count; k++) {
                new_trans_psis[read_id] += sygn * GetTransScore(*it1, *it2);
                if (write) read_trans_als_[read_id] += sygn;
                it2++;
              }
              it1++;
            }
          } else { //  (it2->read_id > read, (i.e. we didn't hit any good reads)
            while (it1 != als1.end() && it1->read_id < it2->read_id) it1++;
          }
        }

        counter++;
        if (counter % 1000 == 0) cout << "\rtrans outer " << (sygn == 1 ? "added ":"removed ") << counter << " / " << paths1.size() * paths2.size();
      }
    }
    cout << "\rtrans outer " << (sygn == 1 ? "added ":"removed ") << counter << " / " << paths1.size() * paths2.size() << endl;
  }

  // eval new cis weighting constant
  double new_cis_constant_ = (1.00 - translocation_prob_) / prob_change.new_paths_length;
  cout << "New cis constant " << new_cis_constant_ << endl;

  // eval new trans_constant
  double new_trans_constant = 0;
  vector<double> multis;
  multis.reserve(prob_change.new_paths.size() * (prob_change.new_paths.size() - 1) / 2);
  for (int i = 0; i < (int)prob_change.new_paths.size() - 1; i++) {
    for (int j = i+1; j < (int)prob_change.new_paths.size(); j++) {
      multis.push_back(1.00 * prob_change.new_paths[i].GetLength() * prob_change.new_paths[j].GetLength());
    }
  }
  sort(multis.begin(), multis.end());
  for (auto x: multis) new_trans_constant += x;
  new_trans_constant = translocation_prob_ / new_trans_constant;

  cout << "New trans constant " << new_trans_constant << endl;

  counter = 0;
  for (int i = 0; i < (int)read_set_->size(); i++) {
    const double read_prob = new_cis_constant_ * new_cis_phis[i] + new_trans_constant * new_trans_psis[i];
    new_total_prob += GetRealReadProbability(read_prob, i);
    counter += 1;
    if (counter % 1000 == 0) cout << "\rEvaled total prob for " << counter << " reads";
  }
  cout << "\rEvaled total prob for " << counter << " reads" << endl;

  if (write) {
    read_cis_phis_ = new_cis_phis;
    read_trans_psis_ = new_trans_psis;
    total_log_prob_ = new_total_prob;
    cis_constant_ = new_cis_constant_;
    trans_constant_ = new_trans_constant;
  }

  return new_total_prob;
}
double HICReadProbabilityCalculator::GetRealReadProbability(double prob, int read_id) const {
  return max(log(max(0.0, prob)), GetMinLogProbability((*read_set_)[read_id].first.size() + (*read_set_)[read_id].second.size()));
}
int HICReadProbabilityCalculator::GetUnalignedReadsCount() {
  int res = 0;
  for (int i = 0; i < read_set_->size(); i++) {
    const double read_prob = cis_constant_ * read_cis_phis_[i] + trans_constant_ * read_trans_psis_[i];
    if (log(max(0.0, read_prob)) <= GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size())) res += 1;
    //if (read_cis_als_[i] + read_trans_als_[i] == 0) res++;
  }
  return res;
}
int HICReadProbabilityCalculator::GetCompletelyUnalignedReadsCount() {
  int res = 0;
  for (int i = 0; i < read_set_->size(); i++) {
   if (read_cis_als_[i] + read_trans_als_[i] == 0) res++;
  }
  return res;
}
