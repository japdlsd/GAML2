#ifndef PATH_ALIGNER_H__
#define PATH_ALIGNER_H__

#include <algorithm>
#include <unordered_map>
#include <map>
#include "read_set.h"
#include "hash_util.h"
#include "path.h"


// Handles caching of alignments


class SingleShortReadPathAligner {
 public:
  SingleShortReadPathAligner() {}
  explicit SingleShortReadPathAligner(SingleShortReadSet<>* single_short_read_set): single_short_read_set_(single_short_read_set) {}
  vector<SingleReadAlignment> GetAlignmentsForPath(const Path& p) {
    auto it = cache_.find(p);
    if (it != cache_.end()) return it->second;
    auto res = single_short_read_set_->GetAlignments(p.ToString(true));
    sort(res.begin(), res.end());
    cache_[p] = res;
    return res;
  };
  void RemovePathFromCache(const Path& p) {
    cache_.erase(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    for (auto p: paths) {
      cache_.erase(p);
    }
  }

  void GetAlignmentsForPathByReference(const Path& p, vector<SingleReadAlignment>& output) {
    auto it = cache_.find(p);
    output.clear();

    if (it != cache_.end()) {
      output.insert(output.end(), it->second.begin(), it->second.end());
    }
    else {
      auto res = single_short_read_set_->GetAlignments(p.ToString(true));
      sort(res.begin(), res.end());
      cache_[p] = res;
      output.insert(output.end(), res.begin(), res.end());
    }
  }

  SingleShortReadSet<>* single_short_read_set_;
 private:
  unordered_map<Path, vector<SingleReadAlignment>> cache_;

};

class SingleShortReadPathAlignerVector {
 public:
  explicit SingleShortReadPathAlignerVector() {}
  explicit SingleShortReadPathAlignerVector(SingleShortReadSet<>* single_short_read_set) : single_short_read_set_(single_short_read_set) {}
  //explicit SingleShortReadPathAligner(SingleShortReadSet<StandardReadIndex>* single_short_read_set) : single_short_read_set_(single_short_read_set) {}


  vector<SingleReadAlignment> GetAlignmentsForPath(const Path& p);

  SingleShortReadSet<>* single_short_read_set_;

 private:

  vector<SingleReadAlignment> GetAlignmentForPathNoCache(const Path& p);
  vector<SingleReadAlignment> GetAlignmentForPathWithCache(const Path& p);

  const static int MAX_CACHE_SIZE = 10000;
  long long next_usage_t_ = 0;
  vector<tuple<Path, vector<SingleReadAlignment>, long long>> cache_;
  void InsertAlignmentForPath(const Path& p, vector<SingleReadAlignment>& al) {
    int pos = GetAlignmentPos(p);
    if (pos == -1) {
      if (cache_.size() >= MAX_CACHE_SIZE) RemoveAllAlignments();
      sort(al.begin(), al.end());
      cache_.push_back(make_tuple(p, al, next_usage_t_++));
    }
  }


  void RemoveAllAlignments() {
    cache_.clear();
    next_usage_t_ = 0;
  }

  int GetAlignmentPos(const Path& p) const{
    int res = -1;
    for (int i = 0; i < (int)cache_.size(); i++) {
      const auto t = cache_[i];
      if (get<0>(t).IsSameNoReverse(p)) {
        res = i;
        break;
      }
    }
    return res;
  }

  vector<SingleReadAlignment> GetCachedAlignmentByPos(int pos) {
    assert(0 <= pos && pos < (int)cache_.size());
    auto &t = cache_[pos];
    get<2>(t) = next_usage_t_++;
    return get<1>(t);
  }
};

class PairedReadPathAligner {
 public:
  explicit PairedReadPathAligner() {}
  explicit PairedReadPathAligner(ShortPairedReadSet<>* paired_read_set): paired_read_set_(paired_read_set) {
    left_aligner_ = SingleShortReadPathAligner(&(paired_read_set_->reads_1_));
    right_aligner_ = SingleShortReadPathAligner(&(paired_read_set_->reads_2_));
  }

  vector<PairedReadAlignment> GetAlignmentsForPath(const Path& p);
  // part = 0 or 1 (0 for left, 1 for right)
  vector<SingleReadAlignment> GetPartAlignmentsForPath(const Path& p, int part);

  void RemovePathFromCache(const Path& p) {
    left_aligner_.RemovePathFromCache(p);
    right_aligner_.RemovePathFromCache(p);
    cache_.erase(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    left_aligner_.RemovePathsFromCache(paths);
    right_aligner_.RemovePathsFromCache(paths);
    for (auto &p: paths) {
      cache_.erase(p);
    }
  }

  ShortPairedReadSet<>* paired_read_set_;

  SingleShortReadPathAligner left_aligner_, right_aligner_;


 private:
  unordered_map<Path, vector<PairedReadAlignment>> cache_;
  //unordered_map<Path, int> cache_uncovered_bases_;


};

/*
class PairedReadPathAligner {
 public:
  explicit PairedReadPathAligner() {}
  explicit PairedReadPathAligner(ShortPairedReadSet<>* paired_read_set): paired_read_set_(paired_read_set) {
    left_aligner_ = SingleShortReadPathAligner(&(paired_read_set_->reads_1_));
    right_aligner_ = SingleShortReadPathAligner(&(paired_read_set_->reads_2_));
  }

  vector<PairedReadAlignment> GetAlignmentsForPath(const Path& p);
  // part = 0 or 1 (0 for left, 1 for right)
  vector<SingleReadAlignment> GetPartAlignmentsForPath(const Path& p, int part);


  void RemovePathFromCache(const Path& p) {
    left_aligner_.RemovePathFromCache(p);
    right_aligner_.RemovePathFromCache(p);
    // @TODO dopisat to pre vrchnu uroven cache
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    left_aligner_.RemovePathsFromCache(paths);
    right_aligner_.RemovePathsFromCache(paths);
    // @TODO dopisat to pre vrchnu uroven cache
  }

  ShortPairedReadSet<>* paired_read_set_;

  SingleShortReadPathAligner left_aligner_, right_aligner_;


 private:
  vector<PairedReadAlignment> GetAlignmentForPathNoCache(const Path& p);
  vector<PairedReadAlignment> GetAlignmentForPathWithCache(const Path& p);

  const static int MAX_CACHE_SIZE = 10000;
  long long next_usage_t_ = 0;
  vector<tuple<Path, vector<PairedReadAlignment>, long long>> cache_;
  void InsertAlignmentForPath(const Path& p, vector<PairedReadAlignment>& al) {
    int pos = GetAlignmentPos(p);
    if (pos == -1) {
      if (cache_.size() >= MAX_CACHE_SIZE) RemoveAllAlignments();
      sort(al.begin(), al.end());
      cache_.push_back(make_tuple(p, al, next_usage_t_++));
    }
  }

  void RemoveAllAlignments() {
    cache_.clear();
    next_usage_t_ = 0;
  }

  int GetAlignmentPos(const Path& p) const{
    int res = -1;
    for (int i = 0; i < (int)cache_.size(); i++) {
      const auto t = cache_[i];
      if (get<0>(t).IsSameNoReverse(p)) {
        res = i;
        break;
      }
    }
    return res;
  }

  vector<PairedReadAlignment> GetCachedAlignmentByPos(int pos) {
    assert(0 <= pos && pos < (int)cache_.size());
    auto &t = cache_[pos];
    get<2>(t) = next_usage_t_++;
    return get<1>(t);
  }

};*/
class HICReadPathAligner {
  // @TODO implement HIC stuff
 public:
  explicit HICReadPathAligner() {}
  explicit HICReadPathAligner(HICReadSet<>* hic_read_set, int binsize): hic_read_set_(hic_read_set), binsize_(binsize) {
    left_aligner_ = SingleShortReadPathAligner(&(hic_read_set_->reads_1_));
    right_aligner_ = SingleShortReadPathAligner(&(hic_read_set_->reads_2_));
  }

  //vector<HICReadAlignment> GetAlignmentsForPath(const Path& p);
  // part = 0 or 1 (0 for left, 1 for right)
  vector<SingleReadAlignment> GetPartAlignmentsForPath(const Path& p, int part);

  double eval_lambda(const Path& p);

  void RemovePathFromCache(const Path& p) {
    left_aligner_.RemovePathFromCache(p);
    right_aligner_.RemovePathFromCache(p);

    lambda_cache_.erase(p);
  }
  void RemovePathsFromCache(const vector<Path>& paths) {
    left_aligner_.RemovePathsFromCache(paths);
    right_aligner_.RemovePathsFromCache(paths);
    for (auto &p: paths) {
      lambda_cache_.erase(p);
    }
  }

  HICReadSet<>* hic_read_set_;
  int binsize_;

  SingleShortReadPathAligner left_aligner_, right_aligner_;
 private:
  unordered_map<Path, double> lambda_cache_;
};

#endif
