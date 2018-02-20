#ifndef PATH_ALIGNER_H__
#define PATH_ALIGNER_H__

#include <algorithm>
#include "read_set.h"
#include "hash_util.h"
#include "path.h"

// Handles caching of alignments
// @TODO implement caching (i8)

class SingleShortReadPathAligner {
 public:
  explicit SingleShortReadPathAligner() {}
  explicit SingleShortReadPathAligner(SingleShortReadSet<>* single_short_read_set) : single_short_read_set_(single_short_read_set) {}
  //explicit SingleShortReadPathAligner(SingleShortReadSet<StandardReadIndex>* single_short_read_set) : single_short_read_set_(single_short_read_set) {}


  vector<SingleReadAlignment> GetAlignmentsForPath(const Path& p);

  SingleShortReadSet<>* single_short_read_set_;

 private:
  const static int MAX_CACHE_SIZE = 0;
  long long next_usage_t_ = 0;
  vector<tuple<Path, vector<SingleReadAlignment>, long long>> cache_;
  void InsertAlignmentForPath(const Path& p, vector<SingleReadAlignment>& al) {
    int pos = GetAlignmentPos(p);
    if (pos == -1) {
      sort(al.begin(), al.end());
      cache_.push_back(make_tuple(p, al, next_usage_t_++));
      while (cache_.size() > MAX_CACHE_SIZE) {
        int oldest_pos = 0;
        for (int i = 0; i < (int)cache_.size(); i++) {
          const auto &t = cache_[i];
          const auto &o = cache_[oldest_pos];
          if (get<2>(t) < get<2>(o)) {
            oldest_pos = i;
          }
        }
        swap(cache_[oldest_pos], cache_.back());
        cache_.pop_back();
      }
    }
  }

  int GetAlignmentPos(const Path& p) const{
    int res = -1;
    for (int i = 0; i < (int)cache_.size(); i++) {
      const auto t = cache_[i];
      if (get<0>(t).IsSame(p)) {
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
    // @TODO initialize single aligners
    left_aligner_ = SingleShortReadPathAligner(&(paired_read_set_->reads_1_));
    right_aligner_ = SingleShortReadPathAligner(&(paired_read_set_->reads_2_));
  }

  vector<PairedReadAlignment> GetAlignmentsForPath(const Path& p);
  // part = 0 or 1 (0 for left, 1 for right)
  vector<SingleReadAlignment> GetPartAlignmentsForPath(const Path& p, int part);

  ShortPairedReadSet<>* paired_read_set_;

  SingleShortReadPathAligner left_aligner_, right_aligner_;

 private:
  const static int MAX_CACHE_SIZE = 0;
  long long next_usage_t_ = 0;
  vector<tuple<Path, vector<PairedReadAlignment>, long long>> cache_;
  void InsertAlignmentForPath(const Path& p, vector<PairedReadAlignment>& al) {
    int pos = GetAlignmentPos(p);
    if (pos == -1) {
      sort(al.begin(), al.end());
      cache_.push_back(make_tuple(p, al, next_usage_t_++));
      while (cache_.size() > MAX_CACHE_SIZE) {
        int oldest_pos = 0;
        for (int i = 0; i < (int)cache_.size(); i++) {
          const auto &t = cache_[i];
          const auto &o = cache_[oldest_pos];
          if (get<2>(t) < get<2>(o)) {
            oldest_pos = i;
          }
        }
        swap(cache_[oldest_pos], cache_.back());
        cache_.pop_back();
      }
    }
  }

  int GetAlignmentPos(const Path& p) const{
    int res = -1;
    for (int i = 0; i < (int)cache_.size(); i++) {
      const auto t = cache_[i];
      if (get<0>(t).IsSame(p)) {
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

};

#endif
