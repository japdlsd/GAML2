#include "path_aligner.h"

vector<SingleReadAlignment> SingleShortReadPathAlignerVector::GetAlignmentsForPath(const Path& p) {
  return GetAlignmentForPathWithCache(p);
}
vector<SingleReadAlignment> SingleShortReadPathAlignerVector::GetAlignmentForPathNoCache(const Path &p) {
  vector<SingleReadAlignment> ret;
  string genome = p.ToString(true);
  ret = single_short_read_set_->GetAlignments(genome);
  return ret;
}
vector<SingleReadAlignment> SingleShortReadPathAlignerVector::GetAlignmentForPathWithCache(const Path &p) {
  vector<SingleReadAlignment> ret;

  const int pos = GetAlignmentPos(p);
  if (pos == -1) {
    string genome = p.ToString(true);
    ret = single_short_read_set_->GetAlignments(genome);
    InsertAlignmentForPath(p, ret);
  }
  else {
    ret = GetCachedAlignmentByPos(pos);
  }

  return ret;
}
vector<PairedReadAlignment> PairedReadPathAligner::GetAlignmentsForPath(const Path &p) {
  auto it = cache_.find(p);
  if (it != cache_.end()) {
    //cout << "cache hit" << endl;
    return it->second;
  }
  vector<PairedReadAlignment> res;

  auto als1 = left_aligner_.GetAlignmentsForPath(p);
  //cout << "left als: " << als1.size() << endl;
  auto als2 = right_aligner_.GetAlignmentsForPath(p);
  //cout << " right als: " << als2.size() << " " << endl;
  // assume they are sorted

  if (!als1.empty() && !als2.empty()) {
    auto it1 = als1.begin();
    auto it2 = als2.begin();

    while (it1 != als1.end() && it2 != als2.end()) {
      //cout << "big cycle" << endl;
      while (it2 != als2.end() && it2->read_id < it1->read_id) it2++;
      if (it2 == als2.end()) {
        //cout << "it2 has ended" << endl;
        break;
      }
      else if (it2->read_id == it1->read_id) {
        int righties_count = 0;

        const int read_id = it1->read_id;

        while (it2 != als2.end() && it2->read_id == read_id) {
          righties_count += 1;
          it2++;
        }
        //cout << "righties size: " << righties.size() << endl;
        while (it1 != als1.end() && it1->read_id == read_id) {
          it2 -= righties_count;
          for (int i = 0; i < righties_count; i++) {
            //cout << "eval orientation" << endl;
            const pair<string, int> characteristics = eval_orientation(*it1, (int)paired_read_set_->reads_1_[read_id].size(), *it2, (int)paired_read_set_->reads_2_[read_id].size());
            res.emplace_back(*it1, *it2, characteristics.first, characteristics.second);
            it2++;
          }
          it1++;
        }
      }
      else { //  (it2->read_id > read, (i.e. we didn't hit any good reads)
        while (it1 != als1.end() && it1->read_id < it2->read_id) it1++;
      }
    }
  }
  
  cache_[p] = res;
  return res;
}


vector<SingleReadAlignment> PairedReadPathAligner::GetPartAlignmentsForPath(const Path &p, int part) {
  vector<SingleReadAlignment> ret;
  if (part == 0) {
    ret = left_aligner_.GetAlignmentsForPath(p);
  }
  else if (part == 1) {
    ret = right_aligner_.GetAlignmentsForPath(p);
  }
  return ret;
}
vector<SingleReadAlignment> HICReadPathAligner::GetPartAlignmentsForPath(const Path &p, int part) {
  vector<SingleReadAlignment> ret;
  if (part == 0) {
    ret = left_aligner_.GetAlignmentsForPath(p);
  }
  else if (part == 1) {
    ret = right_aligner_.GetAlignmentsForPath(p);
  }
  return ret;
}
double HICReadPathAligner::eval_lambda(const Path &p) {
  auto it = lambda_cache_.find(p);
  if (it != lambda_cache_.end()) return it->second;

  int total_count = 0;
  double res = 0;

  const auto als1 = GetPartAlignmentsForPath(p, 0);
  if (als1.empty()) return 0;
  const auto als2 = GetPartAlignmentsForPath(p, 1);
  if (als2.empty()) return 0;

  const auto &reads_1_ = *(left_aligner_.single_short_read_set_);
  const auto &reads_2_ = *(right_aligner_.single_short_read_set_);

  auto it1 = als1.begin();
  auto it2 = als2.begin();

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
          const pair<string,int> or_ins = eval_orientation(*it1, (int)reads_1_[read_id].size(), *it2, (int)reads_2_[read_id].size());
          const int insert_length = or_ins.second;
          if (total_count == 0) {
            total_count = 1;
            res = insert_length;
          }
          else {
            res = res * total_count + insert_length;
            total_count += 1;
            res /= total_count;
          }
          it2++;
        }
        it1++;
      }
    }
    else { //  (it2->read_id > read, (i.e. we didn't hit any good reads)
      while (it1 != als1.end() && it1->read_id < it2->read_id) it1++;
    }
  }
  lambda_cache_[p] = res;
  return res;
}
