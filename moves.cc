#include "moves.h"
//#include "read_probability_calculator.h"
#include <algorithm>
#include <cassert>
#include "graph.h"
#include <utility>
#include <vector>
#include "path.h"
#include "node.h"
#include <unordered_set>
#include <unordered_map>

int FindPathWithSameEnding(vector<Path>& paths, int pi) {
  for (size_t i = 0; i < paths.size(); i++) {
    auto &p = paths[i];
    auto rp = p.GetReverse();
    if ((int)i == pi) continue;
    if (p[0] == paths[pi].back()) {
      return i;
    }
    if (rp[0] == paths[pi].back()) {
      paths[i] = rp;
      return i;
    }
  }
  return -1;
}

void MergePaths(vector<Path>& paths, int pi, int same_end) {
  assert(paths[pi].back() == paths[same_end][0]);
  paths[pi].AppendPath(paths[same_end], 1);
  paths[pi].history_ += PATH_EXTEND_RANDOMLY;
  swap(paths[same_end], paths[paths.size()-1]);
  paths.pop_back();
}

bool ExtendPathsRandomly(const vector<Path>& paths, vector<Path>& out_paths,
                         const MoveConfig& config) {
  int pi = rand()%paths.size();
  out_paths = paths;
  random_shuffle(out_paths.begin(), out_paths.end());
  if (rand()%2 == 1) {
    out_paths[pi].Reverse();
  }
  if (!out_paths[pi].ExtendRandomly(config.big_node_threshold,
                                    config.rand_extend_step_threshold,
                                    config.rand_extend_distance_threshold)) {
    return false;
  }
  int same_end = FindPathWithSameEnding(out_paths, pi);
  if (same_end == -1) {
    //return true;
    // @TODO add to config
    return (rand() % 10 == 0);
  }

  MergePaths(out_paths, pi, same_end);
  return true;
}

bool BreakPaths(const vector<Path>& paths, vector<Path>& out_paths,
                const MoveConfig& config) {
  int pi = rand()%paths.size();
  out_paths = paths;
  if (out_paths[pi].size() < 2) {
    return false;
  }
  int break_pos = 1+rand()%(out_paths[pi].size()-1);
  Path p2 = out_paths[pi].CutAt(break_pos, config.big_node_threshold);
  out_paths[pi].history_ = out_paths[pi].history_ + PATH_BREAK;
  p2.history_ = out_paths[pi].history_;
  out_paths.push_back(p2);
  return true;
}

int chooseWeightedRandomly(const vector<int> &weights) {
  int total_weight = 0;
  for (int w: weights) total_weight += w;
  if (total_weight == 0) return -1;

  int pref_weight = 0;
  int res = -1;

  const int rw = rand() % total_weight;

  while (pref_weight <= rw) {
    res++;
    pref_weight += weights[res];
  }

  return res;
}

Path SampleRandomWalk(Node* start, const unordered_set<int>& allowed_nodes, const unordered_set<int>& desired_targets, const int MAX_LENGTH) {
  vector<Node*> nodes;
  nodes.push_back(start);
  int length = 0;

  vector<Node*> available;
  while (length < MAX_LENGTH && desired_targets.count(nodes.back()->id_) == 0) {
    Node* x = nodes.back();
    available.clear();
    for (auto nx: x->next_) {
      //cerr << nx->id_ << " ";
      if (allowed_nodes.count(nx->id_)) {
        available.push_back(nx);
      }
    }
    if (available.empty()) break;
    const int choice = rand() % available.size();
    //cerr << "choice: " << choice << endl;
    Node* chosen = available[choice];
    nodes.push_back(chosen);
    length += chosen->str_.size();
  }
  auto ret = Path(nodes);
  //cerr << "(UN?)POSSIBLE PATH: " << ret.ToDebugString() << endl;
  return ret;
}

bool JoinWithAdvicePaired(const vector<Path>& paths, vector<Path>& out_paths,
                    const MoveConfig& config, PairedReadProbabilityCalculator& pc, GlobalProbabilityCalculator& pc_global) {
  cerr << "JoinWithAdvicePaired()" << endl;

  if (paths.empty()) return false;
  // choose path to extend (by joining another path) randomly
  // @TODO maybe take into consideration amount of mismatched paired reads
  const int start_path_id = rand() % (int)paths.size();
  const Path& start_path = paths[start_path_id];

  cerr << "start_path_id: " << start_path_id << endl;
  cerr << "start_path: " << start_path.ToDebugString() << endl;

  // exclude nondisjoint paths
  vector<int> disjoint_path_ids;
  // const auto disjoint_length = (int)(pc.mean_distance_ + pc.std_distance_ * 3);
  const auto disjoint_length = (int)(pc.mean_distance_ + pc.std_distance_ * 3);
  for (int i = 0; i < (int)paths.size(); i++) {
    //cerr << "checking path[" << i << "] for disjointness with given path" << endl;
    //cerr << paths[i].ToDebugString() << endl;
    if (start_path.isPartlyDisjoint(paths[i], disjoint_length)) {
      disjoint_path_ids.push_back(i);
    }
  };
  if (disjoint_path_ids.empty()) return false;

  cerr << "Disjoint path ids:";
  for (int i: disjoint_path_ids) cerr << " " << i;
  cerr << endl;

  // align paired reads to the start path
  vector<SingleReadAlignment> start_left_alignments = pc.path_aligner_.GetPartAlignmentsForPath(start_path, 0);
  //sort(start_left_alignments.begin(), start_left_alignments.end());
  vector<SingleReadAlignment> start_right_alignments = pc.path_aligner_.GetPartAlignmentsForPath(start_path, 1);
  //sort(start_right_alignments.begin(), start_right_alignments.end());

  //cerr << "paired reads aligned to the start path" << endl;


  unordered_map<int,int> start_left_count, start_right_count;
  for (auto al: start_left_alignments){
    start_left_count[al.read_id] = 1 + start_left_count[al.read_id];
  }
  for (auto al: start_right_alignments) {
    start_right_count[al.read_id] = 1 + start_right_count[al.read_id];
  }

  vector<int> path_score(paths.size(), 0);
  vector<SingleReadAlignment> als1, als2;
  for (int i: disjoint_path_ids) {
    const Path& p = paths[i];
    als1 = pc.path_aligner_.GetPartAlignmentsForPath(p, 0);
    als2 = pc.path_aligner_.GetPartAlignmentsForPath(p, 1);

    for (auto al1: als1) {
      auto it = start_right_count.find(al1.read_id);
      if (it != start_right_count.end()) {
        path_score[i] += it->second;
      }
    }
    for (auto al2: als2) {
      auto it = start_left_count.find(al2.read_id);
      if (it != start_left_count.end()) {
        path_score[i] += it->second;
      }
    }
  }
  for (auto &x: path_score) {
    if (x < pc.lower_bound_advice_count_) x = 0;
  }

  cerr << "DISJOINT PATHS WITH SCORES:" << endl;
  for (int path_id: disjoint_path_ids) {
    if (path_score[path_id] > 0)  cerr << path_score[path_id] << "\t!! " << paths[path_id].ToDebugString() << endl;
    else cerr << "-\t   " << paths[path_id].ToDebugString() << endl;
  }

  cerr << "PATH SCORES:" << endl;
  for (auto x: path_score) cerr << x << "\t"; cerr << endl;

  // choose randomly (based on score) the walk to join (with the orientation and complementarity)
  const int target_path_id = chooseWeightedRandomly(path_score);
  cerr << "Chosen target path: " << target_path_id;
  if (target_path_id == -1) return false;
  const Path& target_path = paths[target_path_id];
  cerr << " :: " << target_path.ToDebugString() << endl;

  // sample (randomly) possible connections
  // choose the best

  // start path:  xb ----> xe. xbr := xb->rc_, xer := xe->rc_
  // target path: yb ----> ye. ybr := yb->rc_, yer := ye->rc_
  Node* xb = start_path.nodes_.front();
  Node* xe = start_path.nodes_.back();
  //Node* xbr = xb->rc_;
  Node* xer = xe->rc_;
  Node* yb = target_path.nodes_.front();
  Node* ye = target_path.nodes_.back();
  //Node* ybr = yb->rc_;
  Node* yer = ye->rc_;

  // 4 possibilities: xe ---> (yb | yer) ;  ye ---> ( xb | xer)
  Graph* G = xb->graph_;
  const vector<Node*> yb_drb = G->GetDrainageBasinForNode(yb);
  const vector<Node*> yer_drb = G->GetDrainageBasinForNode(yer);

  const vector<Node*> xb_drb = G->GetDrainageBasinForNode(xb);
  const vector<Node*> xer_drb = G->GetDrainageBasinForNode(xer);

  //cerr << "yb_drb\t:: "; for(auto x: yb_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "yer_drb\t:: "; for(auto x: yer_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "xb_drb\t:: "; for(auto x: xb_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "xer_drb\t:: "; for(auto x: xer_drb) cerr << x->id_ << " "; cerr << endl;

  unordered_set<int> first_pool_ids, second_pool_ids;

  for (const auto &A: {yb_drb, yer_drb}) {
    if (find(A.begin(), A.end(), xe) != A.end()) {
      for (auto n: A) {
        first_pool_ids.insert(n->id_);
      }
    }
  }

  for (const auto &A: {xb_drb, xer_drb}) {
    if (find(A.begin(), A.end(), ye) != A.end()) {
      for (auto n: A) {
        second_pool_ids.insert(n->id_);
      }
    }
  }

  //cerr << "first pool ids: [size: " << first_pool_ids.size() << "] "; for(int x: first_pool_ids) cerr << x << " "; cerr << endl;
  //cerr << "second pool ids: [size: " << second_pool_ids.size() << "] "; for(int x: second_pool_ids) cerr << x << " "; cerr << endl;

  vector<Path> possible_connections;

  // @TODO add to the config
  const int SAMPLE_NUM = 300;
  const int MAX_CONN_LENGTH = disjoint_length * 3; // bases, not nodes

  if (!first_pool_ids.empty()) {
    unordered_set<int> desired_targets({yb->id_, yer->id_});
    for (int num = 0; num < SAMPLE_NUM; num++) {
      Path p = SampleRandomWalk(xe, first_pool_ids, desired_targets, MAX_CONN_LENGTH);
      if (p.back() == yb || p.back() == yer) {
        bool is_new_path = true;
        for (auto ex_p: possible_connections) {
          if (ex_p.IsSameNoReverse(p)) {
            is_new_path = false;
            break;
          }
        }
        if (is_new_path) possible_connections.push_back(p);
      }
    }
  }
  if (!second_pool_ids.empty()) {
    unordered_set<int> desired_targets({xb->id_, xer->id_});
    for (int num = 0; num < SAMPLE_NUM; num++) {
      Path p = SampleRandomWalk(ye, second_pool_ids, desired_targets, MAX_CONN_LENGTH);
      cout << "\r" << num+1 << " paths sampled";
      if (p.back() == xb || p.back() == xer) {
        bool is_new_path = true;
        for (const auto ex_p: possible_connections) {
          if (ex_p.IsSameNoReverse(p)) {
            is_new_path = false;
            break;
          }
        }
        if (is_new_path) possible_connections.push_back(p);
      }
    }
  }
  cout << endl;

  //cerr << "possible connections: " << endl;
  //for (auto &p: possible_connections) {
  //  cerr << p.ToDebugString() << endl;
  //}

  if (possible_connections.empty()) return false;

  Path best_path;
  double best_score = -10000000000;

  out_paths.clear();
  out_paths.reserve(paths.size() - 1);
  for (int i = 0; i < (int)paths.size(); i++) {
    if (i != start_path_id && i != target_path_id) {
      out_paths.push_back(paths[i]);
    }
  }


  for (auto &p: possible_connections) {
    vector<Node*> cut_nodes(p.nodes_.begin() + 1, p.nodes_.end() - 1);
    Path inter_path(cut_nodes);

    Path np;

    if (p[0] == xe) {
      if (p.nodes_.back() == yb) {
        np = Path(start_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(target_path);
      }
      if (p.nodes_.back() == yer) {
        np = Path(start_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(target_path.GetReverse());
      }
    }
    if (p[0] == ye) {
      if (p.nodes_.back() == xb) {
        np = Path(target_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(start_path);
      }
      if (p.nodes_.back() == xer) {
        np = Path(target_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(start_path.GetReverse());
      }
    }
    np.history_ = start_path.history_ + PATH_JOIN;

    out_paths.push_back(np);
    ProbabilityChanges pp_change;
    double score = pc_global.GetPathsProbability(out_paths, pp_change);
    if (score > best_score) {
      pc_global.RemovePathFromCache(best_path);
      best_path = np;
    }
    else {
      pc_global.RemovePathFromCache(np);
    }
    out_paths.pop_back();
  }
  //out_paths.pop_back(); <-- it removed some valid path not meant to be deleted
  out_paths.push_back(best_path);
  return true;
}

bool JoinWithAdviceHic(const vector<Path>& paths, vector<Path>& out_paths,
                       const MoveConfig& config, HICReadProbabilityCalculator& pc, GlobalProbabilityCalculator& pc_global) {
  cerr << "JoinWithAdviceHIC()" << endl;

  if (paths.empty()) return false;
  // choose path to extend (by joining another path) randomly
  // @TODO maybe take into consideration amount of mismatched paired reads
  const int start_path_id = rand() % (int)paths.size();
  const Path& start_path = paths[start_path_id];

  //cerr << "start_path_id: " << start_path_id << endl;
  //cerr << "start_path: " << start_path.ToDebugString() << endl;

  // exclude nondisjoint paths
  vector<int> disjoint_path_ids;
  // const auto disjoint_length = (int)(pc.mean_distance_ + pc.std_distance_ * 3);
  const auto disjoint_length = 300;
  for (int i = 0; i < (int)paths.size(); i++) {
    //cerr << "checking path[" << i << "] for disjointness with given path" << endl;
    //cerr << paths[i].ToDebugString() << endl;
    if (start_path.isPartlyDisjoint(paths[i], disjoint_length)) {
      disjoint_path_ids.push_back(i);
    }
  };
  if (disjoint_path_ids.empty()) return false;

  //cerr << "Disjoint path ids:";
  //for (int i: disjoint_path_ids) cerr << " " << i;
  //cerr << endl;

  // align paired reads to the start path
  vector<SingleReadAlignment> start_left_alignments = pc.path_aligner_.GetPartAlignmentsForPath(start_path, 0);
  //sort(start_left_alignments.begin(), start_left_alignments.end());
  vector<SingleReadAlignment> start_right_alignments = pc.path_aligner_.GetPartAlignmentsForPath(start_path, 1);
  //sort(start_right_alignments.begin(), start_right_alignments.end());

  //cerr << "paired reads aligned to the start path" << endl;


  unordered_map<int,int> start_left_count, start_right_count;
  for (auto al: start_left_alignments){
    start_left_count[al.read_id] = 1 + start_left_count[al.read_id];
  }
  for (auto al: start_right_alignments) {
    start_right_count[al.read_id] = 1 + start_right_count[al.read_id];
  }

  vector<int> path_score(paths.size(), 0);
  vector<SingleReadAlignment> als1, als2;
  for (int i: disjoint_path_ids) {
    const Path& p = paths[i];
    als1 = pc.path_aligner_.GetPartAlignmentsForPath(p, 0);
    als2 = pc.path_aligner_.GetPartAlignmentsForPath(p, 1);

    for (auto al1: als1) {
      auto it = start_right_count.find(al1.read_id);
      if (it != start_right_count.end()) {
        path_score[i] += it->second;
      }
    }
    for (auto al2: als2) {
      auto it = start_left_count.find(al2.read_id);
      if (it != start_left_count.end()) {
        path_score[i] += it->second;
      }
    }
  }
  for (auto &x: path_score) {
    if (x < pc.lower_bound_advice_count_) x = 0;
  }


  //cerr << "DISJOINT PATHS WITH SCORES:" << endl;
  //for (int path_id: disjoint_path_ids) {
  //  if (path_score[path_id] > 0)  cerr << path_score[path_id] << "\t!! " << paths[path_id].ToDebugString() << endl;
  //  else cerr << "-\t   " << paths[path_id].ToDebugString() << endl;
  //}

  // choose randomly (based on score) the walk to join (with the orientation and complementarity)
  const int target_path_id = chooseWeightedRandomly(path_score);
  if (target_path_id == -1) return false;
  const Path& target_path = paths[target_path_id];
  //cerr << "Chosen target path: " << target_path_id << " :: " << target_path.ToDebugString() << endl;

  // sample (randomly) possible connections
  // choose the best

  // start path:  xb ----> xe. xbr := xb->rc_, xer := xe->rc_
  // target path: yb ----> ye. ybr := yb->rc_, yer := ye->rc_
  Node* xb = start_path.nodes_.front();
  Node* xe = start_path.nodes_.back();
  //Node* xbr = xb->rc_;
  Node* xer = xe->rc_;
  Node* yb = target_path.nodes_.front();
  Node* ye = target_path.nodes_.back();
  //Node* ybr = yb->rc_;
  Node* yer = ye->rc_;

  // 4 possibilities: xe ---> (yb | yer) ;  ye ---> ( xb | xer)
  Graph* G = xb->graph_;
  const vector<Node*> yb_drb = G->GetDrainageBasinForNode(yb);
  const vector<Node*> yer_drb = G->GetDrainageBasinForNode(yer);

  const vector<Node*> xb_drb = G->GetDrainageBasinForNode(xb);
  const vector<Node*> xer_drb = G->GetDrainageBasinForNode(xer);

  //cerr << "yb_drb\t:: "; for(auto x: yb_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "yer_drb\t:: "; for(auto x: yer_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "xb_drb\t:: "; for(auto x: xb_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "xer_drb\t:: "; for(auto x: xer_drb) cerr << x->id_ << " "; cerr << endl;

  unordered_set<int> first_pool_ids, second_pool_ids;

  for (const auto &A: {yb_drb, yer_drb}) {
    if (find(A.begin(), A.end(), xe) != A.end()) {
      for (auto n: A) {
        first_pool_ids.insert(n->id_);
      }
    }
  }

  for (const auto &A: {xb_drb, xer_drb}) {
    if (find(A.begin(), A.end(), ye) != A.end()) {
      for (auto n: A) {
        second_pool_ids.insert(n->id_);
      }
    }
  }

  //cerr << "first pool ids: [size: " << first_pool_ids.size() << "] "; for(int x: first_pool_ids) cerr << x << " "; cerr << endl;
  //cerr << "second pool ids: [size: " << second_pool_ids.size() << "] "; for(int x: second_pool_ids) cerr << x << " "; cerr << endl;

  vector<Path> possible_connections;

  // @TODO add to the config
  const int SAMPLE_NUM = 300;
  const int MAX_CONN_LENGTH = disjoint_length * 3; // bases, not nodes

  if (!first_pool_ids.empty()) {
    unordered_set<int> desired_targets({yb->id_, yer->id_});
    for (int num = 0; num < SAMPLE_NUM; num++) {
      Path p = SampleRandomWalk(xe, first_pool_ids, desired_targets, MAX_CONN_LENGTH);
      if (p.back() == yb || p.back() == yer) {
        bool is_new_path = true;
        for (auto ex_p: possible_connections) {
          if (ex_p.IsSameNoReverse(p)) {
            is_new_path = false;
            break;
          }
        }
        if (is_new_path) possible_connections.push_back(p);
      }
    }
  }
  if (!second_pool_ids.empty()) {
    unordered_set<int> desired_targets({xb->id_, xer->id_});
    for (int num = 0; num < SAMPLE_NUM; num++) {
      Path p = SampleRandomWalk(ye, second_pool_ids, desired_targets, MAX_CONN_LENGTH);
      cout << "\r" << num+1 << " paths sampled";
      if (p.back() == xb || p.back() == xer) {
        bool is_new_path = true;
        for (const auto ex_p: possible_connections) {
          if (ex_p.IsSameNoReverse(p)) {
            is_new_path = false;
            break;
          }
        }
        if (is_new_path) possible_connections.push_back(p);
      }
    }
  }
  cout << endl;

  //cerr << "possible connections: " << endl;
  //for (auto &p: possible_connections) {
  //  cerr << p.ToDebugString() << endl;
  //}

  if (possible_connections.empty()) return false;

  Path best_path;
  double best_score = -10000000000;

  out_paths.clear();
  out_paths.reserve(paths.size() - 1);
  for (int i = 0; i < (int)paths.size(); i++) {
    if (i != start_path_id && i != target_path_id) {
      out_paths.push_back(paths[i]);
    }
  }


  for (auto &p: possible_connections) {
    vector<Node*> cut_nodes(p.nodes_.begin() + 1, p.nodes_.end() - 1);
    Path inter_path(cut_nodes);

    Path np;

    if (p[0] == xe) {
      if (p.nodes_.back() == yb) {
        np = Path(start_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(target_path);
      }
      if (p.nodes_.back() == yer) {
        np = Path(start_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(target_path.GetReverse());
      }
    }
    if (p[0] == ye) {
      if (p.nodes_.back() == xb) {
        np = Path(target_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(start_path);
      }
      if (p.nodes_.back() == xer) {
        np = Path(target_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(start_path.GetReverse());
      }
    }
    np.history_ = start_path.history_ + PATH_JOIN;

    out_paths.push_back(np);
    ProbabilityChanges pp_change;
    double score = pc_global.GetPathsProbability(out_paths, pp_change);
    if (score > best_score) {
      pc_global.RemovePathFromCache(best_path);
      best_path = np;
    }
    else {
      pc_global.RemovePathFromCache(np);
    }
    out_paths.pop_back();
  }
  //out_paths.pop_back(); <-- it removed some valid path not meant to be deleted
  out_paths.push_back(best_path);
  return true;
}

bool JoinWithAdvice(const vector<Path>& paths, vector<Path>& out_paths,
                    const MoveConfig& config, GlobalProbabilityCalculator& pc) {
  // choose dataset for advices. if no available, then return false
  vector<int> available_paired_read_sets;
  for (int i = 0; i < (int)pc.paired_read_calculators_.size(); i++) {
    auto &p = pc.paired_read_calculators_[i];
    if (p.first.use_as_advice_) {
      available_paired_read_sets.push_back(i);
    }
  }

  for (int i = 0; i < (int)pc.hic_read_calculators_.size(); i++) {
    auto &p = pc.hic_read_calculators_[i];
    if (p.first.use_as_advice_) {
      available_paired_read_sets.push_back(i + (int)pc.paired_read_calculators_.size());
    }
  }

  if (available_paired_read_sets.empty()) {
    return false;
  }

  const int which_dataset_to_choose = rand() % (int)available_paired_read_sets.size();
  if (which_dataset_to_choose < pc.paired_read_calculators_.size()) {
    return JoinWithAdvicePaired(paths,
                                out_paths,
                                config,
                                pc.paired_read_calculators_[available_paired_read_sets[which_dataset_to_choose]].first,
                                pc);
  }
  else {
    return JoinWithAdviceHic(paths,
                             out_paths,
                             config,
                             pc.hic_read_calculators_[available_paired_read_sets[which_dataset_to_choose-(int)pc.paired_read_calculators_.size()]].first,
                             pc);
  }
  // @TODO add HiC choice
}

int chooseRandomPositionByNode(Path& p, const int node_id) {
  static const int BLANK = -1000*1000*1000;
  int res = BLANK;
  int count_good_nodes = 0;
  bool do_reverse = false;
  for (int i = 0; i < (int)p.nodes_.size(); i++) {
    const auto &n = p.nodes_[i];
    if (n->id_ == node_id) {
      count_good_nodes++;
      if (rand()%count_good_nodes == 0) { // to guarantee that every node is picked with equal probability
        res = i;
        do_reverse = false;
      }
      break;
    }
    else if (n->rc_->id_ == node_id) {
      count_good_nodes++;
      if (rand()%count_good_nodes == 0) { // to guarantee that every node is picked with equal probability
        res = (int)p.nodes_.size() - i - 1;
        do_reverse = true;
      }

      break;
    }
  }
  assert(res != BLANK);
  if (do_reverse) p.Reverse();
  return res;
}

pair<int, int> findLargestCommonSubWalk(const Path& p1, const int p1_start_pos, const Path& p2, const int p2_start_pos) {
  int b_inter = 0;
  int e_inter = 0;

  while (p1_start_pos - b_inter >= 0 && p2_start_pos - b_inter >= 0 &&
      p1.nodes_[p1_start_pos - b_inter]->id_ == p2.nodes_[p2_start_pos - b_inter]->id_) {
    b_inter++;
  }

  while (p1_start_pos + e_inter < p1.size() && p2_start_pos + e_inter < p2.size() &&
      p1.nodes_[p1_start_pos + e_inter]->id_ == p2.nodes_[p2_start_pos + e_inter]->id_) {
    e_inter++;
  }
  return make_pair(b_inter, e_inter);
}

vector<Node*> truncateSmallNodes(const vector<Node*>& nodes, const int big_node_threshold) {
  vector<Node*> res;
  int start = 0;
  while (start < (int)nodes.size() && nodes[start]->str_.size() < big_node_threshold) start++;
  int end = (int)nodes.size();
  while (end-1 >= start && nodes[end-1]->str_.size() < big_node_threshold) end--;
  if (start < end) {
    res.insert(res.end(), nodes.begin() + start, nodes.begin() + end);
  }
  return res;
}

bool UntangleCrossedPaths(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator, bool aggressive=false) {
  //cerr << "UntangleCrossedPaths()" << endl;
  if (paths.size() < 2) return false;

  // find nodes that are in more than one path

  // node_id, [path_id, pos]
  static unordered_map<int, vector<pair<int, int>>> paths_by_nodes;
  paths_by_nodes.clear();

  for (int i = 0; i < (int) paths.size(); i++) {
    const auto &p = paths[i];
    for (int j = 0; j < (int) p.nodes_.size(); j++) {
      auto n = p.nodes_[j];
      if (n->str_.size() < 50) continue; // we look only at decently big nodes
      if (n->id_ > n->rc_->id_) n = n->rc_;
      const int id = n->id_;
      if (paths_by_nodes.count(id) == 0) {
        paths_by_nodes[id] = vector<pair<int, int>>({make_pair(i, j)});
      } else {
        paths_by_nodes[id].emplace_back(i, j);
      }
    }
  }

  // debug
  if (0) {
    cerr << "PATHS BY NODES: " << endl;
    for (auto n: paths_by_nodes) {
      cerr << n.first << ": ";
      for (auto p: n.second) {
        cerr << "(" << p.first << ", " << p.second << ") ";
      }
      cerr << endl;
    }
  }

  vector<int> candidate_nodes;
  vector<unordered_set<int>> candidate_paths;
  unordered_set<int> path_buf;
  for (const auto x: paths_by_nodes) {
    const int nid = x.first;
    path_buf.clear();
    for (const auto p: x.second) {
      path_buf.insert(p.first);
    }
    if (path_buf.size() >= 2) {
      candidate_nodes.push_back(nid);
      candidate_paths.push_back(path_buf);
    }
  }

  if (candidate_nodes.empty()) return false;
  // debug
  if (0) {
    cerr << "CANDIDATE NODES: " << endl;
    for (auto x: candidate_nodes) cerr << x << " ";
    cerr << endl;
  }

  const int r_pos = rand() % (int) candidate_nodes.size();
  const int inter_node_id = candidate_nodes[r_pos];
  vector<int> inter_paths(candidate_paths[r_pos].begin(), candidate_paths[r_pos].end());
  const int r1_pos = rand() % (int) inter_paths.size();
  swap(inter_paths[0], inter_paths[r1_pos]);
  const int r2_pos = 1 + rand() % ((int) (inter_paths).size() - 1);
  swap(inter_paths[1], inter_paths[r2_pos]);

  const int p1_id = inter_paths[0];
  const int p2_id = inter_paths[1];

  Path p1 = paths[p1_id];
  Path p2 = paths[p2_id];

  int p1_start_pos = chooseRandomPositionByNode(p1, inter_node_id);
  int p2_start_pos = chooseRandomPositionByNode(p2, inter_node_id);

  int p1_label = 0;
  int p2_label = 1;

  cerr << "INTER NODE ID: " << inter_node_id << endl;
  cerr << "PATH1: " << p1.ToDebugString() << endl;
  cerr << "PATH2: " << p2.ToDebugString() << endl;
  cerr << "p1_start_pos: " << p1_start_pos << ", p2_start_pos: " << p2_start_pos << endl;

  pair<int, int> be = findLargestCommonSubWalk(p1, p1_start_pos, p2, p2_start_pos);
  cerr << "COMMON SUBWALK: " << be.first << " " << be.second << endl;

  // create new hypothesis
  vector<vector<Path>> added_paths;

  for (int i = 0; i < 2; i++) {
    // cut the second, keep the first
    swap(p1, p2);
    swap(p1_start_pos, p2_start_pos);
    swap(p1_label, p2_label);
    vector<Path> add;

    vector<Node *> a(p2.nodes_.begin(), p2.nodes_.begin() + (p2_start_pos - be.first + 1));
    vector<Node *> b(p2.nodes_.begin() + (p2_start_pos + be.second), p2.nodes_.end());
    a = truncateSmallNodes(a, config.big_node_threshold);
    b = truncateSmallNodes(b, config.big_node_threshold);

    if (!a.empty()) add.emplace_back(a, p2.history_ + PATH_UNTANGLE);
    if (!b.empty()) add.emplace_back(b, p2.history_ + PATH_UNTANGLE);
    add.push_back(p1);

    added_paths.push_back(add);

    // aggressive mode: remove cut paths
    if (aggressive) {
      add.clear();
      add.push_back(p1);
      added_paths.push_back(add);
    }
  }

  for (int i = 0; i < 2; i++) {
    swap(p1, p2);
    swap(p1_start_pos, p2_start_pos);
    swap(p1_label, p2_label);
    // start with the first, finish with second, cut the rest
    vector<Path> add;
    //vector<int> remove;

    vector<Node *> a(p1.nodes_.begin(), p1.nodes_.begin() + (p1_start_pos + be.second));
    a.insert(a.end(), p2.nodes_.begin() + (p2_start_pos + be.second), p2.nodes_.end());
    a = truncateSmallNodes(a, config.big_node_threshold);

    vector<Node *> b(p2.nodes_.begin(), p2.nodes_.begin() + (p2_start_pos - be.first + 1));
    b = truncateSmallNodes(b, config.big_node_threshold);
    vector<Node *> c(p1.nodes_.begin() + (p1_start_pos + be.second), p1.nodes_.end());
    c = truncateSmallNodes(c, config.big_node_threshold);

    if (!a.empty()) add.emplace_back(a, p1.history_ + PATH_UNTANGLE);
    if (!b.empty()) add.emplace_back(b, p2.history_ + PATH_UNTANGLE);
    if (!c.empty()) add.emplace_back(c, p2.history_ + PATH_UNTANGLE);

    added_paths.push_back(add);
    // aggressive mode: remove the rest
    if (aggressive) {
      add.clear();
      if (!a.empty()) add.emplace_back(a, p1.history_ + PATH_UNTANGLE);
    }
  }
  {
    // switch first end with second and vice versa
    vector<Path> add;
    //vector<int> remove;

    vector<Node *> a(p1.nodes_.begin(), p1.nodes_.begin() + (p1_start_pos + be.second));
    a.insert(a.end(), p2.nodes_.begin() + (p2_start_pos + be.second), p2.nodes_.end());
    a = truncateSmallNodes(a, config.big_node_threshold);

    vector<Node *> b(p2.nodes_.begin(), p2.nodes_.begin() + (p2_start_pos + be.second));
    b.insert(b.end(), p1.nodes_.begin() + (p1_start_pos + be.second), p1.nodes_.end());
    b = truncateSmallNodes(b, config.big_node_threshold);

    if (!a.empty()) add.emplace_back(a, p1.history_ + PATH_UNTANGLE);
    if (!b.empty()) add.emplace_back(b, p2.history_ + PATH_UNTANGLE);

    added_paths.push_back(add);

  }

  // debug
  if (0) {
    for (int k = 0; k < (int) added_paths.size(); k++) {
      // debug
      {
        cerr << "added_paths: " << k << endl;
        for (int i = 0; i < (int) added_paths[k].size(); i++) {
          cerr << added_paths[k][i].ToDebugString() << endl;
        }
      }
    }
  }

  /*
  vector<Path> new_paths(paths);
  swap(new_paths[new_paths.size() - 1], new_paths[p1_id]);
  swap(new_paths[new_paths.size() - 2], new_paths[p2_id]);
  new_paths.pop_back();
  new_paths.pop_back();
  */
  vector<Path> new_paths;
  new_paths.reserve(paths.size() + 1);
  for (auto &p: paths) {
    if (!p.IsSame(p1) && !p.IsSame(p2)) new_paths.push_back(p);
  }
  //cerr << "old: ";
  //cerr << PathsToDebugString(paths);
  //cerr << "new without changed: ";
  //cerr << PathsToDebugString(new_paths);
  //cerr << "new_paths.size(): " << new_paths.size() << "; paths.size(): " << paths.size() << endl;
  //assert(new_paths.size() + 2 == paths.size()); <- it doesn't have to hold

  if (aggressive) {
    auto random_k = rand() % added_paths.size();
    new_paths.insert(new_paths.end(), added_paths[random_k].begin(), added_paths[random_k].end());
    out_paths = new_paths;
  }
  else {
    double best_prob = -100000000;
    int best_k = 0;

    const auto other_paths_size = new_paths.size();

    for (int k = 0; k < (int) added_paths.size(); k++) {
      new_paths.resize(other_paths_size);
      new_paths.insert(new_paths.end(), added_paths[k].begin(), added_paths[k].end());
      // debug
      if (0) {
        cerr << "added_paths: \n";
        for (int i = 0; i < (int) added_paths[k].size(); i++) {
          cerr << added_paths[k][i].ToDebugString() << endl;
        }
      }
      ProbabilityChanges pp_global;
      double prob = probability_calculator.GetPathsProbability(new_paths, pp_global);
      if (prob > best_prob) {
        best_prob = prob;
        best_k = k;
      }
    }

    for (int k = 0; k < (int) added_paths.size(); k++) {
      if (k != best_k) probability_calculator.RemovePathsFromCache(added_paths[k]);
    }

    const int removed_paths_length = (int) p1.GetLength() + (int) p2.GetLength();
    int added_paths_length = 0;
    for (auto &p: added_paths[best_k]) added_paths_length += p.GetLength();

    if (0) {
      cerr << "removed paths: \n";
      cerr << "[L:" << p1.GetLength() << "]\t" << p1.ToDebugString() << endl;
      cerr << "[L:" << p2.GetLength() << "]\t" << p2.ToDebugString() << endl;

      cerr << "best added_paths: \n";
      for (int i = 0; i < (int) added_paths[best_k].size(); i++) {
        cerr << "[L:" << added_paths[best_k][i].GetLength() << "]\t" << added_paths[best_k][i].ToDebugString() << endl;
      }
    }

    //assert(removed_paths_length >= added_paths_length - 2 * (p1[0]->graph_->k_ - 1));

    new_paths.resize(other_paths_size);
    new_paths.insert(new_paths.end(), added_paths[best_k].begin(), added_paths[best_k].end());
    out_paths = new_paths;
  }
  return true;
}

string MakeMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator,
              bool& accept_higher_prob) {
  // @TODO add probs of moves into config
  vector<int> ratios = {5, 5, 20, 20};

  int move_type = 0;
  while (true) {
    out_paths.clear();
    pair<bool, string> res = TryMove(paths, out_paths, config, probability_calculator, accept_higher_prob, ratios, move_type);
    if (res.first) {
      return res.second;
    }
    else {
      if (ratios[move_type] > 1) ratios[move_type] -= 1;
    }
  }
}

pair<bool, string> TryMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator,
             bool& accept_higher_prob, vector<int>& ratios, int& move_type) {
  move_type = chooseWeightedRandomly(ratios);

  if (move_type == 0) {
    accept_higher_prob = false;
    cout << "ExtendPathsRandomly()" << endl;
    return make_pair(ExtendPathsRandomly(paths, out_paths, config), PATH_EXTEND_RANDOMLY);
  }
  if (move_type == 1) {
    accept_higher_prob = true;
    cout << "BreakPaths()" << endl;
    return make_pair(BreakPaths(paths, out_paths, config), PATH_BREAK);
  }
  // @TODO Joining with advice move (high priority)
  if (move_type == 2) {
    cout << "JoinWithAdvice()" << endl;
    accept_higher_prob = false; // @TODO check with Usama if correct
    return make_pair(JoinWithAdvice(paths, out_paths, config, probability_calculator), PATH_JOIN);
  }
  if (move_type == 3) {
    cout << "UntangleCrossedPaths()" << endl;
    accept_higher_prob = false;
    return make_pair(UntangleCrossedPaths(paths, out_paths, config, probability_calculator), PATH_UNTANGLE);
  }

  // @TODO add Local improvement move (low priority)
  // @TODO add Repeat optimization move (low priority)
  // @TODO Repeat interchange move (low priority)
  return make_pair(false, "");
}
