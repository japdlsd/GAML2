#ifndef MOVES_H__
#define MOVES_H__

#include "path.h"
#include "read_probability_calculator.h"

class MoveConfig {
 public:
  int big_node_threshold;
  int rand_extend_step_threshold;
  int rand_extend_distance_threshold;
  MoveConfig() :
    big_node_threshold(500),
    rand_extend_step_threshold(50),
    rand_extend_distance_threshold(1000)
    {}
};

string MakeMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config,
                GlobalProbabilityCalculator& probability_calculator, bool& accept_higher_prob);
pair<bool, string> TryMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config,
                           GlobalProbabilityCalculator& probability_calculator, bool& accept_higher_prob, vector<int>& ratios, int& move_type);

bool BreakPaths(const vector<Path>& paths, vector<Path>& out_paths,
                const MoveConfig& config);

bool UntangleCrossedPaths(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator, bool aggressive);

#endif
