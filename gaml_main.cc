#include "graph.h"
#include "path.h"
#include "read_set.h"
#include "read_probability_calculator.h"
#include "moves.h"
#include "config.pb.h"
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <csignal>

bool NOT_SIGINTED = true;
void signalHandler(int signum) {
  NOT_SIGINTED = false;
}

void DisconnectPoorlyConnectedPaths(GlobalProbabilityCalculator& pc, const Config& gaml_config, vector<Path>& paths, ofstream& iter_log, int start_iter) {
  MoveConfig move_conf;
  ProbabilityChanges prob_changes;
  double old_prob = pc.GetPathsProbability(paths, prob_changes);
  auto total_size = prob_changes.getLength();
  vector<Path> out_paths;
  for (int it_num = 1; it_num <= gaml_config.postprocessing_break_iter_num(); it_num++) {
    out_paths.clear();

    int move_type = rand() % 2;
    if (move_type == 0) {
      BreakPaths(paths, out_paths, move_conf);
    }
    else {
      UntangleCrossedPaths(paths, out_paths, move_conf, pc, true);
    }
    double new_prob = pc.GetPathsProbability(out_paths, prob_changes);

    if (old_prob - new_prob < gaml_config.postprocessing_break_threshold()) {
      pc.CommitProbabilityChanges(prob_changes);
      paths = out_paths;
      old_prob = new_prob;
      total_size = prob_changes.getLength();
    }

    {
      iter_log << it_num + start_iter <<"," << old_prob << "," << total_size << "," << pc.GetUncoveredBasesCount()
               <<  "," << paths.size() << ", \'" << (move_type == 0 ? PATH_BREAK : PATH_UNTANGLE) << "\'," << 0 << "," << pc.GetUnalignedReadsLog() << endl;
    }
  }
}

void PerformOptimization(GlobalProbabilityCalculator& probability_calculator,
                         const Config& gaml_config, vector<Path>& paths) {
  // @TODO add random seed to config? (i14)
  const int random_seed = 47;
  cout << "random seed: " << random_seed << endl;
  default_random_engine generator(random_seed);

  ProbabilityChanges prob_changes;
  double old_prob = probability_calculator.GetPathsProbability(paths, prob_changes);
  cout << "starting probability: " << old_prob << endl;

  probability_calculator.CommitProbabilityChanges(prob_changes);

  //cout << "INITIAL PATHS:" << endl;
  //cout << PathsToDebugString(paths) << endl;

  ofstream iter_log(gaml_config.iter_log_file());
  iter_log.precision(15);

  MoveConfig move_config;
  string move_type;
  long long total_size = prob_changes.getLength();

  const int normal_iter_num = gaml_config.num_iterations();
  const int finish_iter_num = gaml_config.finishing_iterations();
  const int total_iter_num = normal_iter_num + finish_iter_num;

  int it_num;
  for (it_num = 1; it_num <= total_iter_num && NOT_SIGINTED; it_num++) {
    double T;
    if (it_num <= normal_iter_num) {
      T = gaml_config.t0() / log(it_num / gaml_config.n_divisor() + 1);
    }
    else {
      T = 0;
    }

    cout.precision(15);
    cout << "ITERATION: " << it_num << "\t T: " << T << "\t PROB: " << old_prob << "\t SIZE: " << total_size <<  endl;
    cout << "UNALIGNED:" << probability_calculator.GetUnalignedReadsDebug() << endl;
    cout << "PROB HISTS:"  << endl << probability_calculator.GetProbHists() << endl;
    cout.precision(6);

    vector<Path> new_paths;
    bool accept_high_prob;
    move_type = MakeMove(paths, new_paths, move_config, probability_calculator, accept_high_prob);
    double new_prob = probability_calculator.GetPathsProbability(new_paths, prob_changes);
    cout << "PROPOSED PROB: " << new_prob;

    bool accept = false;
    if (new_prob > old_prob) {
      accept = true;
      cout << " better than old; ";
    } else if (accept_high_prob && T > 0) {
      double prob = exp((new_prob - old_prob) / T);
      cerr << "(prob: " << prob << ")\t";
      uniform_real_distribution<double> dist(0.0, 1.0);
      double samp = dist(generator);
      if (samp < prob) {
        cout << " worse, but temperature; ";
        accept = true;
      }
    }

    if (accept) {
      cout << "\tACCEPT";
      if (old_prob < new_prob){
        // debug
        for (auto &np: new_paths) {
          bool is_new = true;
          for (auto &p: paths) {
            if (p.IsSame(np)) {
              is_new = false;
              break;
            }
          }
          if (is_new) np.history_ += PATH_BETTER;
        }
      }
      old_prob = new_prob;

      paths = new_paths;
      probability_calculator.CommitProbabilityChanges(prob_changes);
      total_size = prob_changes.getLength();

      // debug log
      if ( !prob_changes.paired_read_changes.empty() ) {
        cout << "\nACCEPTED ADDINGS: " << "\n";
        cout << PathsToDebugString(prob_changes.paired_read_changes[0].added_paths) << "\n";
        cout << "ACCEPTED REMOVALS: " << "\n";
        cout << PathsToDebugString(prob_changes.paired_read_changes[0].removed_paths) << endl;
      }
    }
    else {
      if (!probability_calculator.paired_read_calculators_.empty()) {
        probability_calculator.RemovePathsFromCache(prob_changes.paired_read_changes[0].added_paths);
      }
      else if (!probability_calculator.hic_read_calculators_.empty()) {
        probability_calculator.RemovePathsFromCache(prob_changes.hic_read_changes[0].added_paths);
      }
      else if (!probability_calculator.single_read_calculators_.empty()) {
        probability_calculator.RemovePathsFromCache(prob_changes.single_read_changes[0].added_paths);
      }

    }

    // continual output

    if (it_num % gaml_config.output_flush_freq() == 0) {
      ofstream of(gaml_config.output_file());
      PathsToFasta(paths, of);
    }

    {
      iter_log << it_num <<"," << old_prob << "," << total_size << "," << probability_calculator.GetUncoveredBasesCount() <<  "," << paths.size() << ", \'" << move_type << "\'," << T << "," << probability_calculator.GetUnalignedReadsLog() << endl;
    }
  }

  DisconnectPoorlyConnectedPaths(probability_calculator, gaml_config, paths, iter_log, it_num);

  // final output
  ofstream of(gaml_config.output_file());
  PathsToFasta(paths, of);
  iter_log.close();
}

void global_run_logging(const string& log_filename, int argc, char** argv) {
  fstream fs;
  fs.open(log_filename, fstream::out | fstream::app);
  fs << "Running GAML2..." << endl;
  const auto ctt = time(0);
  fs << "TIMESTAMP: " << asctime(localtime(&ctt)); // contains \n symbol
  fs << "COMMAND: ";
  for (int i = 0; i < argc; i++) {
    fs << argv[i] << " ";
  }
  fs << endl;

  fstream cfg_file;
  cfg_file.open(argv[1], fstream::in);
  string buff = "x";
  while (cfg_file >> buff) {
    fs << "\t" << buff << endl;
  }

  cfg_file.close();
  fs.close();
}

int main(int argc, char** argv) {
  // Adding correct handling of Ctrl+C signal (SIGINT)
  signal(SIGINT, signalHandler);

  ifstream config_file(argv[1]);
  google::protobuf::io::IstreamInputStream config_stream(&config_file);

  Config gaml_config;
  google::protobuf::TextFormat::Parse(&config_stream, &gaml_config);

  global_run_logging(gaml_config.global_run_log_file(), argc, argv);

  cout << gaml_config.starting_graph() << endl; 

  Graph *g = LoadGraph(gaml_config.starting_graph());
  cout << "Loaded graph with " << g->nodes_.size() << " nodes" << endl;

  GlobalProbabilityCalculator probability_calculator(gaml_config);

  vector<Path> paths = BuildPathsFromSingleNodes(g->GetBigNodes(gaml_config.big_nodes_threshold()));

  cout << PathsToDebugString(paths) << endl;

  PerformOptimization(probability_calculator, gaml_config, paths);
}
