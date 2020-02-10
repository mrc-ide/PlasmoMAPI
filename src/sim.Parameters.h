
#pragma once

#include <Rcpp.h>
#include <tuple>

#include "Sampler_v3.h"

//------------------------------------------------
// parameters of individual-based simulation model
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // genetic parameters
  int L;
  double prob_cotransmission;
  
  // epidemiological parameters
  double a;
  double mu;
  int u;
  int v;
  int g;
  std::vector<double> prob_infection;
  int n_prob_infection;
  std::vector<double> duration_infection;
  int n_duration_infection;
  double infectivity;
  int max_innoculations;
  
  // deme parameters
  int H;
  std::vector<int> seed_infections;
  std::vector<int> M_vec;
  int n_demes;
  
  // migration
  std::vector<std::tuple<int, int, double>> mig_list;
  int n_mig_list;
  
  // demography
  std::vector<double> life_table;
  std::vector<double> age_death;
  std::vector<double> age_stable;
  int n_age;
  
  // run parameters
  std::vector<int> time_out;
  int n_time_out;
  int max_time;
  bool report_progress;
  
  // misc parameters
  double prob_v_death;  // daily probability of mosquito death
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
  // methods
  void print_summary();
};
