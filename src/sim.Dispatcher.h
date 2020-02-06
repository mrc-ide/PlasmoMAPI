
#pragma once

#include <set>

#include "sim.Parameters.h"
#include "sim.Host.h"
#include "sim.Mosquito.h"
#include "misc_v4.h"
#include "array.h"

//------------------------------------------------
// class defining individual-based simulation model
class Dispatcher  {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointers to inputs
  Parameters* param_ptr;
  Rcpp::Function* update_progress_ptr;
  Rcpp::List* args_progress_ptr;
  
  // make local copies of some parameters
  int n_demes;
  int max_time;
  double a;
  int v;
  double prob_v_death;
  int H;
  double infectivity;
  double mu;
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  Sampler sampler_duration_infection;
  
  // scheduler objects
  std::vector<std::set<int>> schedule_death;
  std::vector<std::vector<std::pair<int, int>>> schedule_Eh_to_Ih;
  std::vector<std::vector<std::pair<int, int>>> schedule_Ih_to_Sh;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_recovery;
  
  // counts of host types
  int H_total;
  std::vector<int> Sh;
  std::vector<int> Eh;
  std::vector<int> Ih;
  
  // population of human hosts
  std::vector<Host> host_pop;
  int next_host_ID;
  
  // store the integer index of hosts in each deme
  array_2d_int host_index;
  std::vector<std::vector<int>> host_infective_index;
  
  // counts of mosquito types
  int M_total;
  std::vector<int> Sv;
  std::vector<int> Ev;
  std::vector<int> Iv;
  
  // number of mosquitoes at various stages
  std::vector<int> n_Ev_death_new;
  array_2d_int n_Ev_death;
  array_2d_int n_Ev_to_Iv;
  
  // population of mosquitoes
  std::vector<Mosquito> mosq_pop;
  
  // store integer index of mosquitoes at various stages
  array_2d_int Sv_index;
  array_3d_int Ev_death;
  array_3d_int Ev_to_Iv;
  array_2d_int Iv_index;
  
  // objects for storing results
  array_3d_double daily_values;
  array_5d_int genotypes;
  array_4d_int indlevel_data;
  
  // misc
  std::vector<double> EIR;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher(Parameters &parameters, Rcpp::Function &update_progress, Rcpp::List &args_progress);
  
  // methods
  void simulate();
  
};