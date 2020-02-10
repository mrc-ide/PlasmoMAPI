
#include "sim.Parameters.h"
#include "misc_v9.h"

using namespace std;

//------------------------------------------------
// constructor
Parameters::Parameters(const Rcpp::List &args) {
  
  // genetic parameters
  L = rcpp_to_int(args["L"]);
  prob_cotransmission = rcpp_to_double(args["prob_cotransmission"]);
  
  // epidemiological parameters
  a = rcpp_to_double(args["a"]);
  mu = rcpp_to_double(args["mu"]);
  u = rcpp_to_int(args["u"]);
  v = rcpp_to_int(args["v"]);
  g = rcpp_to_int(args["g"]);
  prob_infection = rcpp_to_vector_double(args["prob_infection"]);
  n_prob_infection = int(prob_infection.size());
  duration_infection = rcpp_to_vector_double(args["duration_infection"]);
  n_duration_infection = int(duration_infection.size());
  infectivity = rcpp_to_double(args["infectivity"]);
  max_innoculations = rcpp_to_int(args["max_innoculations"]);
  
  // deme parameters
  H = rcpp_to_int(args["H"]);
  seed_infections = rcpp_to_vector_int(args["seed_infections"]);
  M_vec = rcpp_to_vector_int(args["M"]);
  n_demes = int(M_vec.size());
  
  // get migration list from matrix
  vector<vector<double>> mig_matrix = rcpp_to_matrix_double(args["mig_matrix"]);
  for (int i=0; i<(n_demes-1); ++i) {
    for (int j=(i+1); j<n_demes; ++j) {
      if (mig_matrix[i][j] > 0) {
        mig_list.push_back(make_tuple(i, j, mig_matrix[i][j]));
      }
    }
  }
  n_mig_list = mig_list.size();
  
  // demography
  life_table = rcpp_to_vector_double(args["life_table"]);
  age_death = rcpp_to_vector_double(args["age_death"]);
  age_stable = rcpp_to_vector_double(args["age_stable"]);
  n_age = int(age_stable.size());
  
  // run parameters
  time_out = rcpp_to_vector_int(args["time_out"]);
  n_time_out = int(time_out.size());
  max_time = max(time_out);
  report_progress = rcpp_to_bool(args["report_progress"]);
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  
}

//------------------------------------------------
// print summary
void Parameters::print_summary() {
  
  // genetic parameters
  print("-- genetic parameters --");
  print("L: ", L);
  print("prob_cotransmission: ", prob_cotransmission);
  
  // epidemiological parameters
  print("-- epidemiological parameters --");
  print("a: ", a);
  print("mu: ", mu);
  print("u: ", u);
  print("v: ", v);
  print("g: ", g);
  print("prob_infection: ");
  print_vector(prob_infection);
  print("duration_infection: ");
  print_vector(duration_infection);
  print("infectivity: ", infectivity);
  print("max_innoculations: ", max_innoculations);
  print("");
  
  // deme parameters
  print("-- deme parameters --");
  print("H: ", H);
  print("seed_infections: ");
  print_vector(seed_infections);
  print("M_vec: ");
  print_vector(M_vec);
  print("n_demes: ", n_demes);
  print("");
  
  // demography
  print("-- demography --");
  print("life_table: ");
  print_vector(life_table);
  print("age_death: ");
  print_vector(age_death);
  print("age_stable: ");
  print_vector(age_stable);
  print("n_age: ", n_age);
  print("");
  
  // migration
  print("-- migration --");
  for (int i=0; i<n_mig_list; ++i) {
    print(get<0>(mig_list[i]), get<1>(mig_list[i]), get<2>(mig_list[i]));
  }
  print("");
  
  // run parameters
  print("-- run parameters --");
  print("t_out: ");
  print_vector(time_out);
  print("n_time_out: ", n_time_out);
  print("max_time: ", max_time);
  print("");
  
  // misc parameters
  print("-- misc parameters --");
  print("prob_v_death: ", prob_v_death);
  print("");
}
