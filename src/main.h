
#pragma once

#include "misc_v9.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif


//------------------------------------------------
// assign edges to hexes based on intersection
// [[Rcpp::export]]
Rcpp::List assign_map_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);

//------------------------------------------------
// run main analysis
// [[Rcpp::export]]
Rcpp::List pm_analysis_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);

