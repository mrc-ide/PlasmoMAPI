#include "main.h"

using namespace std;

//------------------------------------------------
// square a vector of values
Rcpp::List square_cpp(Rcpp::List args) {
  
  // print message to console
  Rcpp::Rcout << "running C++ dummy1_cpp function\n";
  
  // get inputs from Rcpp format to base C++ format
  vector<double> x = Rcpp::as<vector<double>>(args("x"));
  
  // square values
  for (int i = 0; i < int(x.size()); ++i) {
    x[i] *= x[i];
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x_squared") = x);
  return ret;
}
