
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(const double a = 0, const double b = 1.0);

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(const double p);

//------------------------------------------------
// draw from Binomial(n, p) distribution
int rbinom1(const int n, const double p);

//------------------------------------------------
// draw from Geometric(p) distribution
int rgeom1(const double p);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal 
// probability. Works on positive or negative values of a or b, and works 
// irrespective of which of a or b is larger.
int sample2(const int a, const int b);

//------------------------------------------------
// take in vector of integers (by reference). Resample these without replacement
void reshuffle(std::vector<int> &x);

//------------------------------------------------
// resample a vector of integers without replacement while always sampling from within specified groups
void reshuffle_group(std::vector<int> &x, std::vector<std::vector<int>> &group_list);