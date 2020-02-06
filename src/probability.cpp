
#include "probability.h"

using namespace std;

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1() {
  return R::runif(0,1);
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(const double a, const double b) {
  return R::runif(a,b);
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(const double p) {
  return R::rbinom(1, p);
}

//------------------------------------------------
// draw from Binomial(n, p) distribution
int rbinom1(const int n, const double p) {
  if (n == 0 || p >= 1) {
    return n;
  }
  return R::rbinom(n, p);
}

//------------------------------------------------
// draw from Geometric(p) distribution, with mean (1-p)/p
int rgeom1(const double p) {
  return R::rgeom(p);
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal 
// probability. Works on positive or negative values of a or b, and works 
// irrespective of which of a or b is larger.
int sample2(const int a, const int b) {
  if (a < b) {
    return floor(runif1(a, b+1));
  } else {
    return floor(runif1(b, a+1));
  }
}

//------------------------------------------------
// take in vector of integers (by reference). Resample these without replacement
void reshuffle(vector<int> &x) {
  int i, rnd1, tmp1;
  int n = x.size();
  for (i = 0; i < n; i++) {
    rnd1 = runif1(i, n);	// draw random index from i to end of vector. Note that although runif returns a double, by forcing to int we essentially round this value down to nearest int.
    tmp1 = x[rnd1];		    // temporarily store current value of vector at this position
    x[rnd1] = x[i];			// swap for value at position i
    x[i] = tmp1;			// complete the swap
  }
}

//------------------------------------------------
// resample a vector of integers without replacement while always sampling from within specified groups
void reshuffle_group(vector<int> &x, vector<vector<int>> &group_list) {
  for (int j=0; j<int(group_list.size()); ++j) {
    int rnd1, tmp1, old_pos, new_pos;
    int n = group_list[j].size();
    for (int i=0; i<n; i++) {
      rnd1 = runif1(i, n);
      old_pos = group_list[j][i];
      new_pos = group_list[j][rnd1];
      tmp1 = x[new_pos];
      x[new_pos] = x[old_pos];
      x[old_pos] = tmp1;
    }
  }
}