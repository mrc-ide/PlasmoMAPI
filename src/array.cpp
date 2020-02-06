
#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <iostream>

#include "array.h"
#include "misc_v4.h"

using namespace std;


// ==== TWO-DIMENSIONAL ARRAY OF INTS ============

//------------------------------------------------
// constructors
array_2d_int::array_2d_int(int d1, int d2, int x) {
  init(d1, d2, x);
}

//------------------------------------------------
// initialise
void array_2d_int::init(int d1, int d2, int x) {
  this->d1 = d1;
  this->d2 = d2;
  arr = vector<vector<int>>(d1, vector<int>(d2,x));
}

//------------------------------------------------
// print
void array_2d_int::print() {
  print_matrix(arr);
}


// ==== THREE-DIMENSIONAL ARRAY OF INTS ==========

//------------------------------------------------
// constructors
array_3d_int::array_3d_int(int d1, int d2, int d3, int x) {
  init(d1, d2, d3, x);
}

//------------------------------------------------
// initialise
void array_3d_int::init(int d1, int d2, int d3, int x) {
  this->d1 = d1;
  this->d2 = d2;
  this->d3 = d3;
  arr = vector<vector<vector<int>>>(d1, vector<vector<int>>(d2, vector<int>(d3, x)));
}

//------------------------------------------------
// print
void array_3d_int::print() {
  print_array(arr);
}


// ==== THREE-DIMENSIONAL ARRAY OF DOUBLES =======

//------------------------------------------------
// constructors
array_3d_double::array_3d_double(int d1, int d2, int d3, double x) {
  init(d1, d2, d3, x);
}

//------------------------------------------------
// initialise
void array_3d_double::init(int d1, int d2, int d3, double x) {
  this->d1 = d1;
  this->d2 = d2;
  this->d3 = d3;
  arr = vector<vector<vector<double>>>(d1, vector<vector<double>>(d2, vector<double>(d3, x)));
}

//------------------------------------------------
// print
void array_3d_double::print() {
  print_array(arr);
}


// ==== FOUR-DIMENSIONAL ARRAY OF INTS ===========

//------------------------------------------------
// constructors
array_4d_int::array_4d_int(int d1, int d2, int d3, int d4, int x) {
  init(d1, d2, d3, d4, x);
}

//------------------------------------------------
// initialise
void array_4d_int::init(int d1, int d2, int d3, int d4, int x) {
  this->d1 = d1;
  this->d2 = d2;
  this->d3 = d3;
  this->d4 = d4;
  arr = vector<vector<vector<vector<int>>>>(d1, vector<vector<vector<int>>>(d2, vector<vector<int>>(d3, vector<int>(d4, x))));
}


// ==== FIVE-DIMENSIONAL ARRAY OF INTS ===========

//------------------------------------------------
// constructors
array_5d_int::array_5d_int(int d1, int d2, int d3, int d4, int d5, int x) {
  init(d1, d2, d3, d4, d5, x);
}

//------------------------------------------------
// initialise
void array_5d_int::init(int d1, int d2, int d3, int d4, int d5, int x) {
  this->d1 = d1;
  this->d2 = d2;
  this->d3 = d3;
  this->d4 = d4;
  this->d5 = d5;
  arr = vector<vector<vector<vector<vector<int>>>>>(d1, vector<vector<vector<vector<int>>>>(d2, vector<vector<vector<int>>>(d3, vector<vector<int>>(d4, vector<int>(d5, x)))));
}
