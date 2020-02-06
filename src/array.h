
#pragma once

#include <vector>

//------------------------------------------------
// 2d array of integers
class array_2d_int {
public:
  
  // array and dimensions
  std::vector<std::vector<int>> arr;
  int d1;
  int d2;
  
  // constructors
  array_2d_int(int d1 = 0, int d2 = 0, int x = 0);
  
  // operators
  std::vector<int> & operator [](int i) {return arr[i];}
  
  // methods
  void init(int d1 = 0, int d2 = 0, int x = 0);
  void print();
};

//------------------------------------------------
// 3d array of integers
class array_3d_int {
public:
  
  // array and dimensions
  std::vector<std::vector<std::vector<int>>> arr;
  int d1;
  int d2;
  int d3;
  
  // constructors
  array_3d_int(int d1 = 0, int d2 = 0, int d3 = 0, int x = 0);
  
  // operators
  std::vector<std::vector<int>> & operator [](int i) {return arr[i];}
  
  // methods
  void init(int d1 = 0, int d2 = 0, int d3 = 0, int x = 0);
  void print();
};

//------------------------------------------------
// 3d array of doubles
class array_3d_double {
public:
  
  // array and dimensions
  std::vector<std::vector<std::vector<double>>> arr;
  int d1;
  int d2;
  int d3;
  
  // constructors
  array_3d_double(int d1 = 0, int d2 = 0, int d3 = 0, double x = 0.0);
  
  // operators
  std::vector<std::vector<double>> & operator [](int i) {return arr[i];}
  
  // methods
  void init(int d1 = 0, int d2 = 0, int d3 = 0, double x = 0.0);
  void print();
};

//------------------------------------------------
// 4d array of integers
class array_4d_int {
public:
  
  // array and dimensions
  std::vector<std::vector<std::vector<std::vector<int>>>> arr;
  int d1;
  int d2;
  int d3;
  int d4;
  
  // constructors
  array_4d_int(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int x = 0);
  
  // operators
  std::vector<std::vector<std::vector<int>>> & operator [](int i) {return arr[i];}
  
  // methods
  void init(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int x = 0);
};

//------------------------------------------------
// 5d array of integers
class array_5d_int {
public:
  
  // array and dimensions
  std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> arr;
  int d1;
  int d2;
  int d3;
  int d4;
  int d5;
  
  // constructors
  array_5d_int(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int d5 = 0, int x = 0);
  
  // operators
  std::vector<std::vector<std::vector<std::vector<int>>>> & operator [](int i) {return arr[i];}
  
  // methods
  void init(int d1 = 0, int d2 = 0, int d3 = 0, int d4 = 0, int d5 = 0, int x = 0);
};