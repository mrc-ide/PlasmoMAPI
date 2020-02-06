
#pragma once

#include "sim.Parameters.h"

#include <vector>

//------------------------------------------------
// forward-declare Host class
class Host;

//------------------------------------------------
// class defining mosquito
class Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  
  // parameters copied over
  int L;
  int max_innoculations;
  
  // genotypes
  std::vector<int> haplotype1;
  std::vector<int> haplotype2;
  int n_haplotypes;
  std::vector<std::vector<int>> products;
  int n_products;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito() {};
  
  // methods
  void init(Parameters* param_ptr);
  void new_infection(Host* host_ptr);
  void denovo_infection();
  void death();
  std::vector<int> get_product();
  
};
