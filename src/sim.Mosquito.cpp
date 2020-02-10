
#include "sim.Mosquito.h"
#include "sim.Host.h"
#include "probability_v10.h"
#include "misc_v9.h"

using namespace std;

//------------------------------------------------
// constructors
void Mosquito::init(Parameters* param_ptr) {
  L = param_ptr->L;
  max_innoculations = param_ptr->max_innoculations;
  n_haplotypes = 0;
  n_products = 0;
}

//------------------------------------------------
// copy over genotype from infectious host
void Mosquito::new_infection(Host* host_ptr) {
  
  // choose two haplotypes at random to copy from host
  int rnd1 = sample2(1,host_ptr->n_infective_haplotypes_total);
  int rnd2 = sample2(1,host_ptr->n_infective_haplotypes_total);
  
  // copy first haplotype from host
  int tmp1 = 0;
  for (int i=0; i<max_innoculations; ++i) {
    tmp1 += host_ptr->n_infective_haplotypes[i];
    if (tmp1 >= rnd1) {
      haplotype1 = host_ptr->haplotypes[i][tmp1-rnd1];
      n_haplotypes++;
      break;
    }
  }
  
  // if second copy is identical to first then do not bother copying over
  if (rnd1 == rnd2) {
    return;
  }
  
  // otherwise copy second haplotype from host
  tmp1 = 0;
  for (int i=0; i<max_innoculations; ++i) {
    tmp1 += host_ptr->n_infective_haplotypes[i];
    if (tmp1 >= rnd2) {
      haplotype2 = host_ptr->haplotypes[i][tmp1-rnd2];
      n_haplotypes++;
      break;
    }
  }
  
  return;
}

//------------------------------------------------
// initialise with single random haplotype
void Mosquito::denovo_infection() {
  haplotype1 = vector<int>(L);
  for (int i=0; i<L; ++i) {
    haplotype1[i] = rbernoulli1(0.5) ? 1 : 0;
  }
  n_haplotypes = 1;
}

//------------------------------------------------
// death
void Mosquito::death() {
  haplotype1.clear();
  haplotype2.clear();
  n_haplotypes = 0;
  products.clear();
  n_products = 0;
}

//------------------------------------------------
// get product of recombination
vector<int> Mosquito::get_product() {
  
  // if single haplotype then return this clone
  if (n_haplotypes == 1) {
    return haplotype1;
  }
  
  // otherwise determine whether we need a novel product
  bool new_product = rbernoulli1(1 - 0.25*n_products);
  
  // if so then produce new product and return
  if (new_product) {
    
    // recombination
    products.push_back(haplotype1);
    for (int i=0; i<L; ++i) {
      if (rbernoulli1(0.5)) {
        products[n_products][i] = haplotype2[i];
      }
    }
    n_products++;
    return products[n_products-1];
  }
  
  // otherwise return existing product, chosen at random
  int rnd1 = sample2(0,n_products-1);
  return products[rnd1];
  
}
