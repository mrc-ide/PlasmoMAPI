
#include "sim.Host.h"
#include "probability_v10.h"
#include "misc_v9.h"

using namespace std;

//------------------------------------------------
// initialise host
void Host::init(int index, int &ID, int deme,
                vector<int> &Sh, vector<int> &Eh, vector<int> &Ih,
                vector<vector<int>> &host_infective_index,
                vector<set<int>> &schedule_death,
                vector<vector<pair<int, int>>> &schedule_Eh_to_Ih,
                vector<vector<pair<int, int>>> &schedule_Ih_to_Sh,
                vector<vector<pair<int, int>>> &schedule_infective,
                vector<vector<pair<int, int>>> &schedule_infective_recovery,
                Sampler &sampler_age_stable, Sampler &sampler_age_death, Sampler &sampler_duration_infection,
                Parameters &parameters) {
  
  // identifiers
  this->index = index;
  this->ID = ID++;
  home_deme = deme;
  this->deme = deme;
  
  // pointers
  Sh_ptr = &Sh;
  Eh_ptr = &Eh;
  Ih_ptr = &Ih;
  schedule_death_ptr = &schedule_death,
  schedule_Eh_to_Ih_ptr = &schedule_Eh_to_Ih;
  schedule_Ih_to_Sh_ptr = &schedule_Ih_to_Sh;
  schedule_infective_ptr = &schedule_infective;
  schedule_infective_recovery_ptr = &schedule_infective_recovery;
  host_infective_index_ptr = &host_infective_index;
  age_stable_ptr = &sampler_age_stable;
  age_death_ptr = &sampler_age_death;
  duration_infection_ptr = &sampler_duration_infection;
  param_ptr = &parameters;
  
  // copy over some parameters for convenience
  L = param_ptr->L;
  max_innoculations = param_ptr->max_innoculations;
  n_age = param_ptr->n_age;
  max_time = param_ptr->max_time;
  u = param_ptr->u;
  g = param_ptr->g;
  prob_cotransmission = param_ptr->prob_cotransmission;
  
  // indices relating to global distributions
  prob_infection_index = 0;
  
  // draw age from demography distribution
  int age_years = age_stable_ptr->draw();
  int extra_days = sample2(1, 365);
  int age_days = age_years*365 + extra_days;
  
  // draw duration of life from demography distribution looking forward from
  // current age. This is tricky, as we must account for the fact that if we
  // are already part way into an age group we have a reduced probability of
  // dying within that age group.
  int life_days = 0;
  double prop_year_remaining = 1 - extra_days/365.0;
  double prob_die_this_year = param_ptr->life_table[age_years]*prop_year_remaining;
  if (rbernoulli1(prob_die_this_year) || age_years == (n_age-1)) {
    life_days = age_years*365 + sample2(extra_days, 365);
  } else {
    for (int i=(age_years+1); i<n_age; ++i) {
      if (rbernoulli1(param_ptr->life_table[i])) {
        life_days = i*365 + sample2(1, 365);
        break;
      }
    }
  }
  
  // convert to final birth and death days
  birth_day = -age_days;
  death_day = life_days - age_days;
  if (death_day == 0) {  // in the unlikely even that due to die on day zero, delay death by one day
    death_day++;
  }
  
  // add death_day to scheduler
  if (death_day <= max_time) {
    (*schedule_death_ptr)[death_day].insert(index);
  }
  
  // initialise innoculation objects
  innoc_active = vector<bool>(max_innoculations, false);
  innoc_status_asexual = vector<Status_asexual>(max_innoculations, Inactive_asexual);
  innoc_status_sexual = vector<Status_sexual>(max_innoculations, Inactive_sexual);
  
  // initiliase haplotypes
  haplotypes = vector<vector<vector<int>>>(max_innoculations);
  n_infective_haplotypes = vector<int>(max_innoculations);
  n_infective_haplotypes_total = 0;
  
  // initialise innoculation counts
  n_latent = 0;
  n_infected = 0;
  n_infective = 0;
  
}

//------------------------------------------------
// death
void Host::death(int &ID, int t) {
  
  // drop from infective list if necessary
  if (n_infective > 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
  // update deme counts to reflect death
  if (n_infected > 0) {
    (*Ih_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  } else if (n_latent > 0) {
    (*Eh_ptr)[deme]--;
    (*Sh_ptr)[deme]++;
  }
  
  // new unique ID
  this->ID = ID++;
  
  // make current deme home deme
  home_deme = deme;
  
  // reset indices relating to global distributions
  prob_infection_index = 0;
  
  // date of birth
  this->birth_day = t;
  
  // draw life duration from demography distribution
  int life_years = age_death_ptr->draw();
  int life_days = life_years*365 + sample2(1, 365);
  death_day = birth_day + life_days;
  
  // add new death_day to scheduler
  if (death_day <= max_time) {
    (*schedule_death_ptr)[death_day].insert(index);
  }
  
  // reset innoculation objects
  fill(innoc_active.begin(), innoc_active.end(), false);
  fill(innoc_status_asexual.begin(), innoc_status_asexual.end(), Inactive_asexual);
  fill(innoc_status_sexual.begin(), innoc_status_sexual.end(), Inactive_sexual);
  
  // reset haplotypes
  haplotypes = vector<vector<vector<int>>>(max_innoculations);
  fill(n_infective_haplotypes.begin(), n_infective_haplotypes.end(), 0);
  n_infective_haplotypes_total = 0;
  
  // reset innoculation counts
  n_latent = 0;
  n_infected = 0;
  n_infective = 0;
  
}

//------------------------------------------------
// new infection
void Host::new_infection(Mosquito &mosq, int t) {
  
  // update prob_infection_index irrespective of whether infection takes hold
  update_prob_infection();
  
  // return if already at max_innoculations
  if (get_n_innoculations() == max_innoculations) {
    return;
  }
  
  // update deme counts
  if (get_n_asexual() == 0) {
    (*Sh_ptr)[deme]--;
    (*Eh_ptr)[deme]++;
  }
  
  // update counts
  n_latent++;
  
  // get next free innoculation slot
  int this_slot = 0;
  for (int i=0; i<max_innoculations; ++i) {
    if (!innoc_active[i]) {
      break;
    }
    this_slot++;
  }
  if (this_slot == max_innoculations) {
    Rcpp::stop("could not find free innoculation slot");
  }
  
  // add new innoculation
  innoc_active[this_slot] = true;
  innoc_status_asexual[this_slot] = Liverstage_asexual;
  
  // copy over products of recombination. If mosquito holds a single haplotype
  // then copy this over (clonal expansion). Otherwise copy over potentially
  // multiple recombinant haplotypes.
  if (mosq.n_haplotypes == 1) {
    haplotypes[this_slot].emplace_back(mosq.get_product());
  } else {
    double p = 1.0;
    for (int i = 0; i < 4; ++i) {
      if (rbernoulli1(p)) {
        haplotypes[this_slot].emplace_back(mosq.get_product());
      } else {
        break;
      }
      p *= prob_cotransmission;
    }
  }
  
  // draw duration of infection
  int duration_infection = duration_infection_ptr->draw() + 1;
  
  // get times of future events
  int t1 = t + u;                           // begin bloodstage
  int t2 = t + u + duration_infection;      // end bloodstage
  int t3 = t + u + g;                       // begin infective
  int t4 = t + u + g + duration_infection;  // end infective
  
  // schedule move to Ih
  if (t1 < death_day && t1 <= max_time) {
    (*schedule_Eh_to_Ih_ptr)[t1].emplace_back(index, this_slot);
  }
  
  // schedule bloodstage recovery
  if (t2 < death_day && t2 <= max_time) {
    (*schedule_Ih_to_Sh_ptr)[t2].emplace_back(index, this_slot);
  }
  
  // schedule begin infective
  if (t3 < death_day && t3 <= max_time) {
    (*schedule_infective_ptr)[t3].emplace_back(index, this_slot);
  }
  
  // schedule end infective
  if (t4 < death_day && t4 <= max_time) {
    (*schedule_infective_recovery_ptr)[t4].emplace_back(index, this_slot);
  }
  
}

//------------------------------------------------
// de-novo infection
void Host::denovo_infection() {
  
  // generating starting genotype in a dummy mosquito
  Mosquito dummy_mosquito;
  dummy_mosquito.init(param_ptr);
  dummy_mosquito.denovo_infection();
  
  // carry out infection
  new_infection(dummy_mosquito, 0);
}

//------------------------------------------------
// move from Eh state to Ih
void Host::Eh_to_Ih(int this_slot) {
  
  // update deme counts
  if (n_infected == 0) {
    (*Eh_ptr)[deme]--;
    (*Ih_ptr)[deme]++;
  }
  
  // update status
  innoc_status_asexual[this_slot] = Bloodstage_asexual;
  
  // update counts
  n_latent--;
  n_infected++;
  
}

//------------------------------------------------
// move from Ih state to Sh
void Host::Ih_to_Sh(int this_slot) {
  
  // update deme counts
  if (n_infected == 1) {
    (*Ih_ptr)[deme]--;
    if (n_latent == 0) {
      (*Sh_ptr)[deme]++;
    } else {
      (*Eh_ptr)[deme]++;
    }
  }
  
  // update status
  innoc_status_asexual[this_slot] = Inactive_asexual;
  
  // update counts
  n_infected--;
}

//------------------------------------------------
// begin infective period
void Host::begin_infective(int this_slot) {
  
  // update host status
  innoc_status_sexual[this_slot] = Active_sexual;
  
  // update host counts
  n_infective++;
  
  // update haplotype counts
  n_infective_haplotypes[this_slot] = haplotypes[this_slot].size();
  n_infective_haplotypes_total += n_infective_haplotypes[this_slot];
  
  // if newly infective then add to infectives list
  if (n_infective == 1) {
    (*host_infective_index_ptr)[deme].push_back(index);
  }
  
}

//------------------------------------------------
// end infective period
void Host::end_infective(int this_slot) {
  
  // update host status
  innoc_status_sexual[this_slot] = Inactive_sexual;
  innoc_active[this_slot] = false;
  
  // update host counts
  n_infective--;
  
  // clear heplotypes
  haplotypes[this_slot].clear();
  n_infective_haplotypes_total-= n_infective_haplotypes[this_slot];
  n_infective_haplotypes[this_slot] = 0;
  
  // if no longer infective then drop from infectives list
  if (n_infective == 0) {
    erase_remove((*host_infective_index_ptr)[deme], index);
  }
  
}

//------------------------------------------------
// update probabilty of infection
void Host::update_prob_infection() {
  if (prob_infection_index < (param_ptr->n_prob_infection-1)) {
    prob_infection_index++;
  }
}

//------------------------------------------------
// get total number of innoculations
int Host::get_n_innoculations() {
  return sum_bool(innoc_active);
}

//------------------------------------------------
// get total number of asexual stage innoculations
int Host::get_n_asexual() {
  return n_latent + n_infected;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  return param_ptr->prob_infection[prob_infection_index];
}
