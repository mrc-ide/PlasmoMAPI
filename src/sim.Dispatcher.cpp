
#include "sim.Dispatcher.h"
#include "probability_v10.h"
#include "array.h"

#include <tuple>

using namespace std;

//------------------------------------------------
// constructor
Dispatcher::Dispatcher(Parameters &parameters, Rcpp::Function &update_progress, Rcpp::List &args_progress) {
  
  // pointers to inputs
  param_ptr = &parameters;
  update_progress_ptr = &update_progress;
  args_progress_ptr = &args_progress;
  
  // make local copies of some parameters
  n_demes = param_ptr->n_demes;
  max_time = param_ptr->max_time;
  a = param_ptr->a;
  v = param_ptr->v;
  prob_v_death = param_ptr->prob_v_death;
  H = param_ptr->H;
  infectivity = param_ptr->infectivity;
  mu = param_ptr->mu;
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(param_ptr->age_stable, 1000);
  sampler_age_death = Sampler(param_ptr->age_death, 1000);
  sampler_duration_infection = Sampler(param_ptr->duration_infection, 1000);
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  schedule_death = vector<set<int>>(max_time);
  schedule_Eh_to_Ih = vector<vector<pair<int, int>>>(max_time);
  schedule_Ih_to_Sh = vector<vector<pair<int, int>>>(max_time);
  schedule_infective = vector<vector<pair<int, int>>>(max_time);
  schedule_infective_recovery = vector<vector<pair<int, int>>>(max_time);
  
  // counts of host types
  H_total = n_demes*H;
  Sh = vector<int>(n_demes, H);
  Eh = vector<int>(n_demes);
  Ih = vector<int>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(H_total);
  next_host_ID = 0;
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only
  host_index = array_2d_int(n_demes);
  host_infective_index = vector<vector<int>>(n_demes);
  for (int k = 0; k < n_demes; ++k) {
    host_index[k] = seq_int(k*H, (k+1)*H-1);
    reshuffle(host_index[k]);
  }
  
  // initialise hosts
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < H; ++i) {
      int this_host = host_index[k][i];
      host_pop[this_host].init(this_host, next_host_ID, k, Sh, Eh, Ih, host_infective_index,
                               schedule_death, schedule_Eh_to_Ih, schedule_Ih_to_Sh,
                               schedule_infective, schedule_infective_recovery,
                               sampler_age_stable, sampler_age_death, sampler_duration_infection,
                               parameters);
    }
  }
  
  // seed initial infections
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < param_ptr->seed_infections[k]; i++) {
      host_pop[host_index[k][i]].denovo_infection();
    }
  }
  
  // counts of mosquito types
  M_total = sum(param_ptr->M_vec);
  Sv = param_ptr->M_vec;
  Ev = vector<int>(n_demes);
  Iv = vector<int>(n_demes);
  
  // number of mosquitoes at various stages
  n_Ev_death_new = vector<int>(v);
  n_Ev_death = array_2d_int(n_demes, v);
  n_Ev_to_Iv = array_2d_int(n_demes, v);
  
  // population of mosquitoes
  mosq_pop = vector<Mosquito>(M_total);
  
  // initialise mosquitoes
  for (int i = 0; i < M_total; ++i) {
    mosq_pop[i].init(param_ptr);
  }
  
  // store integer index of mosquitoes at various stages
  Sv_index = array_2d_int(n_demes);
  int M_cum = 0;
  for (int k = 0; k < n_demes; ++k) {
    Sv_index[k] = seq_int(M_cum, M_cum + param_ptr->M_vec[k] - 1);
    M_cum += param_ptr->M_vec[k];
  }
  Ev_death = array_3d_int(n_demes, v);
  Ev_to_Iv = array_3d_int(n_demes, v);
  Iv_index = array_2d_int(n_demes);
  
  // objects for storing results
  daily_values = array_3d_double(n_demes, max_time);
  genotypes = array_5d_int(n_demes, param_ptr->n_time_out);
  indlevel_data = array_4d_int(n_demes, param_ptr->n_time_out);
  
  // misc
  EIR = vector<double>(n_demes);
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  // start message
  if (param_ptr->report_progress) {
    print("Running simulation");
  }
  
  // initialise indices
  int ringtime = 0;
  int index_time_out = 0;
  
  // vector for randomly changing the order in which migration is applied
  vector<int> mig_order = seq_int(0, param_ptr->n_mig_list - 1);
  
  // loop through daily time steps
  for (int t = 0; t < max_time; ++t) {
    
    // report progress
    if (param_ptr->report_progress) {
      (*update_progress_ptr)(*args_progress_ptr, "pb", t + 1, max_time);
    }
    
    // skip all dynamics if t=0. Ensures that user-defined starting conditions
    // are the first values stored in output
    if (t > 0) {
      
      // update ring buffer index
      ringtime = (ringtime == v - 1) ? 0 : ringtime + 1;
      
      
      //-------- MIGRATION --------
  #ifdef FOO
      reshuffle(mig_order);
      for (int i = 0; i < param_ptr->n_mig_list; ++i) {
        
        // draw number of migrants and skip over if zero
        double prob_migration = get<2>(param_ptr->mig_list[mig_order[i]]);
        int n_migrants = rbinom1(H, prob_migration);
        if (n_migrants == 0) {
          continue;
        }
        
        // get demes to swap between
        int deme1 = get<0>(param_ptr->mig_list[mig_order[i]]);
        int deme2 = get<1>(param_ptr->mig_list[mig_order[i]]);
        
        // loop through number of migrants
        for (int j=0; j<n_migrants; ++j) {
          
          // draw migrant index in both demes
          int rnd1 = sample2(0,H-1);
          int rnd2 = sample2(0,H-1);
          int index1 = host_index[deme1][rnd1];
          int index2 = host_index[deme2][rnd2];
          
          // swap host indices
          host_index[deme1][rnd1] = index2;
          host_index[deme2][rnd2] = index1;
          
          // move infectives as needed
          vector<int>::iterator it = find(host_infective_index[deme1].begin(), host_infective_index[deme1].end(), index1);
          if (it != host_infective_index[deme1].end()) {
            int tmp1 = distance(host_infective_index[deme1].begin(), it);
            host_infective_index[deme2].push_back(host_infective_index[deme1][tmp1]);
            quick_erase(host_infective_index[deme1], tmp1);
          }
          it = find(host_infective_index[deme2].begin(), host_infective_index[deme2].end(), index2);
          if (it != host_infective_index[deme2].end()) {
            int tmp1 = distance(host_infective_index[deme2].begin(), it);
            host_infective_index[deme1].push_back(host_infective_index[deme2][tmp1]);
            quick_erase(host_infective_index[deme2], tmp1);
          }
          
          // update host properties
          host_pop[index1].deme = deme2;
          host_pop[index2].deme = deme1;
          
          // update deme counts
          if (host_pop[index1].get_n_asexual() == 0) {
            Sh[deme1]--;
            Sh[deme2]++;
          } else if (host_pop[index1].n_infected == 0) {
            Eh[deme1]--;
            Eh[deme2]++;
          } else {
            Ih[deme1]--;
            Ih[deme2]++;
          }
          if (host_pop[index2].get_n_asexual() == 0) {
            Sh[deme2]--;
            Sh[deme1]++;
          } else if (host_pop[index2].n_infected == 0) {
            Eh[deme2]--;
            Eh[deme1]++;
          } else {
            Ih[deme2]--;
            Ih[deme1]++;
          }
          
        }
      }  // end migration i loop
  #endif
      
      // loop through demes
      for (int k = 0; k < n_demes; ++k) {
        
        
        //-------- MOSQUITO EVENTS --------
        
        // carry out Ev death
        int Ev_death_size = Ev_death[k][ringtime].size();
        if (Ev_death_size > 0) {
          push_back_multiple(Sv_index[k], Ev_death[k][ringtime]);
          Ev_death[k][ringtime].clear();
          Ev[k] -= Ev_death_size;
          Sv[k] += Ev_death_size;
        }
        
        // carry out Iv death
        int Iv_death = rbinom1(Iv[k], prob_v_death);
        for (int i = 0; i < Iv_death; ++i) {
          int rnd1 = sample2(0, Iv[k]-1-i);
          int this_mosq = Iv_index[k][rnd1];
          mosq_pop[this_mosq].death();
          Sv_index[k].push_back(this_mosq);
          quick_erase(Iv_index[k], rnd1);
        }
        Sv[k] += Iv_death;
        Iv[k] -= Iv_death;
        
        // move Ev into Iv
        int Ev_to_Iv_size = Ev_to_Iv[k][ringtime].size();
        if (Ev_to_Iv_size > 0) {
          push_back_multiple(Iv_index[k], Ev_to_Iv[k][ringtime]);
          Ev_to_Iv[k][ringtime].clear();
          Ev[k] -= Ev_to_Iv_size;
          Iv[k] += Ev_to_Iv_size;
        }
        
        // draw number of new bites on infective hosts. Note, the rate is already
        // scaled by infectivity, meaning bites that do not lead to mosquito
        // infection are ignored automatically
        int n_host_infective = host_infective_index[k].size();                  // number of infective hosts in this deme
        double rate_v_infected = a*infectivity*n_host_infective/double(H);      // rate of mosquito biting infective host and becoming infected
        double prob_v_infected_or_death = 1 - exp(-(rate_v_infected + mu));     // probability of mosquito becoming infected or dying (competing hazards)
        double prob_v_infected = rate_v_infected/(rate_v_infected + mu);        // relative probability of mosquito becoming infected vs. dying
        int v_infected_or_death = rbinom1(Sv[k], prob_v_infected_or_death);     // number of susceptible mosquitoes becoming infected or dying
        int v_infected = rbinom1(v_infected_or_death, prob_v_infected);         // number of susceptible mosquitoes becoming infected
        
        // loop through infective bites
        for (int i = 0; i < v_infected; ++i) {
          
          // choose mosquito at random from susceptibles
          int rnd1 = sample2(0, Sv[k]-1);
          int this_mosq = Sv_index[k][rnd1];
          
          // drop from Sv
          quick_erase(Sv_index[k], rnd1);
          
          // the majority of new mosquito infections will die in lag phase. If so
          // then no need to store genotype, instead simply add to Ev_death object
          // to schedule move back to susceptible state (equivalent to death) at
          // future time point.
          int v_time_death = rgeom1(prob_v_death) + 1;
          if (v_time_death <= v) {
            
            // schedule death
            Ev_death[k][(ringtime+v_time_death) % v].push_back(this_mosq);
            
          } else {
            
            // choose host at random from infectives
            int rnd2 = sample2(0, host_infective_index[k].size() - 1);
            int this_host = host_infective_index[k][rnd2];
            
            // infect mosquito
            mosq_pop[this_mosq].new_infection(&host_pop[this_host]);
            
            // schedule move to Iv
            Ev_to_Iv[k][ringtime].push_back(this_mosq);
          }
          
          // update deme counts
          Sv[k]--;
          Ev[k]++;
          
        } // end loop through infective bites
        
        //-------- NEW HUMAN EVENTS --------
        
        // get number of new infectious bites on humans
        EIR[k] = a*Iv[k]/double(H);
        double prob_h_infectious_bite = 1 - exp(-EIR[k]);            // probability of new infectious bite per host
        int h_infectious_bite = rbinom1(H, prob_h_infectious_bite);  // total number of new infectious bites
        
        // apply new infectious bites
        for (int i = 0; i < h_infectious_bite; ++i) {
          
          // choose host at random
          int rnd1 = sample2(0, H - 1);
          int this_host = host_index[k][rnd1];
          
          // determine whether infectious bite is successful
          if (rbernoulli1(host_pop[this_host].get_prob_infection())) {
            
            // infect from random mosquito
            int rnd2 = sample2(0, Iv[k] - 1);
            int this_mosq = Iv_index[k][rnd2];
            host_pop[this_host].new_infection(mosq_pop[this_mosq], t);
          }
          
        }  // end loop through infectious bites
        
      } // end loop through demes
      
      
      //-------- SCHEDULED HUMAN EVENTS --------
      
      // scheduled deaths
      for (auto it = schedule_death[t].begin(); it != schedule_death[t].end(); ++it) {
        int this_host = *it;
        host_pop[this_host].death(next_host_ID, t);
      }
      // scheduled Eh to Ih
      for (auto it = schedule_Eh_to_Ih[t].begin(); it != schedule_Eh_to_Ih[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Eh_to_Ih(this_slot);
      }
      
      // scheduled Ih to Sh
      for (auto it = schedule_Ih_to_Sh[t].begin(); it != schedule_Ih_to_Sh[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].Ih_to_Sh(this_slot);
      }
      
      // scheduled become infective
      for (auto it = schedule_infective[t].begin(); it != schedule_infective[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].begin_infective(this_slot);
      }
      
      // scheduled infective recovery
      for (auto it = schedule_infective_recovery[t].begin(); it != schedule_infective_recovery[t].end(); ++it) {
        int this_host = it->first;
        int this_slot = it->second;
        host_pop[this_host].end_infective(this_slot);
      }
    }  // end of if statement for skipping first timestep
    
    
    //-------- STORE RESULTS --------
    
    // store daily values
    for (int k = 0; k < n_demes; ++k) {
      daily_values[k][t] = {double(Sh[k]), double(Eh[k]), double(Ih[k]),
                            double(Sv[k]), double(Ev[k]), double(Iv[k]),
                            EIR[k]};
    }
    
    // if one of output times
    if ((t + 1) == param_ptr->time_out[index_time_out]) {
      
      // store genotypes and individual-level data
      for (int i = 0; i < H_total; ++i) {
        if (host_pop[i].n_infected > 0) {
          int this_deme = host_pop[i].deme;
          
          // store bloodstage genotypes
          vector<vector<int>> tmp_mat;
          for (int j = 0; j < param_ptr->max_innoculations; ++j) {
            if (host_pop[i].innoc_status_asexual[j] == Bloodstage_asexual) {
              for (unsigned int i2 = 0; i2 < host_pop[i].haplotypes[j].size(); ++i2) {
                tmp_mat.push_back(host_pop[i].haplotypes[j][i2]);
              }
            }
          }
          genotypes.arr[this_deme][index_time_out].push_back(tmp_mat);
          
          // store individual-level data
          int this_age = floor((t - host_pop[i].birth_day)/double(365.0));
          vector<int> tmp_vec = {host_pop[i].ID, host_pop[i].home_deme, this_age, host_pop[i].n_infected};
          indlevel_data[this_deme][index_time_out].push_back(tmp_vec);
        }
      }
      
      // increment output time index
      index_time_out++;
    }
    
  } // end time loop
  
}
