
#------------------------------------------------
#' @title Simulate genetic data from simple P. falciparum model
#'
#' @description Simulate genetic data from a simple model of P. falciparum
#'   epidemiology and genetics.
#'
#' @param L number of loci. The maximum number of loci is 1000, as at higher
#'   numbers haplotypes begin to exceed integer representation (2^L).
#' @param prob_cotransmission probability of mosquito transmitting multiple
#'   haplotypes to host.
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless otherwise
#'   specified.
#' @param u intrinsic incubation period. The number of days from infection to
#'   blood-stage infection in a human host.
#' @param v extrinsic incubation period. The number of days from infection to
#'   becoming infectious in a mosquito.
#' @param g lag time between human blood-stage infection and production of
#'   gametocytes.
#' @param prob_infection probability a human becomes infected after being bitten
#'   by an infected mosquito.
#' @param duration_infection vector specifying probability distribution of time
#'   (in days) of a malaria episode.
#' @param infectivity probability a mosquito becomes infected after biting an
#'   infective human host.
#' @param max_innoculations maximum number of innoculations that an individual
#'   can hold simultaneously.
#' @param H human population size, which is assumed the same in every deme.
#' @param seed_infections vector specifying the initial number of infected
#'   humans in each deme.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#' @param mig_matrix migration matrix specifing the daily probability of
#'   migrating from each deme to each other deme. Migration must be equal in
#'   both directions, meaning this matrix must be symmetric.
#' @param time_out vector of times (days) at which output is produced.
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console during simuation.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats dgeom
#' @export

sim_falciparum <- function(a = 0.3,
                           p = 0.9,
                           mu = -log(p),
                           u = 12,
                           v = 10,
                           g = 10,
                           prob_infection = seq(0.1,0.01,-0.01),
                           duration_infection = dgeom(1:500, 1/200),
                           infectivity = 1,
                           max_innoculations = 5,
                           H = 1000,
                           seed_infections = 100,
                           M = 1000,
                           mig_matrix = diag(length(M)),
                           L = 24,
                           prob_cotransmission = 0.5,
                           time_out = 100,
                           report_progress = TRUE) {
  
  # check inputs
  assert_single_bounded(a)
  assert_single_bounded(p)
  assert_single_pos(mu)
  assert_single_pos_int(u, zero_allowed = FALSE)
  assert_single_pos_int(v, zero_allowed = FALSE)
  assert_single_pos_int(g, zero_allowed = FALSE)
  assert_vector(prob_infection)
  assert_bounded(prob_infection)
  assert_pos(prob_infection[1], zero_allowed = FALSE)
  assert_vector(duration_infection)
  assert_pos(duration_infection, zero_allowed = TRUE)
  assert_single_bounded(infectivity)
  assert_single_pos_int(max_innoculations, zero_allowed = FALSE)
  assert_single_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(seed_infections, zero_allowed = FALSE)
  assert_leq(seed_infections, H)
  assert_pos_int(M, zero_allowed = FALSE)
  assert_same_length(M, seed_infections)
  n_demes <- length(M)
  #assert_symmetric_matrix(mig_matrix)
  assert_bounded(mig_matrix)
  assert_dim(mig_matrix, c(n_demes, n_demes))
  #assert_eq(rowSums(mig_matrix), rep(1,n_demes))
  assert_single_pos_int(L, zero_allowed = FALSE)
  assert_leq(L, 1000)
  assert_single_bounded(prob_cotransmission)
  assert_vector(time_out)
  assert_pos_int(time_out, zero_allowed = TRUE)
  assert_single_logical(report_progress)
  
  # normalise infection duration distribution
  duration_infection <- duration_infection/sum(duration_infection)
  
  # if any prob_infection is zero then this automatically defines the value of
  # max_innoculations
  if (any(prob_infection == 0)) {
    max_innoculations <- min(max_innoculations, which(prob_infection == 0)[1] - 1)
  }
  
  # read in Mali demography distribution
  mali_demog <- pm_file("mali_demog.rds")
  
  # ---------------------------------------------
  # set up arguments for input into C++
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, max(time_out), initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create argument list
  args <- list(a = a,
               mu = mu,
               u = u,
               v = v,
               g = g,
               prob_infection = prob_infection,
               duration_infection = duration_infection,
               infectivity = infectivity,
               max_innoculations = max_innoculations,
               H = H,
               seed_infections = seed_infections,
               M = M,
               mig_matrix = matrix_to_rcpp(mig_matrix),
               L = L,
               prob_cotransmission = prob_cotransmission,
               life_table = mali_demog$life_table,
               age_death = mali_demog$age_death,
               age_stable = mali_demog$age_stable,
               time_out = time_out,
               report_progress = report_progress)
  
  # ---------------------------------------------
  # run efficient C++ function
  
  output_raw <- sim_falciparum_cpp(args, args_functions, args_progress)
  
  # ---------------------------------------------
  # process raw output
  
  message("processing output")
  
  # get daily values
  daily_values <- mapply(function(x) {
    ret <- rcpp_to_matrix(x)
    ret <- cbind(1:nrow(ret), ret)
    ret <- as.data.frame(ret)
    names(ret) <- c("time", "S", "E", "I", "Sv", "Ev", "Iv", "EIR")
    return(ret)
  }, output_raw$daily_values, SIMPLIFY = FALSE)
  names(daily_values) <- sprintf("deme%s", 1:n_demes)
  
  # get individual level data
  indlevel <- mapply(function(y) {
    ret <- mapply(function(x) {
      if (length(x) == 0) {
        return(NULL)
      } else {
        ret <- as.data.frame(rcpp_to_matrix(x))
        names(ret) <- c("ID", "home_deme", "age", "n_innoculations")
        return(ret)
      }
    }, y, SIMPLIFY = FALSE)
    names(ret) <- paste0("time", time_out)
    return(ret)
  }, output_raw$indlevel_data, SIMPLIFY = FALSE)
  names(indlevel) <- paste0("deme", 1:n_demes)
  
  # add haplotypes to indlevel data
  powers_two <- 2^((L-1):0)
  for (i in 1:n_demes) {
    for (j in 1:length(time_out)) {
      if (!is.null(indlevel[[i]][[j]])) {
        
        # get haplotypes in numeric form (converted from binary)
        haplotype_ID <- mapply(function(y) {
          ret <- mapply(function(x) {
            sum(x*powers_two)
          }, y)
          ret <- unique(ret)
          return(ret)
        }, output_raw$genotypes[[i]][[j]], SIMPLIFY = FALSE)
        
        # store haplotypes and number of haplotypes
        indlevel[[i]][[j]]$haplotype_ID <- haplotype_ID
        indlevel[[i]][[j]]$n_haplotypes <- mapply(length, haplotype_ID)
      }
    }
  }
  
  # create custom class and return
  ret <- list(daily_values = daily_values,
              indlevel = indlevel)
  class(ret) <- "pm_sim"
  return(ret)
}
