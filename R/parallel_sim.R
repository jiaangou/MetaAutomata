#'@title Parallelized stochastic simulation

#' @description A function that calls `stochatic_sim()` to run the replicate simulations in parallel. Same parameters as `stochastic_sim()` but includes arguments for number of replicates and cores to use. Parallelization is done via the `doParallel` package.

#' @param initial_df dataframe specifying the initial conditions with ID, xy coordinates, speices denisities
#' @param r growth rate of species (assumed equivalent for all species)
#' @param aij interspecific interaction strength (aii assumed to be 1)
#' @param K carrying capacity of the environment (competition matrix is scaled according to K)
#' @param delta aggregation parameter where 0 is uniform weighting, >0 high densities are weighted more (<0 = less weight)
#' @param timesteps length of time to simulate dynamics
#' @param disp_rate dispersal rate which can range from 0 (no dispersal) to 1 (maximum dispersal)
#' @param nh_size Manhattan distance specifying size of neighborhood
#' @param nh Type of neighborhood for dispersal to occur. Can be either 'Von Neummann' or 'Moore'
#' @param dd_emi If TRUE, emigration is density-dependent and follows type-II functional response
#' @param torus If TRUE, the lattice is wrapped around to remove edge effects
#' @param n_cores Number of cores to use (default is 2)
#' @param reps Number of replicates to run
#' @param type Type of


#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import tibble
#' @import doParallel
#' @importFrom stats "setNames"
#'



parallel_sim <- function(reps, n_cores = 2, initial_df, aij, delta, r = 3, K = 100, timesteps = 5, disp_rate = 0.1, nh_size = 1,
                         nh = "vonNeumann", dd_emi = FALSE,  torus = TRUE, type = "FORK"){

  #library(doParallel)
  cl <- makeCluster(n_cores, type = type) #FORK type so env does not need to be copied

  registerDoParallel(cl = cl)

  #experiments
  experiments <-  foreach(i = 1:reps) %dopar% {
    set.seed(i);
    stochastic_sim(initial_df, aij, delta, r = 3, K = 100, timesteps = 5, disp_rate = 0.1, nh_size = 1,
                   nh = "vonNeumann", dd_emi = FALSE,  torus = TRUE)
  }
  stopCluster(cl)

  return(experiments)

}
