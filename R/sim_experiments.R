#'@title Simulate experiments

#' @param timesteps number of iterations to run each simulation for
#' @param reps the number of replicates to run for each set of parameters
#' @param L the size of lattice
#' @param seed seed number for reproducibility

#' @import tibble


sim_experiments <- function(timesteps = 10, reps = 5, L, seed = 0.1){

  init <- expand.grid(x = 1:L, y = 1:L, N1 = 10, N2 = 10)%>%
    tibble::rowid_to_column('ID')




}
