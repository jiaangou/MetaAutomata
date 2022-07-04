#'@title Beverton-Holt competition

#' @param r growth rate of species
#' @param alpha competition coefficients
#' @param initial_N initial abundances
#' @param timesteps length of time to simulate dynamics
#' @param discrete whether or not abundances should be discrete; discretization leads to a stochastic model
#' @param final_N whether or not function should return only the final N after all timesteps
#' @export

bh_competition <- function(r = 1.5, alpha = 1, initial_N = 10, timesteps = 30, discrete = FALSE, final_N = FALSE, extinction_thresh = 0.01){

  # #Check arguments
  # if(!(length(r) == ncol(alpha)) == (nrow(alpha) == length(initial_N))){
  #   stop(('check argument dimensions'))
  # }

  #Number of species
  num_spp <- length(r)

  #Initiate empty matrix of length timesteps + 1
  Nt <- matrix(0, nrow = timesteps + 1, ncol = num_spp)

  #print(dim(Nt))

  #Input initial N input first value
  Nt[1,] <- initial_N

  #Iterate to populate Nt vector
  for (t in 1:timesteps){

    if(discrete == TRUE){
      N <- Nt[t,]
      #draw from poisson distr with lambda = deterministic N (these are the next gen individuals)
      lambda_N <- N * r / (1 + alpha%*%N)
      lambda_N[lambda_N < 0] <- 0 # set lambda = 0 if negative

      Nt[t+1,] <- rpois(num_spp, lambda = lambda_N) #draw from distribution

    }else{
      N <- Nt[t,]
      N[N < extinction_thresh] <- 0 #Set N to 0 if N below extinction threshold
      Nt[t+1,] <- N * r / (1 + alpha%*%N)
    }

  }

  if(final_N == TRUE){
    return(Nt[nrow(Nt),])
  }else{
    return(Nt)
  }


}
