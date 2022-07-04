#' @title Weighted frequencies

#' @description Weights neighborhood frequencies according to an exponential function

#' @param N densities of neighborhood cells
#' @param delta coefficient determining strength of frequency-dependence (higher freq patches receive more weights)

#' @return returns the weights (sums up to 1) of each cell
#' @export



weighted_freq <- function(N, delta){

  w <- exp(delta*(N+1))

  w_freq <- w/sum(w)

  return(w_freq)
}

