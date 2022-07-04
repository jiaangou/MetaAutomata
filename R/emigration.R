#'@title Density-(in)dependent Emigration

#' @param N number of individuals
#' @param max saturation value of emigration rate/probability
#' @param a half-saturation constant that controls rate of increase
#' @param threshold minimum value of N in which emigration is possible
#' @param discrete whether or not to return emigrants as discrete individuals, setting discrete == TRUE also makes the model stochastic
#' @param dd whether or not to implement density-dependence, if FALSE then a constant emigration rate (i.e. max rate)
#' @importFrom stats "rbinom"
#' @export


dd_emigration <- function(N, max = 0.3, a = 30, threshold =  0, discrete = FALSE, dd = FALSE){

  if(dd == TRUE){

    #Density-dependent
    N_new <- N - threshold
    h <- max
    EmiProb <- h*N_new/(a+N_new)
    #Set the ones lower than threshold to 0
    if(sum(N < threshold)>0){
      EmiProb[which(N < threshold)] <- 0
    }

    #Density-independent (constant emigration rate = max)
  }else{
    EmiProb <- max
  }


  #Discrete/Stochastic
  #If discrete, return integers
  if(discrete == TRUE){

    y <- rbinom(n = length(N), size = N, prob = EmiProb)

  }else{
    y <- N*EmiProb
  }


  return(y)
}
