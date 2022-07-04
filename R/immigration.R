#' @title Immigration

#' @description Distributes emigrants (via sampling with replacement) into neighborhood patches given their weights

#' @param hood coordinates of the neighborhood of emigrating patch
#' @param emigrants the number of individuals to distribute among the neighborhood
#' @param weights the weights of neighborhood patches for probablistic drawing


#' @return returns neighborhood coordinates and the number of times they were drawn
#' @export

immigration <- function(hood, emigrants, weights){

  #Samples neighborhood coordinates with replacement and weights
  hood%>%
    dplyr::sample_n(size = emigrants, weight = weights, replace = TRUE)%>%
    group_by(`x`,`y`)%>%
    tally()%>%
    ungroup()

}
