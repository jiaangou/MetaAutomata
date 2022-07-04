#'@title Stochastic metacommunity simulation

#' @description Pieces functions together and runs a stochastic simulation of metacommunity dynamics

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


#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import tibble
#' @importFrom stats "setNames"


#Silences error message for using the '.' syntax
#if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "x", "y","ID", "hood", "emigrants", "Species", "competition", "immigrants", "time"))


stochastic_sim <- function(initial_df, aij, delta, r = 3, K = 100, timesteps = 5, disp_rate = 0.1, nh_size = 1,
                           nh = "vonNeumann", dd_emi = FALSE,  torus = TRUE){

  #
  no.spp <- initial_df%>%
    select(-ID, -x, -y)%>%
    ncol()


  #Species names
  sp_names <- paste0('N', 1:no.spp)


  #Construct alpha matrix
  alpha <- matrix(aij, ncol = no.spp, nrow = no.spp)
  diag(alpha) <- 1
  alpha <- alpha/(K/(r-1)) #scale coefficients by env carrying capacity

  #
  r <- rep(r, no.spp)

  #Empty list to store iterations -----
  out <- vector(mode = 'list', length = timesteps+1) #spp densities


  #Set initial data as first element of list
  out[[1]] <- initial_df



  #Neighborhood coordinates
  nh_coords <- initial_df%>%
    rowwise()%>%
    mutate(hood = purrr::map2(.x = x, .y = y,
                              .f = function(x,y)neighborhood(x = x, y = y, df = initial_df,
                                                             nh_size = nh_size,
                                                             neighborhood = nh,
                                                             torus = torus)))%>%
    select(`ID`, `x`, `y`, `hood`)%>%
    mutate(ID = as.character(ID))

  #return(nh_coords)

  #Iterate dynamics ----------
  for(i in 2:(timesteps+1)){

    #=================================
    #One iteration of the simulation -------
    #=================================

    # 1. Competition --------------------

    #apply bh_competition to each cell (single iteration)
    comp <- out[[i-1]]%>%
      select(sp_names)%>%
      apply(1, function(x)bh_competition(r = r, alpha = alpha,
                                         initial_N = x, discrete = TRUE,
                                         timesteps = 1, final_N = TRUE))%>%
      t()%>%
      as_tibble()%>%
      setNames(sp_names)%>%
      cbind(initial_df%>%select(`ID`, `x`,`y`), .)

    #return(comp)
    #out[[i]] <- comp

    # 2. Emigrate --------------------
    # Occurs after competition (fraction of remaining population become emigrants)
    #note: spp have same emigration rates
    emi <- comp%>%
      rowwise()%>%
      mutate_at(vars(sp_names), .funs = function(x)dd_emigration(N = x, dd = dd_emi, discrete = TRUE, max = disp_rate))%>%
      ungroup()%>%
      tidyr::pivot_longer(cols = sp_names, names_to = "Species", values_to = 'emigrants')%>%
      mutate(ID = as.character(ID))%>%
      filter(`emigrants` > 0)

    #return(emi)

    #3. Immigrate --------------------

    #3a. Calculate weights
    #Weights calculated according to the densities of previous timestep (i.e. delayed)
    weights <- nh_coords$hood%>%
      lapply(function(x)x%>%left_join(out[[i-1]]%>%select(-ID),
                                      by = c('x', 'y')))%>%
      lapply(function(x)x%>%
               mutate_at(sp_names, .funs = function(x)weighted_freq(N = x, delta = delta)))%>%
      bind_rows()%>%
      tidyr::pivot_longer(cols = sp_names, names_to = 'Species', values_to = 'weights')%>%
      group_by(`x`, `y`, `Species`)%>%
      group_nest(.key = 'weights')


    #return(weights)

    #3b. Distribute emigrants across neighborhoods
    immi <- emi%>%
      left_join(weights, by = c('x','y','Species'))%>%
      left_join(nh_coords, by = c('ID','x','y'))%>%
      mutate(immigration = purrr::pmap(.l = list(hood = hood, emigrants = emigrants, weights = weights),
                                        .f = function(hood, emigrants, weights)immigration(hood, emigrants, weights = weights$weights)))%>%
      select(`Species`, `immigration`)%>%
      unnest(immigration)%>%
      group_by(`Species`, `x`, `y`)%>%
      summarise(immigrants = sum(n), .groups = 'drop')



    # 4. Combine steps: growth - emigration + immigration --------------------
    out[[i]] <- comp%>%
      pivot_longer(cols = sp_names, names_to = 'Species', values_to = 'competition')%>%
      left_join(emi%>%select(-ID), by = c('x','y','Species'))%>%
      left_join(immi, by = c('x','y','Species'))%>%
      tidyr::replace_na(., list(emigrants = 0, immigrants = 0))%>%
      mutate(Density = `competition` - `emigrants` + `immigrants`)%>%
      pivot_wider(names_from = 'Species', values_from = 'Density', id_cols = c('ID','x','y'))

  }



  #Export data
  out <- out%>%
    bind_rows(.id = 'time')%>%
    mutate(time = as.integer(`time`))



  return(out)

}
