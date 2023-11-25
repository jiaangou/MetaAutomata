#FUNCTIONS ---------------
library(dplyr)

########################
#Auxillary fucntions -----
########################
#Dominance matrix: takes L*L*Spp array as input and outputs a L*L matrix with species dominance (1/0)
dominance_matrix <- function(N_array){

  N_array%>%
    apply(c(1,2), function(x)(x[1]>x[2])%>%as.integer)
}

#Coarse graining Spp densities into a single state variable: 11 both species present, 10 - N1 present, 01 - N2 present, 00 - both extinct
coarse_graining <- function(N_array, states = c(4,3)){

  if(states == 4){
    out <- N_array%>%
      apply(c(1,2), function(x)(x>0)%>%
              as.integer%>%
              paste(collapse = ''))

  }else if(states == 3){
   out <- N_array%>%
      apply(c(1,2), function(x)
      ifelse(sum(x) == 0, 0, which.max(x)))

  }else(stop("select either 3 or 4 for the number of coarse grained states"))


  return(out)
}


#Neighborhood array: takes L*L dominance matrix as input and outputs an array of each spatial position and its 4 adjacent neighbors (L*L*4), if rearrange is true, it takes the array and rotates each slice back to its original orientation
neighborhood_array <- function(matrix, L = NULL, rearrange = FALSE){

  #Set up rotations
  if(is.null(L)){
    L <- dim(matrix)[1]
  }

  minus_one <- c(tail(1:L, 1), head(1:L, -1))
  plus_one <- c(tail(1:L, -1), head(1:L, 1))

  #Rotate matrices ----
  #If rearrange is FALSE, then output neighborhood matrices
  if(rearrange == FALSE){

    out <- array(c(matrix[,minus_one],
                   matrix[,plus_one],
                   matrix[minus_one,],
                   matrix[plus_one,]), dim = c(L, L, 4))
    #return(out)
  }
  #If rearrange is TRUE, rearrange each neighborhood matrix back to its original x y positions
  else{



    out <- array(c(matrix[,plus_one, 1],
                   matrix[,minus_one, 2],
                   matrix[plus_one,,3],
                   matrix[minus_one,,4]), dim = c(L, L, 4))


  }

  return(out)
}


#Sample patches
sample_neighbors <- function(patches = NULL, num_samples, weights) {

  if(is.null(patches)){
    patches <- 1:length(weights)
  }
  #Initiate empty vector
  out <- setNames(rep(0, length(patches)), patches)

  #Sample
  samples <- sample(patches, size = num_samples, replace = TRUE, prob = weights)
  count_table <- table(samples)

  #Fill vector with those that have >0 counts
  out[names(count_table)] <- count_table

  return(out)
}

#Neighborhood entropy (configuration entropy): takes L*L*4 neighborhood array as input and outputs associated entropy or frequency table (if table = TRUE)
neighborhood_entropy <- function(array, table = FALSE){
  out <- apply(array, 3, rbind)%>%
    as.data.frame%>%
    table()%>%
    as.data.frame()
  if(table == FALSE){
    out <- out%>%
      mutate(p = Freq/sum(Freq))%>%
      filter(p > 0)%>%
      summarise(entropy = -sum(p*log(p)))
  }
  return(out)


}


#Nearest neighbor arrangement long form
neighbor_long <- function(matrix){

  cell_states <- matrix%>%
    apply(c(1,2), function(x)(x>0)%>%
            as.integer%>%
            paste(collapse = ''))%>%
    reshape::melt()%>%
    rename('y' = X1, 'x' = X2, 'cell_state' = value)



  neighbor_states <- matrix%>%
    apply(c(1,2), function(x)(x>0)%>%
            as.integer%>%
            paste(collapse = ''))%>%
    neighborhood_array()%>%
    reshape::melt()%>%
    rename('y' = X1, 'x' = X2, 'neighbor' = X3, 'n_state' = value)

  out <- cell_states%>%
    left_join(neighbor_states, by = c('x','y'))


  return(out)

}

#Functions for extracting neighborhood information (their coordinates given size, their values, )
#NOTE: Since coordinates are used multiple times to extract values, it is more efficient to generate the whole neighborhood cooridnates once and use it repeatedely. Thus, I have written two separate functions, one for generating the nh coords and the other for extracting values

#Set up array with neighborhood coordinates for each xy coordinate
neighbor_coords <- function(L, nh_size = 1, vonNeumann = TRUE, no_dimensions = 2){

  #setup neighbor offsets
  offsets <- expand.grid(nx = -nh_size:nh_size, ny = -nh_size:nh_size)%>%
    filter(nx != 0 | ny != 0) #removes origin

  if(vonNeumann == TRUE){
    offsets <- offsets%>%
      filter(abs(nx) + abs(ny) <= nh_size) #slice off corners
  }

  #Setup output array: 4D:  neighbor dimension (x or y), focal row, focal column, n_neighbors
  n_neighbors <- nrow(offsets)
  nh_coord_array <- array(NA, dim = c(no_dimensions, L, L, n_neighbors))

  #Loop through rows and columns of matrix
  for(r in 1:L){
    for(c in 1:L){

      nh_coord_array[,r,c,] <- (offsets%>%t() + c(c,r))%>% # converst offsets to neighbor coordinates
        t()%>%
        apply(1, function(x)ifelse(x < 1, L + x, x))%>% #corrects boundaries so they are periodic
        apply(2, function(x)ifelse(x > L, x - L, x))

      #print(neighbor_mat)
      # nh_coord_array[,r,c,] <- neighbor_mat%>%
      #   apply(1, function(x)matrix[x[2], x[1]])
    }

  }
  return(nh_coord_array)
}

#Remove redundant neighbor pairs since edges are undirected (NOTE: doesn't effect calculations but reduces memory)
undir_pair <- function(coords){

  #melt to long form
  out <- coords%>%
    reshape2::melt()%>%
    rename(dimension = `Var1`)%>%
    mutate(dimension = plyr::mapvalues(dimension, from = c(1,2), to = c('nx','ny')))%>%
    rename(fy = `Var2`, fx = `Var3`, neighbor_id = `Var4`)%>%
    tidyr::pivot_wider(names_from = dimension, values_from = value)%>%
    mutate(f_xy = paste0(fx, fy))%>%
    mutate(n_xy = paste0(nx, ny))%>%
    distinct(smaller = pmin(f_xy, n_xy),
             larger = pmax(f_xy, n_xy),
             .keep_all = TRUE)%>%
    select(fy, fx, ny, nx)

  return(out)
}

#Take neighborhood cooridnates and state values as input and return neighborhood states
nh_vals <- function(nh_coord, vals){

  r_L <- dim(nh_coord)[2]
  c_L <- dim(nh_coord)[3]
  n_neighbors  <- dim(nh_coord)[4]

  #Dimensions: row length, column length, number of neighbors
  nh_vals <- array(NA, dim = c(r_L,
                               c_L,
                               n_neighbors))

  for(r in 1:r_L){
    for(c in 1:c_L){

      nh_vals[r,c,] <- nh_coord[,r,c,]%>%
        apply(2, function(x)vals[x[1], x[2]])


    }
  }

  return(nh_vals)
}



#Take N array as input and return entropies: hx, hy, hxy, mi
array_entropies <- function(array){
  require(entropy)
  h <- array%>%
    apply(3, function(x)entropy(x, unit = 'log2'))

  hxy <- array%>%
    apply(c(1,2), sum)%>%
    entropy(unit = 'log2')

  mi <- sum(h) - hxy

  out <-c(h, hxy, mi)
  names(out) <- c('hx','hy','hxy','mi')
  return(out)

}


#Total competitive effects
# total_competition <- function(N, alpha){
#   out <- N%>%
#     apply(c(1,2), function(x)x%*%alpha)%>%
#     aperm(c(2,3,1))
#   return(out)
# }


#Compute metrics
#Metrics: Regional Ns, Hx, MI
compute_metrics <- function(N, n_coords){

  #Regional frequencies
  regional_N <- apply(N, 3, sum)

  #Coarse-graning
  cg <- N%>%
    coarse_graining(states = 3)

  cg4 <- N%>%
    coarse_graining(states = 4)%>%
    factor(levels = c('00', '10', '01', '11'))%>%
    table%>%
    as.vector()

  #Entropy
  hx <- cg%>%
    table%>%
    infotheo::entropy()

  #mutual info
  mi <- n_coords%>%
    rowwise()%>%
    mutate(f_state = cg[fy,fx])%>%
    mutate(n_state = cg[ny, nx])%>%
    ungroup%>%
    summarise(I = infotheo::mutinformation(f_state, n_state))%>%
    as.numeric

  out <- c(regional_N, hx, mi, cg4)

  return(out)

}


#Iterator ------
immi_iterator <- function(w, surv){

  #weights is a 3d array:
  #survivors is 2d array

  #Output has same dimensions as weights
  out <- array(0, dim = dim(w))
  rows <- dim(surv)[1]
  cols <- dim(surv)[2]

  for(r in 1:rows){
    for(c in 1:cols){

      #sample_neighbors(num_samlpes = survivors[r,c], weights = weights[r,c,]
      out[r,c,] <- sample_neighbors(num_samples = surv[r,c], weights = w[r,c,])

    }
  }

  return(out)

}



########################
# Utiliization functions ----
########################
#Niche overlap measures -----
overlap <- function(P, Q){
  kl <- sum(P*log2(P/Q))
  niche_overlap <- sum(P*Q)/sum(P^2)
  return(c(kl = kl, niche_overlap = niche_overlap))

}

KL <- function(P, Q, mu = TRUE){
  pq <- sum(P*log2(P/Q))
  qp <- sum(Q*log2(Q/P))

  if(mu == TRUE){
    out <- mean(c(pq, qp))
  }else{
    out <- pq
  }

  return(out)
}

niche_overlap <- function(P, Q, mu = TRUE){
  pq <- sum(P*Q)/sum(P^2)
  qp <- sum(P*Q)/sum(Q^2)

  if(mu == TRUE){
    out <- mean(c(pq, qp))
  }else{
    out <- pq
  }

  return(out)
}




####################################
#NEW UPDATED: Compute stats -------
###################################
compute_stats <- function(Nt, Nt_minus1, I){

  #Landscape size
  L <- dim(Nt)[1]

  #Regional densities---------------
  reg_N <- Nt%>%
    apply(3, sum)

  #Information-theoretic metrics -----------
  #3-state coarse-graing (i.e. dominance or exclusion)
  cg3 <- Nt%>%
    coarse_graining(states = 3)

  entropies <- cg3%>%
    neighborhood_array()%>% #Align nearest neighbors: left, bottom, right, top
    apply(3, function(x)x)%>%
    matrix(ncol = 1)%>% #collapse to a single column
    cbind(., matrix(rep(cg3,4), ncol = 1))%>% #align with focal cell by duplicating 4 times and collapsing
    as.data.frame()%>%
    summarise(I = infotheo::mutinformation(V1, V2), Hx = infotheo::entropy(V1))%>%
    mutate_all(infotheo::natstobits)%>%
    as.numeric()

  #Co-occurrence -------------
  co_df <- Nt%>%
    apply(3, function(x)x)%>%
    as.data.frame()%>%
    filter(V1>0 & V2>0)

  co_p <- nrow(co_df)/L^2

  #Niche overlap (MacArthur & Levins 1967) (NOTE: overlap is set to 0 if no co-occurrence)
  if(co_p > 0){
    overlap <- co_df%>%
      mutate_all(.funs = function(x)x/sum(x))%>% #normalize densities to probabilities (sum to 1)
      summarise(overlap = niche_overlap(P = V1, Q = V2, mu = FALSE))%>% #compute niche overlap
      as.numeric
  }else{
    overlap <- 0
  }


  #Spatial stability ---------
  #N <- array(c(Nt_minus1, Nt), dim = c(L, L, 2, 2)) #combine the two time points into an array
  #Coarse-grain states at T-1
  cg3_tminus1 <- Nt_minus1%>%
    coarse_graining(states = 3)

  #Compute probability of "bit flipping"  (1 - # of unchanged states)
  spat_stability <- 1 - (sum(cg3 == cg3_tminus1)/L^2)


  #Immigration assortativity -------
  I_long <- I%>%
    apply(3, function(x)x)%>%
    as.data.frame()

  #Tally state frequencies, normalize to 1, compute KL divergence
  immi_kl <- I_long%>%
    mutate(state = as.factor(cg3))%>%
    tidyr::pivot_longer(cols = 1:2, names_to = 'Species', values_to = 'count')%>%
    group_by(state, Species)%>%
    summarise(count = sum(count), .groups = 'drop')%>%
    tidyr::pivot_wider(names_from = 'Species', values_from = 'count')%>%
    filter(V1 >0 & V2>0)%>%
    mutate(P = V1/sum(V1))%>%
    mutate(Q = V2/sum(V2))%>%
    summarise(KL = KL(P = P, Q = Q, mu = FALSE))%>%
    .$KL


  #Compile
  stats <- c(reg_N, entropies, co_p, overlap, spat_stability, immi_kl)
  names(stats) <- c('r_N1', 'r_N2', 'I', 'Hx', 'co_p', 'overlap', 'spatial_stability', 'Immi_KL')

  return(stats)



}


########################
# Fitness distirbution of immigrants ----
########################

fitness_dist <- function(N_array, I_array, r, A, mu = TRUE){

  f_vals <-  N_array%>%
    apply(c(1,2), function(x)r/(1 + x%*%A))%>%
    apply(1, function(x)x)%>%
    reshape2::melt()%>%
    setNames(c('Patch', 'Spp', 'lambda'))

  I_vals <- I_array%>%
    apply(3, function(x)x)%>%
    reshape2::melt()%>%
    setNames(c('Patch', 'Spp', 'I'))


  lambda <- f_vals%>%
    left_join(I_vals, by = c('Patch','Spp'))


  #If mu is true, return the average lambda. If false, return the distribution of fitness values
  if(mu == TRUE){
    out <- lambda%>%
      group_by(Spp)%>%
      summarise(lambda_mu = weighted.mean(x = lambda, w = I))%>%
      .$lambda_mu

    names(out) <- c('lambda_I_Sp1', 'lambda_I_Sp2')

  }else{
    out <- f_vals%>%
      left_join(I_vals, by = c('Patch','Spp'))%>%
      tidyr::uncount(I)
  }


  return(out)

}



########################
#Spatial stability ----------
########################
spatial_stability <- function(N_array){

  #Coarse-grain
  spatial_cg <- N_array%>%
    apply(4, function(x)x%>%
            coarse_graining(states = 3), simplify = FALSE)
  #Timesteps
  len <- length(spatial_cg)

  #Tally the number of state changes for each patch
  state_change <- purrr::map2(.x = spatial_cg[-t_len], .y = spatial_cg[-1],
                              .f = function(x,y)sum(x == y))%>%
    unlist()

  #return
  return(state_change)

}



########################
#Model functions ----
########################


#1. Competition ------------

competition <- function(N, alpha, r){
  no.spp <- length(r)

  #Competitive effects
  comp <- apply(N, c(1,2), function(x)x%*%alpha)%>%
    aperm(c(2,3,1))

  #Update densities
  N_new <- r*N/(1 + comp)
  #Add stochasticity
  N_out <- apply(N_new, c(1,2), function(x)rpois(no.spp, x))%>%
    aperm(c(2,3,1))

  return(N_out)
}

#2. Emigration ------------
emigration <- function(N, disp){
  out <- apply(N, c(1,2), function(x)rbinom(length(x), size = x, prob = disp))%>%
    aperm(c(2,3,1))
  return(out)
}

#3. Viability  ------------
viability <- function(N, cost){
  out <- apply(N, c(1,2), function(x)rbinom(length(x), size = x, prob = 1 - cost))%>%
    aperm(c(2,3,1))
  return(out)
}

#4. Immigration   (dependencies: neighborhood_array, weighted_freq, immi_iterator) ------------
immi_function <- function(N_array, survivors, delta, no.spp = 2){

  #Get number of species, equal to length of 3rd dimension (not really useful)
  #no.spp <- dim(N_array)[3]

  #Get the neigborhood arrays of each species and their weights: each list element is a LxLx4 array
  spp_nh_array <- lapply(1:no.spp, function(x)neighborhood_array(N_array[,,x]))


  #Calculate the weights of each neighbor
  weights <- lapply(1:no.spp, function(x)neighborhood_array(N_array[,,x]))%>%
    lapply(function(x)x%>%
             apply(c(1,2), function(x)MetaAutomata::weighted_freq(x, delta = delta))%>%
             aperm(c(2,3,1)))

  #Sample each neighbor with prob = weights and then rearrange to stack the same xy cells together
  out <- lapply(1:no.spp, function(s)
    immi_iterator(w = weights[[s]],
                  surv = survivors[,,s]))%>%
    lapply(function(x)x%>%
             neighborhood_array(rearrange = TRUE))%>%
    lapply(function(x)x%>%apply(c(1,2), sum))


  return(out)
}


#5. Invasion fitness: growth rate of an individual if it were to invade a patch
invasion_fitness <- function(N, r, A){
  #N_1 <- N + 1
  #W <- (r*(N_1) / (1 + N_1%*%A)) / N_1
  W <- r / (1 + N%*%A)
  return(W)

}


########################
#Simulator 1 ----
########################
########################
# simulator <- function(N, r = 1.5, aii = 1, aij, K = 50,
#                       d_disp, delta, cost, ext_p = 0,
#                       time = 1, n_coords, cg_states = 3,
#                       metrics = FALSE, record_interval = 1, save_intervals = FALSE){
#
#   #Set up parameters: spatial extent (L), number of species, alpha matrix (A)
#   L <- dim(N)[1]
#   no.spp <- length(r)
#   A =  matrix(aij, nrow = no.spp, ncol = no.spp);diag(A) = aii;A <- A/(K/(r-1)) #scale by carrying capacity
#
#   #Counters -----
#   n_stored <- time%/%record_interval
#   j <- 1 #for counting the number of recorded steps
#
#   #Set up variables to store information
#   if(save_intervals == TRUE){
#     if(metrics == TRUE){
#       metrics_out <- matrix(NA, nrow = n_stored, ncol = no.spp + 9) #regional densities, Hx, MI, timestep
#       colnames(metrics_out) <- c(paste0('r_N', 1:no.spp), '00', '10','01','11', 'I','Hx', 'KL', 'niche_overlap', 'time')
#
#     }else{
#       metrics_out <- matrix(NA, nrow = n_stored, ncol = no.spp + 9) #regional densities, Hx, MI, timestep
#       colnames(metrics_out) <- c(paste0('r_N', 1:no.spp), '00', '10','01','11', 'I','Hx', 'KL', 'niche_overlap', 'time')
#
#       N_t <- array(0, dim = c(L, L, no.spp, n_stored))
#     }
#
#   }
#
#
#
#   ##############################
#   #Iterate dynamics ---------
#   for(i in 1:time){
#
#     #############
#     #Competition
#     N_prime <- N%>%
#       competition(alpha = A, r = r)
#
#     #Emigration
#     E <- N_prime%>%
#       emigration(disp = d_disp)
#
#     #Survival
#     S <- E%>%
#       viability(cost = cost)
#
#     #Immigration
#     I <- immi_function(N_prime, survivors = S, delta = delta)%>%
#       unlist%>%
#       array(., dim = c(L, L, no.spp))
#
#     #Disturbance
#     D_mat <- matrix(rbinom(L^2, 1, prob = 1-ext_p), nrow = L, ncol = L)%>%
#       rep(no.spp)%>%
#       array(., dim = c(L, L, no.spp))
#
#     #Update states
#     N <- (N_prime - E + I) * D_mat
#
#     #Check for extinction status by computing regional Ns and checking if any is 0
#     extinction <- any(apply(N, 3, sum) == 0) #if any N is 0, then extinction == TRUE
#     #############
#
#
#     ######################################
#     #Fitness distribution of immigrants -----------
#     #f_dist <- fitness_dist(N_array = N_prime + I, I_array = I, r = r, A = A)
#     #return(f_dist)
#     ######################################
#
#
#     #############
#     #Compute metrics
#     if(save_intervals == TRUE){
#
#       #Check if timestep is divisible by
#       if(i %% record_interval == 0){
#         if(metrics == TRUE){
#           #Compute only if no extinction has occurred
#           if(extinction == FALSE){
#             #Compute and reocord
#             metrics_out[j,] <- c(compute_stats(N), i)
#             j <- j+1 #update counter
#           }
#         }else{
#           #Record densities
#           N_t[,,,j] <- N
#           j <- j+1 #update counter
#         }
#       }
#     }
#     #############
#
#
#     #Stop iterations if one species goes extinct
#     if(extinction){
#       if(save_intervals == TRUE){
#         #Trim output arrays (remove NAs)
#         if(metrics == TRUE){
#           #metrics_out[j,] <- c(compute_stats(N), i) #Compute extinction step
#           #metrics_out[j,] <- c(compute_metrics(N = N, n_coords = n_coords, cg_states = cg_states), i) #Compute extinction step
#           metrics_out <- metrics_out[1:j,] #remove Nas
#         }else{
#           N_t[,,,j] <- N
#           N_t <- N_t[,,,1:j] #
#         }
#       }
#       break
#     }
#   }
#   ##############################
#
#   if(save_intervals == FALSE){
#     #Compute last time step
#     metrics_out <- c(compute_stats(N), i) #Compute extinction step
#     #metrics_out <- c(compute_metrics(N = N, n_coords = n_coords, cg_states = cg_states), i)
#     names(metrics_out) <- c(paste0('r_N', 1:no.spp), '00', '10','01','11', 'I','Hx', 'KL', 'niche_overlap', 'time')
#     #names(metrics_out) <- c(paste0('r_N', 1:no.spp), 'Hx', 'MI', 'time')
#     #Save last state
#     N_t <- N
#   }
#
#   #Return outputs
#   if(metrics == TRUE){
#     return(metrics_out)
#
#   }else{
#     return(N_t)
#   }
#
# }

########################
#Simulator 2 ----
########################
########################

simulator <- function(N, r = 1.5, aii = 1, aij, K = 50,
                      d_disp, delta, cost_c, ext_p = 0,
                      weight_fitness = FALSE, time = 1,
                      metrics = FALSE, record_interval = 1, save_intervals = FALSE){
  ##############################
  #Setup -------
  ##############################
  #Spatial extent
  L <- dim(N)[1]
  #Number of species
  no.spp <- length(r)
  #Alpha matrix
  A =  matrix(aij, nrow = no.spp, ncol = no.spp);diag(A) = aii;A <- A/(K/(r-1)) #scale by carrying capacity
  #Variables to store iterations
  n_stored <- time%/%record_interval #number of stored is total simulation time modolus intervals to stores
  j <- 1 #for counting the number of recorded steps


  #Set up arrays for recording states
  #Stats to compute: regional densities (r_N1, r_N2), co-occurrence prob (11 / L^2), mutual information
  if(save_intervals == TRUE){
    if(metrics == TRUE){
      metrics_out <- matrix(NA, nrow = n_stored, ncol = no.spp + 7)
      colnames(metrics_out) <- c(paste0('r_N', 1:no.spp), 'I','Hx', 'cooccurrence', 'niche_overlap', 'spatial_stability', 'immi_kl', 'time')

    }else{
      N_t <- array(0, dim = c(L, L, no.spp, n_stored))
    }

  }

  ##############################
  #Iterate dynamics ---------
  for(i in 1:time){

    #Store initial condition (before dynamics: competition, emigration, viability, immigration, extinction)
    #NOTE: this is needed for computing statistics
    N0 <- N

    #############
    #Competition
    N_prime <- N%>%
      competition(alpha = A, r = r)

    #Emigration
    E <- N_prime%>%
      emigration(disp = d_disp)

    #Survival
    cost <- 1 - exp(-cost_c*delta) #calculate dispersal cost ( = 1 - survival prob) given delta and c
    S <- E%>%
      viability(cost = cost)


    #Invasion fitness
    if(weight_fitness == TRUE){
      W <- N_prime%>%
        apply(c(1,2), function(x)invasion_fitness(x, r = r, A = A)-1)%>% #minus 1 offsets the plus 1 constant when weighting by N
        aperm(c(2,3,1))

    }else{
      W <- N_prime
    }


    #Immigration
    I <- immi_function(W, survivors = S, delta = delta)%>%
      unlist%>%
      array(., dim = c(L, L, no.spp))

    #Disturbance
    D_mat <- matrix(rbinom(L^2, 1, prob = 1-ext_p), nrow = L, ncol = L)%>%
      rep(no.spp)%>%
      array(., dim = c(L, L, no.spp))

    #Update states
    N <- (N_prime - E + I) * D_mat

    #Check for extinction status by computing regional Ns and checking if any is 0
    extinction <- any(apply(N, 3, sum) == 0) #if any N is 0, then extinction == TRUE
    #############


    ####################################
    ####################################
    #print(compute_stats(Nt = N, Nt_minus1 = N0, I = I))
    ####################################
    ####################################


    #############
    #Compute metrics
    if(save_intervals == TRUE){

      #Check if timestep is divisible by
      if(i %% record_interval == 0){
        if(metrics == TRUE){
          #Compute only if no extinction has occurred
          if(extinction == FALSE){
            #Compute and reocord
            metrics_out[j,] <- c(compute_stats(Nt = N, Nt_minus1 = N0, I = I), i)
            j <- j+1 #update counter
          }
        }else{
          #Record densities
          N_t[,,,j] <- N
          j <- j+1 #update counter
        }
      }
    }
    #############


    #Stop iterations if one species goes extinct
    if(extinction){
      if(save_intervals == TRUE){
        #Trim output arrays (remove NAs)
        if(metrics == TRUE){
          metrics_out <- metrics_out[1:j,] #remove Nas
        }else{
          N_t[,,,j] <- N
          N_t <- N_t[,,,1:j] #
        }
      }
      break
    }
  }
  ##############################

  if(save_intervals == FALSE){
    #Compute last time step
    metrics_out <- c(compute_stats(Nt = N, Nt_minus1 = N0, I = I), i) #Compute only last time step
    names(metrics_out) <- c(paste0('r_N', 1:no.spp), 'I','Hx', 'cooccurrence', 'niche_overlap', 'spatial_stability', 'immi_kl', 'time')
    #Save last state
    N_t <- N
  }

  #Return outputs
  if(metrics == TRUE){
    return(metrics_out)

  }else{
    return(N_t)
  }

}


