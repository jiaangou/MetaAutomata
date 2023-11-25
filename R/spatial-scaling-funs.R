library(dplyr)
library(ggplot2)
source('R/complete_functions.r')


co_sim <- readRDS("coexist_sim_20OCT23.rds")%>%
  bind_rows()

#Parameter space --------
par_range <- 20
short_range <- 3
replicates <- 5 #number of replicates for each parameter combination (5)
co_pars <- expand.grid(disp_rate = seq(from = 0, to = 0.4, length.out = par_range),
                       delta = seq(from = 0, to = 0.4, length.out = par_range),
                       aij = seq(from = 1, to = 2, length.out = 3),
                       p = seq(from = 0, to = 0.1, length.out = short_range),
                       c = seq(from = 0, to = 1, length.out = short_range),
                       weight_fit = c(FALSE, TRUE))%>%
  slice(rep(1:n(), each = replicates))



#Compute additional stats
coexist_dat <- co_pars%>%
  bind_cols(co_sim)

#Coexistence probability
coexist_dat%>%
  #filter(disp_rate < 0.3, delta < 0.3)%>%
  filter(disp_rate < 0.2, delta < 0.2)%>%
  mutate(coexist = ifelse(time == max(time), 1, 0))%>%
  group_by(disp_rate, delta, aij, p, c, weight_fit)%>%
  summarise(coexist_p = mean(coexist), coexist_time = mean(time))%>%
  filter(weight_fit == FALSE)%>%
  ggplot(aes(x = disp_rate, y = delta))+
  #geom_tile(aes(fill = coexist_p))+
  geom_tile(aes(fill = coexist_time))+
  scale_fill_gradient(low = 'grey', high = 'black')+
  facet_grid(p ~ aij + c)+
  theme_bw()




#Time-series
#space_ts <- readRDS('space_tibble_16OCT23.rds')
space_ts <- readRDS('tseries_22OCT23.rds')
test <- space_ts[[6]]



#Spatial scaling


#Step 1: Sample n random patches without replacement on the LxL grid
grid_sample <- function(n, L){
  check <- n < L^2 #can't draw more than the number of patches
  if(check == FALSE) stop("number of samples (n) cannot exceed total number of patches (L^2)")

  coords <- expand.grid(x = 1:L, y = 1:L)%>%
    sample_n(size = n, replace = FALSE)

  return(coords)
}
rand_patch <- grid_sample(n = 20, L = 50)

#Step 2: Get neighborhood coordinates
neighbor_offsets <- function(nh_size = 1, keep_origin = TRUE){

  offset <- expand.grid(nx = -nh_size:nh_size, ny =-nh_size:nh_size)%>%
    filter(abs(nx) + abs(ny) <= nh_size)

  if(keep_origin == FALSE){
    offset <- offset%>%
      filter(nx != 0 | ny != 0)

  }
  return(offset)
}
offset <- neighbor_offsets(1)

neighbors_coords <- lapply(1:nrow(rand_patch), function(x) offset%>%
         apply(1, function(y)y + as.numeric(rand_patch[x,]))%>%
         t())


#Step 3: get neighborhood states (state = sum species densities across patches)
lapply(neighbors_coords, function(x)x%>%
         apply(1, function(x)test[x[1], x[2],,50])%>%
         rowSums)%>%
  do.call(rbind, .)%>%
  apply(1, function(x)sum(x>0)) #species richness





























