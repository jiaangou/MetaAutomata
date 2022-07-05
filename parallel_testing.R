
#initial condition
L <- 3
init <- expand.grid(x = 1:L, y = 1:L, N1 = 10, N2 = 10)%>%
  tibble::rowid_to_column('ID')


#Parameters
pars <- expand.grid(disp_rate = seq(from = 0.1, to = 1, by = 0.1),
                    aij = seq(from = 1, to = 1.5, length.out = 3),
                    delta = seq(from = 0, to = 0.5, by = 0.1))

pars <- expand.grid(disp_rate = seq(from = 0.1, to = 1, length.out = 2),
                    #aij = seq(from = 1, to = 1.5, length.out = 2),
                    delta = seq(from = 0, to = 0.5, length.out  = 2))

library(furrr)


pars%>%
  mutate(sims = purrr::pmap(.l = list(disp_rate = disp_rate, delta = delta),
                            .f = function(disp_rate, delta)stochastic_sim(initial_df = init, aij = 1.1, disp_rate = disp_rate, delta = delta)))



#foreach parellel method
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
funs <- list(stochastic_sim)
system.time(
  foreach(i = 1:nrow(pars), .packages = 'MetaAutomata') %dopar% funs[[1]](initial_df = init, aij = 1.1, disp_rate = pars$disp_rate[i], delta = pars$delta[i])
)
stopCluster(cl)



out <- pars%>%
  mutate(sims = furrr::future_pmap(.l = list(disp_rate = disp_rate, delta = delta),
                   .f = function(disp_rate, delta)stochastic_sim(initial_df = init, aij = 1.1, disp_rate = disp_rate, delta = delta)))



mutate(immigration = purrr::pmap(.l = list(hood = hood, emigrants = emigrants, weights = weights),
                                 .f = function(hood, emigrants, weights)immigration(hood, emigrants, weights = weights$weights)))%>%
