# loading observed data
source("load_data.R")
data <-load_data(site = "CLU", year = "2019")
samplesize_dat <- data[[1]]
pos_dat <- data[[2]]

# feeding in test samplesizes to the model 
covar_samplesize <-covariate_table(samplesize_dat, 
                                   order = "constant", 
                                   times= "time")
# stochastic step model
source("stoch_step.R")
stochStep <-stochStep()

# deterministic skeleton
source("det_model_skeleton.R")
det_model_skeleton <- det_model_skeleton()

# measure processes, dmeas is likelihood, rmeas is to simulate
source("measure_processes.R")
measure_processes <-measure_processes()
dmeas <- measure_processes[[1]]
rmeas <- measure_processes[[2]]

# initial state values and names of params and states
source("ini_and_params.R")
state_names <- states()
param_names <- params()
init_states <- inits()

# specifying parameter values
params <- c(
  
  # base dynamics params
  beta_val = 0.00416, # infection rate S -> R
  gamma_val = 3.13, # recovery rate I -> R
  omega_val = 0, # immune waning rate R -> S
  omega_m_val = 0.8, # maternal antibody waning rate
  kappa_val = 2350, # carrying capacity
  rho_val = 0, # I -> L 
  epsilon_val = 0, # L -> I 
  mu_val = 0.44, # juvenile maturation rate
  mj_val = 0.50, # juvenile death rate
  m_val = 0.23, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  
  # seasonality params
  c = 18, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.9, # birth pulse time shift
  
  # measuring process params
  zeta = 0.3, # test accuracy
  disp = 1000, # dispersion parameter
  d = 10 # number of bats contributing to the same pool
)

# creating the POMP object
pomp_object <- pomp(data=pos_dat,
                    times="time",t0=734,
                    rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                    rinit=init_states,
                    skeleton=map(det_model_skeleton, delta.t = 1),
                    statenames=state_names,
                    paramnames=param_names,
                    rmeasure=rmeas,
                    dmeasure=dmeas,
                    covar=covar_samplesize
)

# simulating and calculating deterministic trajectory
sim <- simulate(pomp_object,params=params,format = "data.frame", nsim = 100)
x <- trajectory(pomp_object,params=params,format="d")


colnames(pos_dat) <-c("time", "true_pos")
sim_plus_data <- merge(sim, pos_dat, by.sim = "time")
sim_plus_data <- merge(sim_plus_data, samplesize_dat, by.sim_plus_data = "time")
sim_plus_data %>% mutate(true_test_prev = true_pos/samplesize, sim_test_prev = pos/samplesize,
                         sim_model_prev = (In + Ij + Im + If)/
                           (Sn + Sj + Sm + Sf + En + Ej + 
                              Em + Ef + In + Ij + Im + If +
                              Rn + Rj + Rm + Rf + Ma)
) -> sim_plus_data



source("plots.R")
plots <- plots(sim_plus_data, x, sim)
names(plots) <- c("plprev", 
                     "plIm", "plRm", "plSm", 
                     "plIj", "plRj", "plSj", 
                     "plIn", "plRn", "plSn",
                     "plMa")




