pomp_model_fn <-function(site, year, init_states, input_params){
  
# loading observed data
setwd("~/BatsProject/Australian_model")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pomp)

source("load_data.R", local = TRUE)
data <-load_data(site, year)
samplesize_dat <- data[[1]]
pos_dat <- data[[2]]

# pos_dat_fake0 <- rbind(c("time" = -365*50+1, "pos" = 0 ), pos_dat)
samplesize_dat <- rbind(c("time" = -365*12+1, "samplesize" = 0 ), samplesize_dat)
# feeding in test samplesizes to the model 
covar_samplesize <-covariate_table(samplesize_dat, 
                                   order = "constant", 
                                   times= "time")
# stochastic step model
source("stoch_step.R", local = TRUE)
stochStep <-stochStep()

# deterministic skeleton
source("det_model_skeleton.R", local = TRUE)
det_model_skeleton <- det_model_skeleton()

# measure processes, dmeas is likelihood, rmeas is to simulate
source("measure_processes.R", local = TRUE)
measure_processes <-measure_processes()
dmeas <- measure_processes[[1]]
rmeas <- measure_processes[[2]]

# initial state values and names of params and states
source("ini_and_params.R", local = TRUE)
state_names <- states()
param_names <- params()


source("aarons_params.R", local = TRUE)




# creating the POMP object
pomp_object <- pomp(data=pos_dat,
                    times="time",t0=-365*12+1,
                    rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                    rinit=init_states,
                    skeleton=map(det_model_skeleton, delta.t = 1),
                    statenames=state_names,
                    paramnames=param_names,
                    partrans=parameter_trans(log=c("R0","c","s", "phi", "disp", "d", "gamma_val", "omega_m_val", "s_v",
                                                   "phi_v", "rho_val", "epsilon_val"),
                                             logit = c("zeta")),
                    rmeasure=rmeas,
                    dmeasure=dmeas,
                    covar=covar_samplesize,
                    accumvars = "H",
                    cdir=".", cfile="hacking_win_bug"
)

# pop_equ_pomp_model <- pomp(data=data.frame("time" = seq(1, 365*50)),
#                               times="time",t0=0,
#                               rprocess=discrete_time(step.fun=stochStep,delta.t=1),
#                               rinit=model_4_before_equ_ini,
#                               skeleton=map(det_model_skeleton, delta.t = 1),
#                               statenames=state_names,
#                               paramnames=param_names
#  )
# 
#  pop_equ <- trajectory(pop_equ_pomp_model, params=model_4_params,format="d")
# tail(pop_equ, 1)
# plot()

#simulating and calculating deterministic trajectory
sim <- simulate(pomp_object,params=input_params,format = "data.frame", nsim = 100, time = seq(0, 356*3))
x <- trajectory(pomp_object,params=input_params,format="d", time = seq(0, 356*3))

# 
colnames(pos_dat) <-c("time", "true_pos")
sim_plus_data <- merge(sim, pos_dat, by.sim = "time")
sim_plus_data <- merge(sim_plus_data, samplesize_dat, by.sim_plus_data = "time")
sim_plus_data %>% mutate(true_test_prev = true_pos/samplesize, sim_test_prev = pos/samplesize,
                         sim_model_prev = (In + Ij + Im + If)/
                           (Sn + Sj + Sm + Sf + En + Ej +
                              Em + Ef + In + Ij + Im + If +
                              Rn + Rj + Rm + Rf + Ma)
) -> sim_plus_data

# 
# 
source("plots.R")
plots <- plots(sim_plus_data, x, sim)
names(plots) <- c(   "plIm", "plRm", "plEm", "plSm",
                     "plIj", "plRj", "plEj", "plSj",
                     "plIn", "plRn", "plEn", "plSn",
                     "plMa", "plH", "plprev")

# 

return(plots)
}