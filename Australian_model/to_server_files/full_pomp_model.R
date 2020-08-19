#cureently set to model 8

# loading observed data
source("load_data.R", local = TRUE)
data <-load_data("CLU", "all")
samplesize_dat <- data[[1]]
pos_dat <- data[[2]]


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
init_states <- model_8_before_equ_ini

# creating the POMP object
pomp_object <- pomp(data=pos_dat,
                    times="time",t0=-365*12+1,
                    rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                    rinit=init_states,
                    skeleton=map(det_model_skeleton, delta.t = 1),
                    statenames=state_names,
                    paramnames=param_names,
                    partrans=parameter_trans(log=c("R0","c","s", "phi", "disp", "d", "gamma_val","omega_val", "omega_m_val", "s_v",
                                                   "phi_v", "rho_val", "epsilon_val"),
                                             logit = c("zeta")),
                    rmeasure=rmeas,
                    dmeasure=dmeas,
                    covar=covar_samplesize,
                    accumvars = "H",
                    cdir=".", cfile="hacking_win_bug"
)

