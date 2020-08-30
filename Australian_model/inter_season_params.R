## model 8, 
#gamma and phi_v not fixed
# s and phi fixed
## k1, k2 

model_8_params_inter_seasons_lower1 <- c(
  
  # base dynamics params
  model_type = 2, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
  
  R0 = 22.46,
  gamma_val = 0.968, # recovery rate I -> R
  omega_val = 1.345, # immune waning rate R -> S
  omega_m_val = 0.801, # maternal antibody waning rate
  kappa_val = 5000, # carrying capacity
  rho_val = 0, # I -> L
  epsilon_val = 0, # L -> I
  mu_val = 1.37, # juvenile maturation rate
  mj_val = 0.501, # juvenile death rate
  m_val = 0.191, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  
  # seasonality params
  c = 12.586, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.179, # birth pulse time shift
  
  k1 = 0.1,
  k2 = 0.1,
  c_v1 = 1, # seasonal drive scaling factor
  s_v = 37.69, # seasonal drive synchronicity
  phi_v = 1.473, # sesonal drive time shift
  
  # measuring process params
  zeta = 0.046, # test accuracy
  disp = 10, # dispersion parameter
  d = 0.01 # number of bats contributing to the same pool
  
)

model_8_params_inter_seasons_upper1 <- c(
  
  # base dynamics params
  model_type = 2, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
  
  R0 = 49.849,
  gamma_val = 6.653, # recovery rate I -> R
  omega_val = 18.309, # immune waning rate R -> S
  omega_m_val = 0.801, # maternal antibody waning rate
  kappa_val = 5000, # carrying capacity
  rho_val = 0, # I -> L
  epsilon_val = 0, # L -> I
  mu_val = 1.37, # juvenile maturation rate
  mj_val = 0.501, # juvenile death rate
  m_val = 0.191, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  
  # seasonality params
  c = 19.817, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.179, # birth pulse time shift
  
  k1 = 3,
  k2 = 3,
  c_v1 = 1, # seasonal drive scaling factor
  s_v = 199.776, # seasonal drive synchronicity
  phi_v = 2.111 , # sesonal drive time shift
  
  # measuring process params
  zeta = 0.903, # test accuracy
  disp = 1000, # dispersion parameter
  d = 17.558 # number of bats contributing to the same pool
  
)

#----------------------------------------------------#
# fit cv1 as well

model_8_params_inter_seasons_lower2 <- c(
  
  # base dynamics params
  model_type = 2, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
  
  R0 = 22.46,
  gamma_val = 0.968, # recovery rate I -> R
  omega_val = 1.345, # immune waning rate R -> S
  omega_m_val = 0.801, # maternal antibody waning rate
  kappa_val = 5000, # carrying capacity
  rho_val = 0, # I -> L
  epsilon_val = 0, # L -> I
  mu_val = 1.37, # juvenile maturation rate
  mj_val = 0.501, # juvenile death rate
  m_val = 0.191, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  
  # seasonality params
  c = 12.586, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.179, # birth pulse time shift
  
  k1 = 0.1,
  k2 = 0.1,
  c_v1 = 0.1, # seasonal drive scaling factor
  s_v = 37.69, # seasonal drive synchronicity
  phi_v = 1.473, # sesonal drive time shift
  
  # measuring process params
  zeta = 0.046, # test accuracy
  disp = 10, # dispersion parameter
  d = 0.01 # number of bats contributing to the same pool
  
)

model_8_params_inter_seasons_upper2 <- c(
  
  # base dynamics params
  model_type = 2, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
  
  R0 = 49.849,
  gamma_val = 6.653, # recovery rate I -> R
  omega_val = 18.309, # immune waning rate R -> S
  omega_m_val = 0.801, # maternal antibody waning rate
  kappa_val = 5000, # carrying capacity
  rho_val = 0, # I -> L
  epsilon_val = 0, # L -> I
  mu_val = 1.37, # juvenile maturation rate
  mj_val = 0.501, # juvenile death rate
  m_val = 0.191, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  
  # seasonality params
  c = 19.817, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.179, # birth pulse time shift
  
  k1 = 3,
  k2 = 3,
  c_v1 = 3, # seasonal drive scaling factor
  s_v = 199.776, # seasonal drive synchronicity
  phi_v = 2.111 , # sesonal drive time shift
  
  # measuring process params
  zeta = 0.903, # test accuracy
  disp = 1000, # dispersion parameter
  d = 17.558 # number of bats contributing to the same pool
  
)