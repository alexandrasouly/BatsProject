#---------------------------------- MODEL 5 -------------------------------------------
#                      model 5 is SIR with no seasonal forcing 

# posterior parameters from Aaron's draft
model_5_params <- c(
  
  # base dynamics params
  beta_val = 0.06, # infection rate S -> R
  gamma_val = 3.135, # recovery rate I -> R
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
  c = 18.65, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.19, # birth pulse time shift
  
  # measuring process params
  zeta = 0.3, # test accuracy
  disp = 1000, # dispersion parameter
  d = 10 # number of bats contributing to the same pool
  
)

# initial states from Aaron's draft
 before_equ_ini <- Csnippet("
                        
                        Ma = 0;
                        Sn = 675;
                        Sj = 575;
                        Sm = 550;
                        Sf = 550;
                        
                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;
                        
                        In = 0;
                        Ij = 0;
                        Im = 27;
                        If = 27;
                        
                        Rn = 0;
                        Rj = 0;
                        Rm = 0;
                        Rf = 0;
                        
                        ")
 
 # 50 years waiting for equilibria 
 after_equ_ini <- Csnippet("
                        
                        Ma = 501;
                        Sn = 3;
                        Sj = 70;
                        Sm = 3;
                        Sf = 3;
                        
                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;
                        
                        In = 3;
                        Ij = 60;
                        Im = 7;
                        If = 7;
                        
                        Rn = 4;
                        Rj = 255;
                        Rm = 479;
                        Rf = 479;
                        
                        ")
   