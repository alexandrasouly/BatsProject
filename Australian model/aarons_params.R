



#---------------------------------- MODEL 5 -------------------------------------------
#                      model 5 is SIR with no seasonal forcing 

# posterior parameters from Aaron's draft
model_5_params <- c(
  
  # base dynamics params
  model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only 
  
  beta_val = 0.06, # infection rate S -> R
  gamma_val = 3.135, # recovery rate I -> R
  omega_val = 0, # immune waning rate R -> S
  omega_m_val = 0.8, # maternal antibody waning rate
  kappa_val = 2350, # carrying capacity
  rho_val = 0, # I -> L 
  epsilon_val = 0, # L -> I
  

  mu_val = 1.37, # juvenile maturation rate
  mj_val = 0.50, # juvenile death rate
  m_val = 0.23, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  
  # seasonality params
  c = 18.65, # birth pulse scaling factor
  s = 130, # birth pulse synchronicity
  phi = 7.19, # birth pulse time shift
  
  c_v = 1, # seasonal drive scaling factor
  s_v = 0, # seasonal drive synchronicity
  phi_v = 0, # sesonal drive time shift
  
  # measuring process params
  zeta = 0.3, # test accuracy
  disp = 1000, # dispersion parameter
  d = 10 # number of bats contributing to the same pool
  
)



# initial states from Aaron's draft
 model_5_before_equ_ini <- Csnippet("
                        
                        Ma = 0;
                        Sn = 591;
                        Sj = 253;
                        Sm = 753;
                        Sf = 753;
                        
                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;
                        
                        In = 0;
                        Ij = 0;
                        Im = 29;
                        If = 29;
                        
                        Rn = 0;
                        Rj = 0;
                        Rm = 0;
                        Rf = 0;
                        
                        ")
 
 # 50 years waiting for equilibria 
 model_5_after_equ_ini <- Csnippet("
                        
                        Ma = 685;
                        Sn = 5;
                        Sj = 68;
                        Sm = 6;
                        Sf = 6;
                        
                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;
                        
                        In = 8;
                        Ij = 65;
                        Im = 20;
                        If = 20;
                        
                        Rn = 10;
                        Rj = 107;
                        Rm = 696;
                        Rf = 696;
                        
                        ")
 
 
 #---------------------------------- MODEL 6 -------------------------------------------
 #                      model 6 is SIR with seasonal forcing 
 
 # posterior parameters from Aaron's draft
 model_6_params <- c(
   
   # base dynamics params
   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only 
   
   beta_val = 0.0719187, # infection rate S -> R
   gamma_val = 22.449, # recovery rate I -> R
   omega_val = 0, # immune waning rate R -> S
   omega_m_val = 0.801, # maternal antibody waning rate
   kappa_val = 4216, # carrying capacity
   rho_val = 0, # I -> L 
   epsilon_val = 0, # L -> I 
   mu_val = 1.37, # juvenile maturation rate
   mj_val = 0.499, # juvenile death rate
   m_val = 0.187, # adult death rate
   delta_t = 365, # scaling time as days instead of years
   
   # seasonality params
   c = 15.973, # birth pulse scaling factor
   s = 129.841, # birth pulse synchronicity
   phi = 7.181, # birth pulse time shift
   
   c_v = 1, # seasonal drive scaling factor
   s_v = 26.581, # seasonal drive synchronicity
   phi_v = 0.573, # sesonal drive time shift
   
   # measuring process params
   zeta = 0.433, # test accuracy
   disp = 1000, # dispersion parameter
   d = 8.2 # number of bats contributing to the same pool
   
 )
 
 # initial states from Aaron's draft
 model_6_before_equ_ini <- Csnippet("
                        
                        Ma = 0;
                        Sn = 922;
                        Sj = 396;
                        Sm = 1449;
                        Sf = 1449;
                        
                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;
                        
                        In = 0;
                        Ij = 0;
                        Im = 105;
                        If = 105;
                        
                        Rn = 0;
                        Rj = 0;
                        Rm = 0;
                        Rf = 0;
                        
                        ")
 
 # 50 years waiting for equilibria 
 model_6_after_equ_ini <- Csnippet("
                        
                        Ma = 1106;
                        Sn = 2;
                        Sj = 73;
                        Sm = 5;
                        Sf = 5;
                        
                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;
                        
                        In = 9;
                        Ij = 156;
                        Im = 33;
                        If = 33;
                        
                        Rn = 4;
                        Rj = 167;
                        Rm = 1377;
                        Rf = 1377;
                        
                        ")
   