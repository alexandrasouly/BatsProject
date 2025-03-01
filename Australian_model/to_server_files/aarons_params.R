



#---------------------------------- MODEL 5 -------------------------------------------
#                      model 5 is SIR with no seasonal forcing
# posterior parameters from Aaron's draft
# model_5_params <- c(
# 
#   # base dynamics params
#   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
# 
#   R0 =  42.152,
#   gamma_val = 3.135, # recovery rate I -> R
#   omega_val = 0, # immune waning rate R -> S
#   omega_m_val = 0.8, # maternal antibody waning rate
#   kappa_val = 5000, # carrying capacity for Clunes
#   rho_val = 0, # I -> L
#   epsilon_val = 0, # L -> I
# 
# 
#   mu_val = 1.37, # juvenile maturation rate
#   mj_val = 0.50, # juvenile death rate
#   m_val = 0.23, # adult death rate
#   delta_t = 365, # scaling time as days instead of years
# 
#   # seasonality params
#   c = 18.65, # birth pulse scaling factor
#   s = 130, # birth pulse synchronicity
#   phi = 7.19, # birth pulse time shift
# 
#   c_v = 1, # seasonal drive scaling factor
#   s_v = 0, # seasonal drive synchronicity
#   phi_v = 0, # sesonal drive time shift
# 
#   # measuring process params
#   zeta = 0.3, # test accuracy
#   disp = 1000, # dispersion parameter
#   d = 10 # number of bats contributing to the same pool
# 
# )

# #Aaron's 95% confidence interval lower bounds
# model_5_lower <- c(
#   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
#
#   R0 = , # infection rate S -> R
#   gamma_val = 2.094, # recovery rate I -> R
#   omega_val = 0, # immune waning rate R -> S
#   omega_m_val = 0.749, # maternal antibody waning rate
#   kappa_val = 5000, # carrying capacity
#   rho_val = 0, # I -> L
#   epsilon_val = 0, # L -> I
#
#
#   mu_val = 1.37, # juvenile maturation rate
#   mj_val = 0.50, # juvenile death rate
#   m_val = 0.23, # adult death rate
#   delta_t = 365, # scaling time as days instead of years
#
#   # seasonality params
#   c = 13.68, # birth pulse scaling factor
#   s = 112.701, # birth pulse synchronicity
#   phi = 7.17, # birth pulse time shift
#
#   c_v = 1, # seasonal drive scaling factor
#   s_v = 0, # seasonal drive synchronicity
#   phi_v = 0, # sesonal drive time shift
#
#   # measuring process params
#   zeta = 0.261, # test accuracy
#   disp = 1000, # dispersion parameter
#   d = 4 # number of bats contributing to the same pool
#
# )
#
#
# model_5_upper <- c(
#   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
#
#   R0 = , # infection rate S -> R
#   gamma_val = 4.036, # recovery rate I -> R
#   omega_val = 0, # immune waning rate R -> S
#   omega_m_val = 0.859, # maternal antibody waning rate
#   kappa_val = 5000, # carrying capacity
#   rho_val = 0, # I -> L
#   epsilon_val = 0, # L -> I
#
#
#   mu_val = 1.37, # juvenile maturation rate
#   mj_val = 0.50, # juvenile death rate
#   m_val = 0.23, # adult death rate
#   delta_t = 365, # scaling time as days instead of years
#
#   # seasonality params
#   c = 21.359, # birth pulse scaling factor
#   s = 151.075, # birth pulse synchronicity
#   phi = 7.2, # birth pulse time shift
#
#   c_v = 1, # seasonal drive scaling factor
#   s_v = 0, # seasonal drive synchronicity
#   phi_v = 0, # sesonal drive time shift
#
#   # measuring process params
#   zeta = 0.819, # test accuracy
#   disp = 1000, # dispersion parameter
#   d = 19 # number of bats contributing to the same pool
#
#
# )



# initial states from Aaron's draft
 # model_5_before_equ_ini <- Csnippet("
 #                        
 #                        Ma = 0;
 #                        Sn = 1264;
 #                        Sj = 541;
 #                        Sm = 1597;
 #                        Sf = 1597;
 #                        
 #                        En = 0;
 #                        Ej = 0;
 #                        Em = 0;
 #                        Ef = 0;
 #                        
 #                        In = 0;
 #                        Ij = 0;
 #                        Im = 125;
 #                        If = 125;
 #                        
 #                        Rn = 0;
 #                        Rj = 0;
 #                        Rm = 0;
 #                        Rf = 0;

 #                        H = 0;
 #                        
 #                        ")
 # 
 # 50 years waiting for equilibria 
 # model_5_after_equ_ini <- Csnippet("
 #                        
 #                        Ma = 1459;
 #                        Sn = 10;
 #                        Sj = 145;
 #                        Sm = 13;
 #                        Sf = 13;
 #                        
 #                        En = 0;
 #                        Ej = 0;
 #                        Em = 0;
 #                        Ef = 0;
 #                        
 #                        In = 18;
 #                        Ij = 138;
 #                        Im = 42;
 #                        If = 42;
 #                        
 #                        Rn = 21;
 #                        Rj = 227;
 #                        Rm = 1481;
 #                        Rf = 1481;
 #                        
 #                        ")
 # # 
 # 
 # #---------------------------------- MODEL 6 -------------------------------------------
 # #                      model 6 is SIR with seasonal forcing 
 # 
 # # posterior parameters from Aaron's draft
 # model_6_params <- c(
 #   
 #   # base dynamics params
 #   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only 
 #   
 #   R0 = NA, 
 #   gamma_val = 22.449, # recovery rate I -> R
 #   omega_val = 0, # immune waning rate R -> S
 #   omega_m_val = 0.801, # maternal antibody waning rate
 #   kappa_val = 4216, # carrying capacity
 #   rho_val = 0, # I -> L 
 #   epsilon_val = 0, # L -> I 
 #   mu_val = 1.37, # juvenile maturation rate
 #   mj_val = 0.499, # juvenile death rate
 #   m_val = 0.187, # adult death rate
 #   delta_t = 365, # scaling time as days instead of years
 #   
 #   # seasonality params
 #   c = 15.973, # birth pulse scaling factor
 #   s = 129.841, # birth pulse synchronicity
 #   phi = 7.181, # birth pulse time shift
 #   
 #   c_v = 1, # seasonal drive scaling factor
 #   s_v = 26.581, # seasonal drive synchronicity
 #   phi_v = 0.573, # sesonal drive time shift
 #   
 #   # measuring process params
 #   zeta = 0.433, # test accuracy
 #   disp = 1000, # dispersion parameter
 #   d = 8.2 # number of bats contributing to the same pool
 #   
 # )
 # 
 # # initial states from Aaron's draft
 # model_6_before_equ_ini <- Csnippet("
 #                        
 #                        Ma = 0;
 #                        Sn = 922;
 #                        Sj = 396;
 #                        Sm = 1449;
 #                        Sf = 1449;
 #                        
 #                        En = 0;
 #                        Ej = 0;
 #                        Em = 0;
 #                        Ef = 0;
 #                        
 #                        In = 0;
 #                        Ij = 0;
 #                        Im = 105;
 #                        If = 105;
 #                        
 #                        Rn = 0;
 #                        Rj = 0;
 #                        Rm = 0;
 #                        Rf = 0;
 #                        H = 0;
 #                        
 #                        ")
 # 
 # # 50 years waiting for equilibria 
 # model_6_after_equ_ini <- Csnippet("
 #                        
 #                        Ma = 1106;
 #                        Sn = 2;
 #                        Sj = 73;
 #                        Sm = 5;
 #                        Sf = 5;
 #                        
 #                        En = 0;
 #                        Ej = 0;
 #                        Em = 0;
 #                        Ef = 0;
 #                        
 #                        In = 9;
 #                        Ij = 156;
 #                        Im = 33;
 #                        If = 33;
 #                        
 #                        Rn = 4;
 #                        Rj = 167;
 #                        Rm = 1377;
 #                        Rf = 1377;
 #                        
 #                        ")
 # 
 # #---------------------------------- MODEL 4 -------------------------------------------
 #                model 4 is SILI with seasonal forcing and mat immun for L
 # 
 # # posterior parameters from Aaron's draft
 # model_4_params <- c(
 # 
 #   # base dynamics params
 #   model_type = 3, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
 # 
 #   R0 = 7.659, # infection rate S -> R
 #   gamma_val = 0, # recovery rate I -> R
 #   omega_val = 0, # immune waning rate R -> S
 #   omega_m_val = 0.799, # maternal antibody waning rate
 #   kappa_val = 5000, # carrying capacity
 #   rho_val = 59.22, # I -> L
 #   epsilon_val = 3.799, # L -> I
 #   mu_val = 1.37, # juvenile maturation rate
 #   mj_val = 0.5, # juvenile death rate
 #   m_val = 0.187, # adult death rate
 #   delta_t = 365, # scaling time as days instead of years
 # 
 #   # seasonality params
 #   c = 16.068, # birth pulse scaling factor
 #   s = 129.744, # birth pulse synchronicity
 #   phi = 7.18, # birth pulse time shift
 # 
 #   c_v = 1, # seasonal drive scaling factor
 #   s_v = 1.861, # seasonal drive synchronicity
 #   phi_v = 0.034, # sesonal drive time shift
 # 
 #   # measuring process params
 #   zeta = 0.661, # test accuracy
 #   disp = 1000, # dispersion parameter
 #   d = 12 # number of bats contributing to the same pool
 # 
 # )
 # 
 # # initial states from Aaron's draft
 # model_4_before_equ_ini <- Csnippet("
 # 
 #                        Ma = 0;
 #                        Sn = 1097;
 #                        Sj = 469;
 #                        Sm = 1717;
 #                        Sf = 1717;
 # 
 #                        En = 0;
 #                        Ej = 0;
 #                        Em = 0;
 #                        Ef = 0;
 # 
 #                        In = 0;
 #                        Ij = 0;
 #                        Im = 125;
 #                        If = 125;
 # 
 #                        Rn = 0;
 #                        Rj = 0;
 #                        Rm = 0;
 #                        Rf = 0;
 #                        
 #                        H = 0;
 # 
 #                        ")
 # model_4_params_lower <- c(
 #    
 #    # base dynamics params
 #    model_type = 3, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
 #    
 #    R0 = 3.393, # infection rate S -> R
 #    gamma_val = 0, # recovery rate I -> R
 #    omega_val = 0, # immune waning rate R -> S
 #    omega_m_val = 0.799, # maternal antibody waning rate
 #    kappa_val = 5000, # carrying capacity
 #    rho_val = 19.409, # I -> L
 #    epsilon_val = 0.2, # L -> I
 #    mu_val = 1.37, # juvenile maturation rate
 #    mj_val = 0.5, # juvenile death rate
 #    m_val = 0.187, # adult death rate
 #    delta_t = 365, # scaling time as days instead of years
 #    
 #    # seasonality params
 #    c = 12.405, # birth pulse scaling factor
 #    s = 110.155, # birth pulse synchronicity
 #    phi = 7.18, # birth pulse time shift
 #    
 #    c_v = 1, # seasonal drive scaling factor
 #    s_v = 0.906, # seasonal drive synchronicity
 #    phi_v = 0, # sesonal drive time shift
 #    
 #    # measuring process params
 #    zeta = 0.278, # test accuracy
 #    disp = 10, # dispersion parameter
 #    d = 2 # number of bats contributing to the same pool
 #    
 # )
 # 
 # model_4_params_upper <- c(
 #    
 #    # base dynamics params
 #    model_type = 3, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
 #    
 #    R0 = 18.199, # infection rate S -> R
 #    gamma_val = 0, # recovery rate I -> R
 #    omega_val = 0, # immune waning rate R -> S
 #    omega_m_val = 0.799, # maternal antibody waning rate
 #    kappa_val = 5000, # carrying capacity
 #    rho_val = 306.196, # I -> L
 #    epsilon_val = 25.351, # L -> I
 #    mu_val = 1.37, # juvenile maturation rate
 #    mj_val = 0.5, # juvenile death rate
 #    m_val = 0.187, # adult death rate
 #    delta_t = 365, # scaling time as days instead of years
 #    
 #    # seasonality params
 #    c = 19.664, # birth pulse scaling factor
 #    s = 148.909, # birth pulse synchronicity
 #    phi = 7.18, # birth pulse time shift
 #    
 #    c_v = 1, # seasonal drive scaling factor
 #    s_v = 2.655, # seasonal drive synchronicity
 #    phi_v = 0.136, # sesonal drive time shift
 #    
 #    # measuring process params
 #    zeta = 0.999, # test accuracy
 #    disp = 1000, # dispersion parameter
 #    d = 20 # number of bats contributing to the same pool
 #    )
 #    
 #    # #---------------------------------- MODEL 8 -------------------------------------------
 #    # #                      model 8 is SIRS with seasonal forcing 
 #    #
 #    # posterior parameters from Aaron's draft
 #    model_8_params <- c(
 # 
 #      # base dynamics params
 #      model_type = 2, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
 # 
 #      R0 = 35.426,
 #      gamma_val = 4.653, # recovery rate I -> R
 #      omega_val = 5.701, # immune waning rate R -> S
 #      omega_m_val = 0.801, # maternal antibody waning rate
 #      kappa_val = 5000, # carrying capacity
 #      rho_val = 0, # I -> L
 #      epsilon_val = 0, # L -> I
 #      mu_val = 1.37, # juvenile maturation rate
 #      mj_val = 0.501, # juvenile death rate
 #      m_val = 0.191, # adult death rate
 #      delta_t = 365, # scaling time as days instead of years
 # 
 #      # seasonality params
 #      c = 16.301, # birth pulse scaling factor
 #      s = 129.979, # birth pulse synchronicity
 #      phi = 7.179, # birth pulse time shift
 # 
 #      c_v = 1, # seasonal drive scaling factor
 #      s_v = 125.843, # seasonal drive synchronicity
 #      phi_v = 1.75, # sesonal drive time shift
 # 
 #      # measuring process params
 #      zeta = 0.356, # test accuracy
 #      disp = 100, # dispersion parameter
 #      d = 5.221 # number of bats contributing to the same pool
 # 
 #    )
    #
    # initial states from Aaron's draft
    model_8_before_equ_ini <- Csnippet("

                        Ma = 0;
                        Sn = 1097;
                        Sj = 469;
                        Sm = 1717;
                        Sf = 1717;

                        En = 0;
                        Ej = 0;
                        Em = 0;
                        Ef = 0;

                        In = 0;
                        Ij = 0;
                        Im = 125;
                        If = 125;

                        Rn = 0;
                        Rj = 0;
                        Rm = 0;
                        Rf = 0;
                        
                        H = 0;

                        ")
    

    model_8_params_lower <- c(

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
      s = 109.237, # birth pulse synchronicity
      phi = 7.179, # birth pulse time shift

      c_v = 1, # seasonal drive scaling factor
      s_v = 37.69, # seasonal drive synchronicity
      phi_v = 1.473 , # sesonal drive time shift

      # measuring process params
      zeta = 0.046, # test accuracy
      disp = 10, # dispersion parameter
      d = 0.01 # number of bats contributing to the same pool

      )

    model_8_params_upper <- c(

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
      s = 148.619, # birth pulse synchronicity
      phi = 7.179, # birth pulse time shift

      c_v = 1, # seasonal drive scaling factor
      s_v = 199.776, # seasonal drive synchronicity
      phi_v = 2.111 , # sesonal drive time shift

      # measuring process params
      zeta = 0.903, # test accuracy
      disp = 1000, # dispersion parameter
      d = 17.558 # number of bats contributing to the same pool

    )
    
    # # #---------------------------------- MODEL 6 -------------------------------------------
    # #                model 6 is SIR with seasonal forcing
    # # 
    # # # posterior parameters from Aaron's draft
    # model_6_params <- c(
    # 
    #   # base dynamics params
    #   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
    # 
    #   R0 = 13.395, # infection rate S -> R
    #   gamma_val = 22.449, # recovery rate I -> R
    #   omega_val = 0, # immune waning rate R -> S
    #   omega_m_val = 0.801, # maternal antibody waning rate
    #   kappa_val = 5000, # carrying capacity
    #   rho_val = 0, # I -> L
    #   epsilon_val = 0, # L -> I
    #   mu_val = 1.37, # juvenile maturation rate
    #   mj_val = 0.5, # juvenile death rate
    #   m_val = 0.187, # adult death rate
    #   delta_t = 365, # scaling time as days instead of years
    # 
    #   # seasonality params
    #   c = 15.973, # birth pulse scaling factor
    #   s = 129.841, # birth pulse synchronicity
    #   phi = 7.181, # birth pulse time shift
    # 
    #   c_v = 1, # seasonal drive scaling factor
    #   s_v = 26.581, # seasonal drive synchronicity
    #   phi_v = 0.573, # sesonal drive time shift
    # 
    #   # measuring process params
    #   zeta = 0.443, # test accuracy
    #   disp = 100, # dispersion parameter
    #   d = 8.246 # number of bats contributing to the same pool
    # 
    # )
    # 
    # # initial states from Aaron's draft
    # model_6_before_equ_ini <- Csnippet("
    # 
    #                     Ma = 0;
    #                     Sn = 1097;
    #                     Sj = 469;
    #                     Sm = 1717;
    #                     Sf = 1717;
    # 
    #                     En = 0;
    #                     Ej = 0;
    #                     Em = 0;
    #                     Ef = 0;
    # 
    #                     In = 0;
    #                     Ij = 0;
    #                     Im = 125;
    #                     If = 125;
    # 
    #                     Rn = 0;
    #                     Rj = 0;
    #                     Rm = 0;
    #                     Rf = 0;
    #                     
    #                     H = 0;
    # 
    #                     ")
    # model_6_params_lower <- c(
    # 
    #   # base dynamics params
    #   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
    # 
    #   R0 = 7.087, # infection rate S -> R
    #   gamma_val = 11.455, # recovery rate I -> R
    #   omega_val = 0, # immune waning rate R -> S
    #   omega_m_val = 0.801, # maternal antibody waning rate
    #   kappa_val = 5000, # carrying capacity
    #   rho_val = 0, # I -> L
    #   epsilon_val = 0, # L -> I
    #   mu_val = 1.37, # juvenile maturation rate
    #   mj_val = 0.5, # juvenile death rate
    #   m_val = 0.187, # adult death rate
    #   delta_t = 365, # scaling time as days instead of years
    # 
    #   # seasonality params
    #   c = 12.262, # birth pulse scaling factor
    #   s = 110.473, # birth pulse synchronicity
    #   phi = 7.181, # birth pulse time shift
    # 
    #   c_v = 1, # seasonal drive scaling factor
    #   s_v = 5.96, # seasonal drive synchronicity
    #   phi_v = 0.389, # sesonal drive time shift
    # 
    #   # measuring process params
    #   zeta = 0.133, # test accuracy
    #   disp = 10, # dispersion parameter
    #   d = 1.751# number of bats contributing to the same pool
    # 
    # )
    # 
    # model_6_params_upper <- c(
    # 
    #   # base dynamics params
    #   model_type = 1, #  1: SIR, 2:SIRS or 3:SILI, matters for seasonal forcing only
    # 
    #   R0 = 22.842, # infection rate S -> R
    #   gamma_val = 48.417, # recovery rate I -> R
    #   omega_val = 0, # immune waning rate R -> S
    #   omega_m_val = 0.801, # maternal antibody waning rate
    #   kappa_val = 5000, # carrying capacity
    #   rho_val = 0, # I -> L
    #   epsilon_val = 0, # L -> I
    #   mu_val = 1.37, # juvenile maturation rate
    #   mj_val = 0.5, # juvenile death rate
    #   m_val = 0.187, # adult death rate
    #   delta_t = 365, # scaling time as days instead of years
    # 
    #   # seasonality params
    #   c = 19.695, # birth pulse scaling factor
    #   s = 149.831, # birth pulse synchronicity
    #   phi = 7.181, # birth pulse time shift
    # 
    #   c_v = 1, # seasonal drive scaling factor
    #   s_v = 149.089, # seasonal drive synchronicity
    #   phi_v = 0.904, # sesonal drive time shift
    # 
    #   # measuring process params
    #   zeta = 0.925, # test accuracy
    #   disp = 1000, # dispersion parameter
    #   d = 18.817 # number of bats contributing to the same pool
    # )
    # 
