inits<- function(){
  init_states <- Csnippet("
                        
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
  return(init_states)
}

states <- function(){

  state_names <- c("Ma",
                 "Sn","Sj", "Sm", "Sf",
                 "En","Ej", "Em", "Ef",
                 "In","Ij", "Im", "If",
                 "Rn","Rj", "Rm", "Rf")
  
  return(state_names)
}

params <- function(){
  param_names <- c( "beta_val", "gamma_val", "omega_val", "omega_m_val", "kappa_val", "rho_val", 
                "epsilon_val", "mu_val", "mj_val", "m_val",  "delta_t",  "c", "s", "phi",
                "zeta", "disp", "d")
  return(param_names)
}