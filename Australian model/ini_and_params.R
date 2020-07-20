inits<- function(){
  init_states <- Csnippet("
                        
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