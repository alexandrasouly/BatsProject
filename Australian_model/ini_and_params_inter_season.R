

states <- function(){

  state_names <- c("Ma",
                 "Sn","Sj", "Sm", "Sf",
                 "En","Ej", "Em", "Ef",
                 "In","Ij", "Im", "If",
                 "Rn","Rj", "Rm", "Rf",
                 "H")
  
  return(state_names)
}

params <- function(){
  param_names <- c("model_type", "R0", "gamma_val", "omega_val", "omega_m_val", "kappa_val", "rho_val", 
                "epsilon_val", "mu_val", "mj_val", "m_val",  "delta_t",  "c", "s", "phi", "c_v1", "c_v2", "c_v3", "s_v", "phi_v",
                "zeta", "disp", "d")
  return(param_names)
}