stochStep <- function(){


  stochStep <- Csnippet("
  /* -------------------------- seasonal forcing ----------------------------------------- */
  
  double b = c * exp(-s * (pow(cos(M_PI * t / 365 - phi), 2)));
  
  double c_v;
  
                       if( t <= 365 ) {
                          c_v = c_v1;
                       } else if( t >= 366 && t <= 730 ) {
                          c_v = c_v1 * k1;
                       } else {
                          c_v = c_v1 * k2;
                       } 
  
  double sdrive = c_v * exp(-s_v * (pow(cos(M_PI * t / 365 - phi_v), 2)));
  
  /*---------------------------- pre-calculations and scaling params:------------------------------- */
  
  double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If +
    Rn + Rj + Rm + Rf + Ma);
  
  double gamma;
  if (model_type == 1) {
    gamma = sdrive * gamma_val / delta_t;
  } else {
    gamma = gamma_val / delta_t;
  }
  
  double omega;
  if (model_type == 2) {
    omega = sdrive * omega_val / delta_t;
  } else {
    omega = omega_val / delta_t;
  }
  
  double omega_m = omega_m_val / delta_t;
  
  double kappa = kappa_val;
  
  double rho = rho_val / delta_t;
  
  double epsilon;
  if (model_type == 3) {
    epsilon = sdrive * epsilon_val / delta_t;
  } else {
    epsilon = epsilon_val / delta_t;
  }
  
  double mu = mu_val / delta_t;
  double mj = mj_val / delta_t;
  double m = m_val / delta_t;
  
  /* ------------------------ general rates of changes  ---------------------------------------------- */
  
  double r_births = b / delta_t;
  double n_Sbirths = rbinom(Sf + If, 1 - exp(-r_births));
  double n_MaBirths = rbinom(Rf + Ef, 1 - exp(-r_births));
  
  double r_juv_death = mj * (N / kappa);
  double r_death = m * (N / kappa);
  
  double r_N_J = omega_m;
  double r_J_A = mu;
  
  double beta_val =  R0 *((epsilon_val + m_val)*(gamma_val + m_val + rho_val)
                  -epsilon_val*rho_val) / (kappa_val*(epsilon_val + m_val));
  double beta = beta_val / delta_t;
  
  double r_R_S = omega;
  double r_I_E = rho;
  double r_S_I = beta * (If + Im + In + Ij);
  double r_E_I = epsilon;
  double r_I_R = gamma;
  
  /* --------------------- newborn transitions: ------------------------------------------------------ */
  
  /* --- rates of changing state --- */
  
  double r_Ma_out = r_juv_death + r_N_J;
  double r_Sn_out = r_juv_death + r_N_J + r_S_I;
  double r_En_out = r_juv_death + r_N_J + r_E_I;
  double r_In_out = r_juv_death + r_N_J + r_I_R + r_I_E;
  double r_Rn_out = r_juv_death + r_N_J + r_R_S;
  
  /* numbers changing state */
  
  double n_Ma_out = rbinom(Ma, 1 - exp(-r_Ma_out));
  double n_Sn_out = rbinom(Sn, 1 - exp(-r_Sn_out));
  double n_En_out = rbinom(En, 1 - exp(-r_En_out));
  double n_In_out = rbinom(In, 1 - exp(-r_In_out));
  double n_Rn_out = rbinom(Rn, 1 - exp(-r_Rn_out));
  
  double probs_Ma_out[2] = {
    r_juv_death / r_Ma_out,
    r_N_J / r_Ma_out
  };
  double probs_Sn_out[3] = {
    r_juv_death / r_Sn_out,
    r_N_J / r_Sn_out,
    r_S_I / r_Sn_out
  };
  double probs_En_out[3] = {
    r_juv_death / r_En_out,
    r_N_J / r_En_out,
    r_E_I / r_En_out
  };
  double probs_In_out[4] = {
    r_juv_death / r_In_out,
    r_N_J / r_In_out,
    r_I_R / r_In_out,
    r_I_E / r_In_out
  };
  double probs_Rn_out[3] = {
    r_juv_death / r_Rn_out,
    r_N_J / r_Rn_out,
    r_R_S / r_Rn_out
  };
  
  int multi_Ma_out[2];
  int multi_Sn_out[3];
  int multi_En_out[3];
  int multi_In_out[4];
  int multi_Rn_out[3];
  
  rmultinom(n_Ma_out, probs_Ma_out, 2, & multi_Ma_out[0]);
  rmultinom(n_Sn_out, probs_Sn_out, 3, & multi_Sn_out[0]);
  rmultinom(n_En_out, probs_En_out, 3, & multi_En_out[0]);
  rmultinom(n_In_out, probs_In_out, 4, & multi_In_out[0]);
  rmultinom(n_Rn_out, probs_Rn_out, 3, & multi_Rn_out[0]);
  
  double n_Rn_Sn = multi_Rn_out[2];
  double n_In_En = multi_In_out[3];
  double n_Sn_In = multi_Sn_out[2];
  double n_En_In = multi_En_out[2];
  double n_In_Rn = multi_In_out[2];
  
  double n_Rn_Rj = multi_Rn_out[1];
  double n_In_Ij = multi_In_out[1];
  double n_Sn_Sj = multi_Sn_out[1];
  double n_En_Ej = multi_En_out[1];
  double n_Ma_Sj = multi_Ma_out[1];
  
  /* equations of changing state */
  
  if (Ma + (n_MaBirths) - n_Ma_out <= 0) {
    Ma = 0;
  } else {
    Ma = Ma + (n_MaBirths) - n_Ma_out;
  }
  
  if (Sn + (n_Sbirths) - n_Sn_out + n_Rn_Sn <= 0) {
    Sn = 0;
  } else {
    Sn = Sn + (n_Sbirths) - n_Sn_out + n_Rn_Sn;
  }
  
  if (En - n_En_out + n_In_En <= 0) {
    En = 0;
  } else {
    En = En - n_En_out + n_In_En;
  }
  
  if (In - n_In_out + n_Sn_In + n_En_In <= 0) {
    In = 0;
  } else {
    In = In - n_In_out + n_Sn_In + n_En_In;
  }
  
  if (Rn - n_Rn_out + n_In_Rn <= 0) {
    Rn = 0;
  } else {
    Rn = Rn - n_Rn_out + n_In_Rn;
  }
  
  /* ------------------------------ juvenile transitions: ------------------------------------------------------ */
  
  /* --- rates of changing state --- */
  
  double r_Sj_out = r_juv_death + r_J_A + r_S_I;
  double r_Ej_out = r_juv_death + r_J_A + r_E_I;
  double r_Ij_out = r_juv_death + r_J_A + r_I_R + r_I_E;
  double r_Rj_out = r_juv_death + r_J_A + r_R_S;
  
  /* numbers changing state */
  
  double n_Sj_out = rbinom(Sj, 1 - exp(-r_Sj_out));
  double n_Ej_out = rbinom(Ej, 1 - exp(-r_Ej_out));
  double n_Ij_out = rbinom(Ij, 1 - exp(-r_Ij_out));
  double n_Rj_out = rbinom(Rj, 1 - exp(-r_Rj_out));
  
  double probs_Sj_out[4] = {
    r_juv_death / r_Sj_out,
    r_J_A / (2 * r_Sj_out),
    r_J_A / (2 * r_Sj_out),
    r_S_I / r_Sj_out
  };
  double probs_Ej_out[4] = {
    r_juv_death / r_Ej_out,
    r_J_A / (2 * r_Ej_out),
    r_J_A / (2 * r_Ej_out),
    r_E_I / r_Ej_out
  };
  double probs_Ij_out[5] = {
    r_juv_death / r_Ij_out,
    r_J_A / (2 * r_Ij_out),
    r_J_A / (2 * r_Ij_out),
    r_I_R / r_Ij_out,
    r_I_E / r_Ij_out
  };
  double probs_Rj_out[4] = {
    r_juv_death / r_Rj_out,
    r_J_A / (2 * r_Rj_out),
    r_J_A / (2 * r_Rj_out),
    r_R_S / r_Rj_out
  };
  
  int multi_Sj_out[4];
  int multi_Ej_out[4];
  int multi_Ij_out[5];
  int multi_Rj_out[4];
  
  rmultinom(n_Sj_out, probs_Sj_out, 4, & multi_Sj_out[0]);
  rmultinom(n_Ej_out, probs_Ej_out, 4, & multi_Ej_out[0]);
  rmultinom(n_Ij_out, probs_Ij_out, 5, & multi_Ij_out[0]);
  rmultinom(n_Rj_out, probs_Rj_out, 4, & multi_Rj_out[0]);
  
  double n_Rj_Sj = multi_Rj_out[3];
  double n_Ij_Ej = multi_Ij_out[4];
  double n_Sj_Ij = multi_Sj_out[3];
  double n_Ej_Ij = multi_Ej_out[3];
  double n_Ij_Rj = multi_Ij_out[3];
  
  double n_Rj_Rf = multi_Rj_out[1];
  double n_Ij_If = multi_Ij_out[1];
  double n_Sj_Sf = multi_Sj_out[1];
  double n_Ej_Ef = multi_Ej_out[1];
  
  double n_Rj_Rm = multi_Rj_out[2];
  double n_Ij_Im = multi_Ij_out[2];
  double n_Sj_Sm = multi_Sj_out[2];
  double n_Ej_Em = multi_Ej_out[2];
  
  /* equations of changing state */
  
  if (Sj - n_Sj_out + n_Ma_Sj + n_Sn_Sj + n_Rj_Sj <= 0) {
    Sj = 0;
  } else {
    Sj = Sj - n_Sj_out + n_Ma_Sj + n_Sn_Sj + n_Rj_Sj;
  }
  
  if (Ej - n_Ej_out + n_En_Ej + n_Ij_Ej <= 0) {
    Ej = 0;
  } else {
    Ej = Ej - n_Ej_out + n_En_Ej + n_Ij_Ej;
  }
  
  if (Ij - n_Ij_out + n_In_Ij + n_Sj_Ij + n_Ej_Ij <= 0) {
    Ij = 0;
  } else {
    Ij = Ij - n_Ij_out + n_In_Ij + n_Sj_Ij + n_Ej_Ij;
  }
  
  if (Rj - n_Rj_out + n_Rn_Rj + n_Ij_Rj <= 0) {
    Rj = 0;
  } else {
    Rj = Rj - n_Rj_out + n_Rn_Rj + n_Ij_Rj;
  }
  
  /* ------------------------------ adult male transitions: ------------------------------------------------------ */
  
  /* --- rates of changing state --- */
  
  double r_Sm_out = r_death + r_S_I;
  double r_Em_out = r_death + r_E_I;
  double r_Im_out = r_death + r_I_R + r_I_E;
  double r_Rm_out = r_death + r_R_S;
  
  /* numbers changing state */
  
  double n_Sm_out = rbinom(Sm, 1 - exp(-r_Sm_out));
  double n_Em_out = rbinom(Em, 1 - exp(-r_Em_out));
  double n_Im_out = rbinom(Im, 1 - exp(-r_Im_out));
  double n_Rm_out = rbinom(Rm, 1 - exp(-r_Rm_out));
  
  double probs_Sm_out[2] = {
    r_death / r_Sm_out,
    r_S_I / r_Sm_out
  };
  double probs_Em_out[2] = {
    r_death / r_Em_out,
    r_E_I / r_Em_out
  };
  double probs_Im_out[3] = {
    r_death / r_Im_out,
    r_I_R / r_Im_out,
    r_I_E / r_Im_out
  };
  double probs_Rm_out[2] = {
    r_death / r_Rm_out,
    r_R_S / r_Rm_out
  };
  
  int multi_Sm_out[2];
  int multi_Em_out[2];
  int multi_Im_out[3];
  int multi_Rm_out[2];
  
  rmultinom(n_Sm_out, probs_Sm_out, 2, & multi_Sm_out[0]);
  rmultinom(n_Em_out, probs_Em_out, 2, & multi_Em_out[0]);
  rmultinom(n_Im_out, probs_Im_out, 3, & multi_Im_out[0]);
  rmultinom(n_Rm_out, probs_Rm_out, 2, & multi_Rm_out[0]);
  
  double n_Rm_Sm = multi_Rm_out[1];
  double n_Im_Em = multi_Im_out[2];
  double n_Sm_Im = multi_Sm_out[1];
  double n_Em_Im = multi_Em_out[1];
  double n_Im_Rm = multi_Im_out[1];
  
  /* equations of changing state */
  
  if (Sm - n_Sm_out + n_Sj_Sm + n_Rm_Sm <= 0) {
    Sm = 0;
  } else {
    Sm = Sm - n_Sm_out + n_Sj_Sm + n_Rm_Sm;
  }
  
  if (Em - n_Em_out + n_Ej_Em + n_Im_Em <= 0) {
    Em = 0;
  } else {
    Em = Em - n_Em_out + n_Ej_Em + n_Im_Em;
  }
  
  if (Im - n_Im_out + n_Ij_Im + n_Sm_Im + n_Em_Im <= 0) {
    Im = 0;
  } else {
    Im = Im - n_Im_out + n_Ij_Im + n_Sm_Im + n_Em_Im;
  }
  
  if (Rm - n_Rm_out + n_Rj_Rm + n_Im_Rm <= 0) {
    Rm = 0;
  } else {
    Rm = Rm - n_Rm_out + n_Rj_Rm + n_Im_Rm;
  }
  
  /* ------------------------------ adult female transitions: ------------------------------------------------------ */
  
  /* --- rates of changing state --- */
  
  double r_Sf_out = r_death + r_S_I;
  double r_Ef_out = r_death + r_E_I;
  double r_If_out = r_death + r_I_R + r_I_E;
  double r_Rf_out = r_death + r_R_S;
  
  /* numbers changing state */
  
  double n_Sf_out = rbinom(Sf, 1 - exp(-r_Sf_out));
  double n_Ef_out = rbinom(Ef, 1 - exp(-r_Ef_out));
  double n_If_out = rbinom(If, 1 - exp(-r_If_out));
  double n_Rf_out = rbinom(Rf, 1 - exp(-r_Rf_out));
  
  double probs_Sf_out[2] = {
    r_death / r_Sf_out,
    r_S_I / r_Sf_out
  };
  double probs_Ef_out[2] = {
    r_death / r_Ef_out,
    r_E_I / r_Ef_out
  };
  double probs_If_out[3] = {
    r_death / r_If_out,
    r_I_R / r_If_out,
    r_I_E / r_If_out
  };
  double probs_Rf_out[2] = {
    r_death / r_Rf_out,
    r_R_S / r_Rf_out
  };
  
  int multi_Sf_out[2];
  int multi_Ef_out[2];
  int multi_If_out[3];
  int multi_Rf_out[2];
  
  rmultinom(n_Sf_out, probs_Sf_out, 2, & multi_Sf_out[0]);
  rmultinom(n_Ef_out, probs_Ef_out, 2, & multi_Ef_out[0]);
  rmultinom(n_If_out, probs_If_out, 3, & multi_If_out[0]);
  rmultinom(n_Rf_out, probs_Rf_out, 2, & multi_Rf_out[0]);
  
  double n_Rf_Sf = multi_Rf_out[1];
  double n_If_Ef = multi_If_out[2];
  double n_Sf_If = multi_Sf_out[1];
  double n_Ef_If = multi_Ef_out[1];
  double n_If_Rf = multi_If_out[1];
  
  /* equations of changing state */
  
  if (Sf - n_Sf_out + n_Sj_Sf + n_Rf_Sf <= 0) {
    Sf = 0;
  } else {
    Sf = Sf - n_Sf_out + n_Sj_Sf + n_Rf_Sf;
  }
  
  if (Ef - n_Ef_out + n_Ej_Ef + n_If_Ef <= 0) {
    Ef = 0;
  } else {
    Ef = Ef - n_Ef_out + n_Ej_Ef + n_If_Ef;
  }
  
  if (If - n_If_out + n_Ij_If + n_Sf_If + n_Ef_If <= 0) {
    If = 0;
  } else {
    If = If - n_If_out + n_Ij_If + n_Sf_If + n_Ef_If;
  }
  
  if (Rf - n_Rf_out + n_Rj_Rf + n_If_Rf <= 0) {
    Rf = 0;
  } else {
    Rf = Rf - n_Rf_out + n_Rj_Rf + n_If_Rf;
  }
  
  /* ------------------------------ accumvars new infected: ------------------------------------------------------ */
  if (t == 0) {
    H = 0;
  } else {
  H =  n_Sn_In + n_En_In 
        + n_In_Ij + n_Sj_Ij + n_Ej_Ij 
        + n_Ij_Im + n_Sm_Im + n_Em_Im 
        + n_Ij_If + n_Sf_If + n_Ef_If;
  }                      
  ")
  return(stochStep)
    
}
