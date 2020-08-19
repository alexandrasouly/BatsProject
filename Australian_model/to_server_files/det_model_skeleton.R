det_model_skeleton <- function(){
  
  det_model_skeleton <- Csnippet("
                      
                      /* seasonal forcing of birth */

                      double b = c * exp(-s * (pow(cos(M_PI * t / 365 - phi), 2)));
                      
                      double sdrive = c_v * exp(-s_v * (pow(cos(M_PI * t / 365 - phi_v), 2)));
                      
                      /* pre-calculations and scaling params: */
                      
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
                      
                      double Sbirths = b * (Sf + If) / delta_t;
                      double MaBirths = b * (Rf + Ef) / delta_t;
                      
                      double beta_val =  R0 *((epsilon_val + m_val)*(gamma_val + m_val + rho_val)
                                          -epsilon_val*rho_val) / (kappa_val*(epsilon_val + m_val));
                      double beta_v = beta_val / delta_t;
                      
                      double betSN = omega_m + mj * (N / kappa) + beta_v * (If + Im + Ij + In);
                      double betSj = mu + mj * (N / kappa) + (beta_v) * (If + Im + In + Ij);
                      double betSmf = m * (N / kappa) + (beta_v) * (If + Im + In + Ij);
                      double betI = beta_v * (If + Im + In + Ij);
                      
                      /* newborn transitions: */
                      
                      if (Ma + (MaBirths) - (omega_m + mj * (N / kappa)) * Ma <= 0) {
                        DMa = 0;
                      } else {
                        DMa = Ma + (MaBirths) - (omega_m + mj * (N / kappa)) * Ma;
                      }
                      
                      if (Sn + (Sbirths) - (betSN * Sn) + omega * Rn <= 0) {
                        DSn = 0;
                      } else {
                        DSn = Sn + (Sbirths) - (betSN * Sn) + omega * Rn;
                      }
                      
                      if (En - ((omega_m + mj * (N / kappa) + epsilon) * En) + rho * In <= 0) {
                        DEn = 0;
                      } else {
                        DEn = En - ((omega_m + mj * (N / kappa) + epsilon) * En) + rho * In;
                      }
                      
                      if (In - ((omega_m + mj * (N / kappa) + gamma + rho) * In) + (betI * Sn) + epsilon * En <= 0) {
                        DIn = 0;
                      } else {
                        DIn = In - ((omega_m + mj * (N / kappa) + gamma + rho) * In) + (betI * Sn) + epsilon * En;
                      }
                      
                      if (Rn - (omega_m + mj * (N / kappa) + omega) * Rn + gamma * In <= 0) {
                        DRn = 0;
                      } else {
                        DRn = Rn - (omega_m + mj * (N / kappa) + omega) * Rn + gamma * In;
                      }
                      
                      /* juvenile transitions: */
                      
                      if (Sj + omega_m * (Sn + Ma) - betSj * Sj + omega * Rj <= 0) {
                        DSj = 0;
                      } else {
                        DSj = Sj + omega_m * (Sn + Ma) - betSj * Sj + omega * Rj;
                      }
                      
                      if (Ej + omega_m * En - (mu + mj * (N / kappa) + epsilon) * Ej + rho * Ij <= 0) {
                        DEj = 0;
                      } else {
                        DEj = Ej + omega_m * En - (mu + mj * (N / kappa) + epsilon) * Ej + rho * Ij;
                      }
                      
                      if (Ij + omega_m * In - (mu + mj * (N / kappa) + gamma + rho) * Ij + betI * Sj + epsilon * Ej <= 0) {
                        DIj = 0;
                      } else {
                        DIj = Ij + omega_m * In - (mu + mj * (N / kappa) + gamma + rho) * Ij + betI * Sj + epsilon * Ej;
                      }
                      
                      if (Rj + omega_m * Rn - (mu + mj * (N / kappa) + omega) * Rj + gamma * Ij <= 0) {
                        DRj = 0;
                      } else {
                        DRj = Rj + omega_m * Rn - (mu + mj * (N / kappa) + omega) * Rj + gamma * Ij;
                      }
                      
                      /* adult male transitions: */
                      
                      if (Sm + mu * (Sj / 2) - betSmf * Sm + omega * Rm <= 0) {
                        DSm = 0;
                      } else {
                        DSm = Sm + mu * (Sj / 2) - betSmf * Sm + omega * Rm;
                      }
                      
                      if (Em + mu * (Ej / 2) - (m * (N / kappa) + epsilon) * Em + rho * Im <= 0) {
                        DEm = 0;
                      } else {
                        DEm = Em + mu * (Ej / 2) - (m * (N / kappa) + epsilon) * Em + rho * Im;
                      }
                      
                      if (Im + mu * (Ij / 2) - (m * (N / kappa) + gamma + rho) * Im + (betI * Sm) + epsilon * Em <= 0) {
                        DIm = 0;
                      } else {
                        DIm = Im + mu * (Ij / 2) - (m * (N / kappa) + gamma + rho) * Im + (betI * Sm) + epsilon * Em;
                      }
                      
                      if (Rm + mu * (Rj / 2) - (m * (N / kappa) + omega) * Rm + gamma * Im <= 0) {
                        DRm = 0;
                      } else {
                        DRm = Rm + mu * (Rj / 2) - (m * (N / kappa) + omega) * Rm + gamma * Im;
                      }
                      
                      /* adult female transitions: */
                      
                      if (Sf + mu * (Sj / 2) - betSmf * Sf + omega * Rf <= 0) {
                        DSf = 0;
                      } else {
                        DSf = Sf + mu * (Sj / 2) - betSmf * Sf + omega * Rf;
                      }
                      
                      if (Ef + mu * (Ej / 2) - (m * (N / kappa) + epsilon) * Ef + rho * If <= 0) {
                        DEf = 0;
                      } else {
                        DEf = Ef + mu * (Ej / 2) - (m * (N / kappa) + epsilon) * Ef + rho * If;
                      }
                      
                      if (If + mu * (Ij / 2) - (m * (N / kappa) + gamma + rho) * If + (betI * Sf) + epsilon * Ef <= 0) {
                        DIf = 0;
                      } else {
                        DIf = If + mu * (Ij / 2) - (m * (N / kappa) + gamma + rho) * If + (betI * Sf) + epsilon * Ef;
                      }
                      
                      if (Rf + mu * (Rj / 2) - (m * (N / kappa) + omega) * Rf + gamma * If <= 0) {
                        DRf = 0;
                      } else {
                        DRf = Rf + mu * (Rj / 2) - (m * (N / kappa) + omega) * Rf + gamma * If;
                      }
                      
                      if (t == 0){
                      DH = 0;
                      } else {
                      DH = 
                                          betI * Sn + epsilon * En 
                         + omega_m * In  + betI * Sj + epsilon * Ej
                         + mu * (Ij / 2) + betI * Sm + epsilon * Em
                         + mu * (Ij / 2) + betI * Sf + epsilon * Ef;
                      }   
                          
                          ")
  return(det_model_skeleton)
  
}