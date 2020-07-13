# full dynamics SIR/SIRS/SILI
# deterministic equations
# same equs as in Aaron's supp material
# maternal immunity included
# here for now E and R gives maternal immunity as in supp eqs, not github equs

library(pomp)

det_model_skeleton <- Csnippet("

  /* pre-calculations and scaling params: */
  
  double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If 
    + Rn + Rj + Rm + Rf + Ma);
    
  double b = b_val / delta_t;
  double beta_v = beta_val / delta_t;
  double gamma = gamma_val / delta_t;
  double omega = omega_val / delta_t;
  double omega_m = omega_m_val / delta_t;
  double kappa = kappa_val / delta_t;
  double rho = rho_val / delta_t;
  double epsilon = epsilon_val / delta_t;
  double mu = mu_val / delta_t;
  double mj = mj_val / delta_t;
  double m = m_val /delta_t;
    

  double Sbirths = b * (Sf + If);
  double MaBirths = b * (Rf + Ef);
  
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
  

")


init_states <- Csnippet("
  
  Ma = 0;

  Sn = 20;
  Sj = 20;
  Sm = 40;
  Sf = 40;
  
  En = 0;
  Ej = 0;
  Em = 0;
  Ef = 0;
  
  In = 5;
  Ij = 5;
  Im = 10;
  If = 10;
  
  Rn = 0;
  Rj = 0;
  Rm = 0;
  Rf = 0;
  
  
")

# try SIR here: rho, epsilon = 0
params <- c(
  
  b_val = 0.6,
  beta_val = 1,
  gamma_val = 0.5,
  omega_val = 0,
  omega_m_val = 0.4,
  kappa_val = 1000,
  rho_val = 0,
  epsilon_val = 0,
  mu_val = 0.44,
  mj_val = 0.6,
  m_val = 0.2,
  delta_t = 365)

state_names <- c("Ma",
               "Sn","Sj", "Sm", "Sf",
               "En","Ej", "Em", "Ef",
               "In","Ij", "Im", "If",
               "Rn","Rj", "Rm", "Rf"
               )


param_names <- c( "b_val", "beta_val", "gamma_val", "omega_val", "omega_m_val", "kappa_val", "rho_val", 
               "epsilon_val", "mu_val", "mj_val", "m_val",  "delta_t"
               )

# can go by 1, as we specify this is delta_t explicitly
pomp_object <- pomp(data=data.frame(time=seq(0,1,by=1/365),cases=NA),
                 times="time",t0=0,
                 skeleton=map(det_model_skeleton, delta.t = 1/365),
                 rinit=init_states,
                 statenames=state_names,
                 paramnames=param_names
)

x <- trajectory(pomp_object,params=params,format="d")

library(ggplot2)
ggplot(data=x,mapping=aes(x=time,y=If))+geom_line()
ggplot(data=x,mapping=aes(x=time,y=Rf))+geom_line()
ggplot(data=x,mapping=aes(x=time,y=Sf))+geom_line()
