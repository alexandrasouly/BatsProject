library(pomp)
stochStep <- Csnippet("

  /* -------------------------- seasonal forcing of birth ----------------------------------------- */
  
  double b = c * exp( - s * (pow(cos(M_PI * t - phi), 2)));

  /*---------------------------- pre-calculations and scaling params:------------------------------- */
  
  double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If 
    + Rn + Rj + Rm + Rf + Ma);
    
  double beta = beta_val / delta_t;
  double gamma = gamma_val / delta_t;
  double omega = omega_val / delta_t;
  double omega_m = omega_m_val / delta_t;
  
  double kappa = kappa_val;
  
  double rho = rho_val / delta_t;
  double epsilon = epsilon_val / delta_t;
  double mu = mu_val / delta_t;
  double mj = mj_val / delta_t;
  double m = m_val /delta_t;
  
  
  /* ------------------------ general rates of changes  ---------------------------------------------- */
  
  double r_births = b / delta_t; 
  double n_Sbirths = rbinom(Sf + If, 1 - exp( - r_births));
  double n_MaBirths = rbinom(Rf + Ef, 1 - exp( - r_births));
  
  double r_juv_death = mj * (N / kappa);
  double r_death = m * (N / kappa);
  
  double r_N_J = omega_m;
  double r_J_A = mu;
  
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
  double r_In_out = r_juv_death + r_N_J + r_I_R +r_I_E;
  double r_Rn_out = r_juv_death + r_N_J + r_R_S;
  
  /* numbers changing state */
  
  double n_Ma_out = rbinom(Ma, 1 - exp( - r_Ma_out));
  double n_Sn_out = rbinom(Sn, 1 - exp( - r_Sn_out));
  double n_En_out = rbinom(En, 1 - exp( - r_En_out));
  double n_In_out = rbinom(In, 1 - exp( - r_In_out));
  double n_Rn_out = rbinom(Rn, 1 - exp( - r_Rn_out));
  
  double probs_Ma_out[2] = {r_juv_death/r_Ma_out, r_N_J/r_Ma_out};
  double probs_Sn_out[3] = {r_juv_death/r_Sn_out, r_N_J/r_Sn_out, r_S_I/r_Sn_out};
  double probs_En_out[3] = {r_juv_death/r_En_out, r_N_J/r_En_out, r_E_I/r_En_out};
  double probs_In_out[4] = {r_juv_death/r_In_out, r_N_J/r_In_out, r_I_R/r_In_out, r_I_E/r_In_out};
  double probs_Rn_out[3] = {r_juv_death/r_Rn_out, r_N_J/r_Rn_out, r_R_S/r_Rn_out};
  
  int multi_Ma_out[2];
  int multi_Sn_out[3];
  int multi_En_out[3];
  int multi_In_out[4];
  int multi_Rn_out[3];
  
  rmultinom(n_Ma_out, probs_Ma_out, 2, &multi_Ma_out[0]);
  rmultinom(n_Sn_out, probs_Sn_out, 3, &multi_Sn_out[0]);
  rmultinom(n_En_out, probs_En_out, 3, &multi_En_out[0]);
  rmultinom(n_In_out, probs_In_out, 4, &multi_In_out[0]);
  rmultinom(n_Rn_out, probs_Rn_out, 3, &multi_Rn_out[0]);
  
  double n_Rn_Sn = multi_Rn_out[2];
  double n_In_En = multi_In_out[3];
  double n_Sn_In = multi_Sn_out[2];
  double n_En_In = multi_En_out[2];
  double n_In_Rn = multi_In_out[2];
  
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
  double r_Ij_out = r_juv_death + r_J_A + r_I_R +r_I_E;
  double r_Rj_out = r_juv_death + r_J_A + r_R_S;
  
  /* numbers changing state */
  
  double n_Sj_out = rbinom(Sj, 1 - exp( - r_Sj_out));
  double n_Ej_out = rbinom(Ej, 1 - exp( - r_Ej_out));
  double n_Ij_out = rbinom(Ij, 1 - exp( - r_Ij_out));
  double n_Rj_out = rbinom(Rj, 1 - exp( - r_Rj_out));
  
  double probs_Sj_out[3] = {r_juv_death/r_Sj_out, r_J_A/r_Sj_out, r_S_I/r_Sj_out};
  double probs_Ej_out[3] = {r_juv_death/r_Ej_out, r_J_A/r_Ej_out, r_E_I/r_Ej_out};
  double probs_Ij_out[4] = {r_juv_death/r_Ij_out, r_J_A/r_Ij_out, r_I_R/r_Ij_out, r_I_E/r_Ij_out};
  double probs_Rj_out[3] = {r_juv_death/r_Rj_out, r_J_A/r_Rj_out, r_R_S/r_Rj_out};
  
  int multi_Sj_out[3];
  int multi_Ej_out[3];
  int multi_Ij_out[4];
  int multi_Rj_out[3];
  
  rmultinom(n_Sj_out, probs_Sj_out, 3, &multi_Sj_out[0]);
  rmultinom(n_Ej_out, probs_Ej_out, 3, &multi_Ej_out[0]);
  rmultinom(n_Ij_out, probs_Ij_out, 4, &multi_Ij_out[0]);
  rmultinom(n_Rj_out, probs_Rj_out, 3, &multi_Rj_out[0]);
  
  double n_Rj_Sj = multi_Rj_out[2];
  double n_Ij_Ej = multi_Ij_out[3];
  double n_Sj_Ij = multi_Sj_out[2];
  double n_Ej_Ij = multi_Ej_out[2];
  double n_Ij_Rj = multi_Ij_out[2];
  
  /* equations of changing state */
  
  THIS HEEEEREEEEE Add NEWBORNS into all of this
  if (Sj  - n_Sj_out + n_Rj_Sj <= 0) {
    Sj = 0;
  } else {
    Sj = Sj - n_Sj_out + n_Rj_Sj;
  }
  
  if (Ej - n_Ej_out + n_Ij_Ej <= 0) {
    Ej = 0;
  } else {
    Ej = Ej - n_Ej_out + n_Ij_Ej;
  }
  
  if (Ij - n_Ij_out + n_Sj_Ij + n_Ej_Ij <= 0) {
    Ij = 0;
  } else {
    Ij = Ij - n_Ij_out + n_Sj_Ij + n_Ej_Ij;
  }
  
  if (Rj - n_Rj_out + n_Ij_Rj <= 0) {
    Rj = 0;
  } else {
    Rj = Rj - n_Rj_out + n_Ij_Rj;
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

params <- c(
  
  beta_val = 0.02, # infection rate S -> R
  gamma_val = 0.8, # recovery rate I -> R
  omega_val = 0, # immune waning rate R -> S
  omega_m_val = 0.4, # maternal antibody waning rate
  kappa_val = 1000, # carrying capacity
  rho_val = 0, # I -> L 
  epsilon_val = 0, # L -> I 
  mu_val = 0.44, # juvenile maturation rate
  mj_val = 0.8, # juvenile death rate
  m_val = 0.2, # adult death rate
  delta_t = 365, # scaling time as days instead of years
  c = 0, #birth pulse scaling factor
  s = 14.3, #birth pulse synchronicity
  phi = 4.5) #birth pulse time shift

state_names <- c("Ma",
                 "Sn","Sj", "Sm", "Sf",
                 "En","Ej", "Em", "Ef",
                 "In","Ij", "Im", "If",
                 "Rn","Rj", "Rm", "Rf"
)


param_names <- c( "beta_val", "gamma_val", "omega_val", "omega_m_val", "kappa_val", "rho_val", 
                  "epsilon_val", "mu_val", "mj_val", "m_val",  "delta_t",  "c", "s", "phi"
)

pomp_object <- pomp(data=data.frame(time=seq(0,3,by=1/365),cases=NA),
                    times="time",t0=0,
                    rprocess=discrete_time(step.fun=stochStep,delta.t=1/365),
                    rinit=init_states,
                    statenames=state_names,
                    paramnames=param_names
)

sim <- simulate(pomp_object,params=params,format = "data.frame", nsim = 500)

created_plot <-ggplot(NULL)+
  geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sn, group= factor(.id), colour = "Stochastic infected"))+
  labs(title="Stochastic SIR model")