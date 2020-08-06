/* pomp C snippet file: hacking_win_bug */
/* Time: 2020-08-02 21:56:33.779 +0200 */
/* Salt: 71774ECB5E8B8A6CC7F6247A */

#include <C:/Users/alexa/OneDrive/Documents/R/win-library/4.0/pomp/include/pomp.h>
#include <R_ext/Rdynload.h>

 


/* C snippet: 'rinit' */
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define samplesize		(__covars[__covindex[0]])
#define Ma		(__x[__stateindex[0]])
#define Sn		(__x[__stateindex[1]])
#define Sj		(__x[__stateindex[2]])
#define Sm		(__x[__stateindex[3]])
#define Sf		(__x[__stateindex[4]])
#define En		(__x[__stateindex[5]])
#define Ej		(__x[__stateindex[6]])
#define Em		(__x[__stateindex[7]])
#define Ef		(__x[__stateindex[8]])
#define In		(__x[__stateindex[9]])
#define Ij		(__x[__stateindex[10]])
#define Im		(__x[__stateindex[11]])
#define If		(__x[__stateindex[12]])
#define Rn		(__x[__stateindex[13]])
#define Rj		(__x[__stateindex[14]])
#define Rm		(__x[__stateindex[15]])
#define Rf		(__x[__stateindex[16]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 

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

                         
}

#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef samplesize
#undef Ma
#undef Sn
#undef Sj
#undef Sm
#undef Sf
#undef En
#undef Ej
#undef Em
#undef Ef
#undef In
#undef Ij
#undef Im
#undef If
#undef Rn
#undef Rj
#undef Rm
#undef Rf

/* C snippet: 'step.fn' */
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define samplesize		(__covars[__covindex[0]])
#define Ma		(__x[__stateindex[0]])
#define Sn		(__x[__stateindex[1]])
#define Sj		(__x[__stateindex[2]])
#define Sm		(__x[__stateindex[3]])
#define Sf		(__x[__stateindex[4]])
#define En		(__x[__stateindex[5]])
#define Ej		(__x[__stateindex[6]])
#define Em		(__x[__stateindex[7]])
#define Ef		(__x[__stateindex[8]])
#define In		(__x[__stateindex[9]])
#define Ij		(__x[__stateindex[10]])
#define Im		(__x[__stateindex[11]])
#define If		(__x[__stateindex[12]])
#define Rn		(__x[__stateindex[13]])
#define Rj		(__x[__stateindex[14]])
#define Rm		(__x[__stateindex[15]])
#define Rf		(__x[__stateindex[16]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 
  /* -------------------------- seasonal forcing ----------------------------------------- */
  
  double b = c * exp(-s * (pow(cos(M_PI * t / 365 - phi), 2)));
  
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
                        
   
}

#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef samplesize
#undef Ma
#undef Sn
#undef Sj
#undef Sm
#undef Sf
#undef En
#undef Ej
#undef Em
#undef Ef
#undef In
#undef Ij
#undef Im
#undef If
#undef Rn
#undef Rj
#undef Rm
#undef Rf

/* C snippet: 'rmeasure' */
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define samplesize		(__covars[__covindex[0]])
#define Ma		(__x[__stateindex[0]])
#define Sn		(__x[__stateindex[1]])
#define Sj		(__x[__stateindex[2]])
#define Sm		(__x[__stateindex[3]])
#define Sf		(__x[__stateindex[4]])
#define En		(__x[__stateindex[5]])
#define Ej		(__x[__stateindex[6]])
#define Em		(__x[__stateindex[7]])
#define Ef		(__x[__stateindex[8]])
#define In		(__x[__stateindex[9]])
#define Ij		(__x[__stateindex[10]])
#define Im		(__x[__stateindex[11]])
#define If		(__x[__stateindex[12]])
#define Rn		(__x[__stateindex[13]])
#define Rj		(__x[__stateindex[14]])
#define Rm		(__x[__stateindex[15]])
#define Rf		(__x[__stateindex[16]])
#define pos		(__y[__obsindex[0]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
                          
                          double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If 
                                      + Rn + Rj + Rm + Rf + Ma);
                          double I_total =  In + Ij + Im + If;
                          double p = zeta * I_total / N;
                          double p_pool = 1 - pow((1 - p), d);
                          pos = rbetabinom(samplesize, p_pool, disp);
                           
}

#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef samplesize
#undef Ma
#undef Sn
#undef Sj
#undef Sm
#undef Sf
#undef En
#undef Ej
#undef Em
#undef Ef
#undef In
#undef Ij
#undef Im
#undef If
#undef Rn
#undef Rj
#undef Rm
#undef Rf
#undef pos

/* C snippet: 'dmeasure' */
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define samplesize		(__covars[__covindex[0]])
#define Ma		(__x[__stateindex[0]])
#define Sn		(__x[__stateindex[1]])
#define Sj		(__x[__stateindex[2]])
#define Sm		(__x[__stateindex[3]])
#define Sf		(__x[__stateindex[4]])
#define En		(__x[__stateindex[5]])
#define Ej		(__x[__stateindex[6]])
#define Em		(__x[__stateindex[7]])
#define Ef		(__x[__stateindex[8]])
#define In		(__x[__stateindex[9]])
#define Ij		(__x[__stateindex[10]])
#define Im		(__x[__stateindex[11]])
#define If		(__x[__stateindex[12]])
#define Rn		(__x[__stateindex[13]])
#define Rj		(__x[__stateindex[14]])
#define Rm		(__x[__stateindex[15]])
#define Rf		(__x[__stateindex[16]])
#define pos		(__y[__obsindex[0]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
                          
                          double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If 
                                      + Rn + Rj + Rm + Rf + Ma);
                          double I_total =  In + Ij + Im + If;
                          double p = zeta * I_total / N;
                          double p_pool = 1 - pow((1 - p), d);
                          lik = dbetabinom(pos, samplesize, p_pool+0.00001, 
                                           disp, give_log);
                           
}

#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef samplesize
#undef Ma
#undef Sn
#undef Sj
#undef Sm
#undef Sf
#undef En
#undef Ej
#undef Em
#undef Ef
#undef In
#undef Ij
#undef Im
#undef If
#undef Rn
#undef Rj
#undef Rm
#undef Rf
#undef pos
#undef lik

/* C snippet: 'toEst' */
#define samplesize		(__covars[__covindex[0]])
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define T_model_type		(__pt[__parindex[0]])
#define T_R0		(__pt[__parindex[1]])
#define T_gamma_val		(__pt[__parindex[2]])
#define T_omega_val		(__pt[__parindex[3]])
#define T_omega_m_val		(__pt[__parindex[4]])
#define T_kappa_val		(__pt[__parindex[5]])
#define T_rho_val		(__pt[__parindex[6]])
#define T_epsilon_val		(__pt[__parindex[7]])
#define T_mu_val		(__pt[__parindex[8]])
#define T_mj_val		(__pt[__parindex[9]])
#define T_m_val		(__pt[__parindex[10]])
#define T_delta_t		(__pt[__parindex[11]])
#define T_c		(__pt[__parindex[12]])
#define T_s		(__pt[__parindex[13]])
#define T_phi		(__pt[__parindex[14]])
#define T_c_v		(__pt[__parindex[15]])
#define T_s_v		(__pt[__parindex[16]])
#define T_phi_v		(__pt[__parindex[17]])
#define T_zeta		(__pt[__parindex[18]])
#define T_disp		(__pt[__parindex[19]])
#define T_d		(__pt[__parindex[20]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_R0 = log(R0);
	T_c = log(c);
	T_s = log(s);
	T_phi = log(phi);
	T_disp = log(disp);
	T_d = log(d);
	T_gamma_val = log(gamma_val);
	T_omega_m_val = log(omega_m_val);
	T_s_v = log(s_v);
	T_phi_v = log(phi_v);
	T_rho_val = log(rho_val);
	T_epsilon_val = log(epsilon_val);
	T_zeta = logit(zeta); 
}

#undef samplesize
#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef T_model_type
#undef T_R0
#undef T_gamma_val
#undef T_omega_val
#undef T_omega_m_val
#undef T_kappa_val
#undef T_rho_val
#undef T_epsilon_val
#undef T_mu_val
#undef T_mj_val
#undef T_m_val
#undef T_delta_t
#undef T_c
#undef T_s
#undef T_phi
#undef T_c_v
#undef T_s_v
#undef T_phi_v
#undef T_zeta
#undef T_disp
#undef T_d

/* C snippet: 'fromEst' */
#define samplesize		(__covars[__covindex[0]])
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define T_model_type		(__pt[__parindex[0]])
#define T_R0		(__pt[__parindex[1]])
#define T_gamma_val		(__pt[__parindex[2]])
#define T_omega_val		(__pt[__parindex[3]])
#define T_omega_m_val		(__pt[__parindex[4]])
#define T_kappa_val		(__pt[__parindex[5]])
#define T_rho_val		(__pt[__parindex[6]])
#define T_epsilon_val		(__pt[__parindex[7]])
#define T_mu_val		(__pt[__parindex[8]])
#define T_mj_val		(__pt[__parindex[9]])
#define T_m_val		(__pt[__parindex[10]])
#define T_delta_t		(__pt[__parindex[11]])
#define T_c		(__pt[__parindex[12]])
#define T_s		(__pt[__parindex[13]])
#define T_phi		(__pt[__parindex[14]])
#define T_c_v		(__pt[__parindex[15]])
#define T_s_v		(__pt[__parindex[16]])
#define T_phi_v		(__pt[__parindex[17]])
#define T_zeta		(__pt[__parindex[18]])
#define T_disp		(__pt[__parindex[19]])
#define T_d		(__pt[__parindex[20]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	R0 = exp(T_R0);
	c = exp(T_c);
	s = exp(T_s);
	phi = exp(T_phi);
	disp = exp(T_disp);
	d = exp(T_d);
	gamma_val = exp(T_gamma_val);
	omega_m_val = exp(T_omega_m_val);
	s_v = exp(T_s_v);
	phi_v = exp(T_phi_v);
	rho_val = exp(T_rho_val);
	epsilon_val = exp(T_epsilon_val);
	zeta = expit(T_zeta); 
}

#undef samplesize
#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef T_model_type
#undef T_R0
#undef T_gamma_val
#undef T_omega_val
#undef T_omega_m_val
#undef T_kappa_val
#undef T_rho_val
#undef T_epsilon_val
#undef T_mu_val
#undef T_mj_val
#undef T_m_val
#undef T_delta_t
#undef T_c
#undef T_s
#undef T_phi
#undef T_c_v
#undef T_s_v
#undef T_phi_v
#undef T_zeta
#undef T_disp
#undef T_d

/* C snippet: 'skeleton' */
#define model_type		(__p[__parindex[0]])
#define R0		(__p[__parindex[1]])
#define gamma_val		(__p[__parindex[2]])
#define omega_val		(__p[__parindex[3]])
#define omega_m_val		(__p[__parindex[4]])
#define kappa_val		(__p[__parindex[5]])
#define rho_val		(__p[__parindex[6]])
#define epsilon_val		(__p[__parindex[7]])
#define mu_val		(__p[__parindex[8]])
#define mj_val		(__p[__parindex[9]])
#define m_val		(__p[__parindex[10]])
#define delta_t		(__p[__parindex[11]])
#define c		(__p[__parindex[12]])
#define s		(__p[__parindex[13]])
#define phi		(__p[__parindex[14]])
#define c_v		(__p[__parindex[15]])
#define s_v		(__p[__parindex[16]])
#define phi_v		(__p[__parindex[17]])
#define zeta		(__p[__parindex[18]])
#define disp		(__p[__parindex[19]])
#define d		(__p[__parindex[20]])
#define samplesize		(__covars[__covindex[0]])
#define Ma		(__x[__stateindex[0]])
#define Sn		(__x[__stateindex[1]])
#define Sj		(__x[__stateindex[2]])
#define Sm		(__x[__stateindex[3]])
#define Sf		(__x[__stateindex[4]])
#define En		(__x[__stateindex[5]])
#define Ej		(__x[__stateindex[6]])
#define Em		(__x[__stateindex[7]])
#define Ef		(__x[__stateindex[8]])
#define In		(__x[__stateindex[9]])
#define Ij		(__x[__stateindex[10]])
#define Im		(__x[__stateindex[11]])
#define If		(__x[__stateindex[12]])
#define Rn		(__x[__stateindex[13]])
#define Rj		(__x[__stateindex[14]])
#define Rm		(__x[__stateindex[15]])
#define Rf		(__x[__stateindex[16]])
#define DMa		(__f[__stateindex[0]])
#define DSn		(__f[__stateindex[1]])
#define DSj		(__f[__stateindex[2]])
#define DSm		(__f[__stateindex[3]])
#define DSf		(__f[__stateindex[4]])
#define DEn		(__f[__stateindex[5]])
#define DEj		(__f[__stateindex[6]])
#define DEm		(__f[__stateindex[7]])
#define DEf		(__f[__stateindex[8]])
#define DIn		(__f[__stateindex[9]])
#define DIj		(__f[__stateindex[10]])
#define DIm		(__f[__stateindex[11]])
#define DIf		(__f[__stateindex[12]])
#define DRn		(__f[__stateindex[13]])
#define DRj		(__f[__stateindex[14]])
#define DRm		(__f[__stateindex[15]])
#define DRf		(__f[__stateindex[16]])

void __pomp_skelfn (double *__f, const double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
                      
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
                          
                          
                           
}

#undef model_type
#undef R0
#undef gamma_val
#undef omega_val
#undef omega_m_val
#undef kappa_val
#undef rho_val
#undef epsilon_val
#undef mu_val
#undef mj_val
#undef m_val
#undef delta_t
#undef c
#undef s
#undef phi
#undef c_v
#undef s_v
#undef phi_v
#undef zeta
#undef disp
#undef d
#undef samplesize
#undef Ma
#undef Sn
#undef Sj
#undef Sm
#undef Sf
#undef En
#undef Ej
#undef Em
#undef Ef
#undef In
#undef Ij
#undef Im
#undef If
#undef Rn
#undef Rj
#undef Rm
#undef Rf
#undef DMa
#undef DSn
#undef DSj
#undef DSm
#undef DSf
#undef DEn
#undef DEj
#undef DEm
#undef DEf
#undef DIn
#undef DIj
#undef DIm
#undef DIf
#undef DRn
#undef DRj
#undef DRm
#undef DRf

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_hacking_win_bug (DllInfo *info)
{
R_RegisterCCallable("hacking_win_bug", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("hacking_win_bug", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("hacking_win_bug", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("hacking_win_bug", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
R_RegisterCCallable("hacking_win_bug", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("hacking_win_bug", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("hacking_win_bug", "__pomp_to_trans", (DL_FUNC) __pomp_to_trans);
R_RegisterCCallable("hacking_win_bug", "__pomp_from_trans", (DL_FUNC) __pomp_from_trans);
R_RegisterCCallable("hacking_win_bug", "__pomp_skelfn", (DL_FUNC) __pomp_skelfn);
}
