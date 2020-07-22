beta_from_r0 <-function(R0, m, gamma, kappa, epsilon, rho){
  
  beta = R0 *((epsilon + m)*(gamma + m + rho)-epsilon*rho) / (kappa*(epsilon + m))
}


beta<-beta_from_r0(13.395, 0.187, 22.449, 4216, 0, 0)


ini_pop_size <-function(m, mj, mu, kappa, omega_m){
  b <- (omega_m+mj)*((mj+mu)/omega_m)*m/mu
  N<-kappa
  Na=(round(N/(1+m/mu+b/(omega_m+mj))))
  Nj <- round(Na*m/mu)
  Nn <- N - Na - Nj
  params <- c(N, Na, Nj, Nn)
}

pop_size <-ini_pop_size( model_6_params[["m_val"]], model_6_params[["mj_val"]] ,
                         model_6_params[["mu_val"]], model_6_params[["kappa_val"]], model_6_params[["omega_m_val"]])
pop_size

