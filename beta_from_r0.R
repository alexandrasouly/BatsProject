beta_from_r0 <-function(R0, m, gamma, kappa, epsilon, rho){
  
  beta = R0 *((epsilon + m)*(gamma + m + rho)-epsilon*rho) / (kappa*(epsilon + m))
}


beta<-beta_from_r0(42.152, 0.23, 3.135, 2350, 0, 0)


ini_pop_size <-function(m, mj, mu, kappa, omega_m){
  b <- (omega_m+mj)*((mj+mu)/omega_m)*m/mu
  N<-kappa
  Na=(round(N/(1+m/mu+b/(omega_m+mj))))
  Nj <- round(Na*m/mu)
  Nn <- N - Na - Nj
  params <- c(N, Na, Nj, Nn)
}

params <-ini_pop_size(0.23, 0.5, 0.44, 2350, 0.8)
params

