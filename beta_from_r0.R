
ini_pop_size <-function(m, mj, mu, kappa, omega_m){
  b <- (omega_m+mj)*((mj+mu)/omega_m)*m/mu
  N<-kappa
  Sm=(round(N/(1+m/mu+b/(omega_m+mj)))/2)
  Sf=(round(N/(1+m/mu+b/(omega_m+mj)))/2)
  Na <- Sf + Sm
  Sj <- round(Na*m/mu)
  Sn <- N - Na - Sj
  If <- round(N*0.025)
  Im <- round(N*0.025)
  return(list("N"=N, "Sf" = Sf, "Sm"=Sm, "Sj"=Sj, "Sn"=Sn, "If" = If, "Im" =Im))
}

pop_size <-ini_pop_size( model_4_params[["m_val"]], model_4_params[["mj_val"]] ,
                         model_4_params[["mu_val"]], model_4_params[["kappa_val"]], model_4_params[["omega_m_val"]])
pop_size

pop_size <-ini_pop_size( 0.187, 0.5 ,
                         1.37, 5000, 0.799)
pop_size