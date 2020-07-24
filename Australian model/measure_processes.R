measure_processes <- function(){
  dmeas <- Csnippet("
                          
                          double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If 
                                      + Rn + Rj + Rm + Rf + Ma);
                          double I_total =  In + Ij + Im + If;
                          double p = zeta * I_total / N;
                          double p_pool = 1 - pow((1 - p), d);
                          lik = dbetabinom(pos, samplesize, p_pool+0.00001, 
                                           disp, give_log);
                          ")
  
  rmeas <- Csnippet("
                          
                          double N = (Sn + Sj + Sm + Sf + En + Ej + Em + Ef + In + Ij + Im + If 
                                      + Rn + Rj + Rm + Rf + Ma);
                          double I_total =  In + Ij + Im + If;
                          double p = zeta * I_total / N;
                          double p_pool = 1 - pow((1 - p), d);
                          pos = rbetabinom(samplesize, p_pool, disp);
                          ")                          
  return(list(dmeas, rmeas))
}