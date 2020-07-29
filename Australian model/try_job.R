source("full_pomp_model.R", local = TRUE)

library(foreach)
library(iterators)
library(dplyr)
library(ggplot2)
library(ggpubr)

library(doParallel)
registerDoParallel(cores= 8)
library(doRNG)
registerDoRNG(625904618)
library(pomp)



Nmif = 10 
Np = 30
params = model_4_params

# filtering---
  run_time <- system.time({
        mf3<- foreach(i = 1:2,
                  .combine=c, 
                  .packages = "pomp"
        ) %dopar%
        {
          cat(i)
          mif2( 
            pomp_object,
            params = params,
            Nmif = Nmif, 
            Np = Np, 
            cooling.fraction.50=0.5,
            cooling.type="geometric",
            rw.sd=rw.sd(R0=0.02, rho_val=0.02, epsilon_val = 0.02, zeta=0.02,
                        c=0.02, s=0.02, phi=0.02, s_v=0.02, phi_v=0.02, disp=0.02, d=0.01)
          )
        }
  })
  
  save(Nmif, Np, run_time, params,  file = "trying_to_save_data.rda")

