# model 8 SIRS likelihood profile for gamma
# phi is fixed at 7.179
# we saw that gamma looks like it's converging so it makes sense to do a likelihood profile to get the 95% confidence interval.

setwd("~/to_server_files")

library(foreach)
library(iterators)
library(plyr)
library(dplyr)
library(doParallel)
library(doRNG)
library(pomp)

detectCores() -> cores
registerDoParallel(cores= cores-1)
registerDoRNG(625904618)

source("full_pomp_model.R", local = TRUE)
guesses <- read.csv("guesses.csv", row.names=1)


  
# ----------------------------------------- params of the filtering -------------------------------
Nmif = 50
Np = 2500
# ------------------------------------------- filtering ---------------------------------------------
mf3<- foreach(guess=iter(guesses,"row"),
              .combine=c,
              .packages = "pomp"
              
) %dopar%
  
  {   mif2(
    pomp_object,
    params = guess,
    Nmif = Nmif,
    Np = Np,
    cooling.fraction.50=0.5,
    cooling.type="geometric",
    rw.sd=rw.sd(R0=0.02, zeta=0.02, omega_val = 0.02,
                c=0.02, s=0.02, s_v=0.02, phi_v=0.02, disp=0.02, d=0.01)
  )
  }

save(guesses, mf3,  file = "model8_lik_filtering_cont1.rda")
# ------------------------------------------- filtering part2 ---------------------------------------------
mf3iter2<- foreach(mf3item=iter(mf3),
                   .combine=c,
                   .packages = "pomp"
                   
) %dopar%
  {
    continue(mf3item, Nmif = 50)
  }

save(guesses, mf3iter2,  file = "model8_lik_filtering_cont2.rda")

# -------------- calculating the likelihoods properly ------------------------------------------------

lik_list <- foreach(mf3item=iter(mf3iter2),
                    .combine=rbind,
                    .packages = "pomp"
                    
) %dopar%
  {
    evals <- replicate(5,logLik(pfilter(mf3item,Np=4000)))
    ll <- logmeanexp(evals,se=TRUE)
    c(coef(mf3item),loglik=ll[1],loglik=ll[2])
  }

 lik_list <- as.data.frame(lik_list)
# # #staring a new csv file 
lik_list %>%
  arrange(-loglik) %>%
  write.csv("model8_server_likelihoods.csv")




