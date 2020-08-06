source("full_pomp_model.R", local = TRUE)
source("aarons_params.R", local = TRUE)
library(foreach)
library(iterators)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(doParallel)
registerDoParallel(cores= 6)
library(doRNG)
registerDoRNG(625904618)
library(pomp)
# ----------------------------------------- params of the filtering --------------------------------
sobolDesign(
  lower= model_8_params_lower,
  upper= model_8_params_upper,
  nseq=24) -> guesses

guesses <-  rbind(guesses, model_8_params, model_8_params, model_8_params, model_8_params, model_8_params)

Nmif = 20
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
          rw.sd=rw.sd(R0=0.02, zeta=0.02, gamma_val = 0.02, omega_val = 0.02,
                      c=0.02, s=0.02, s_v=0.02, phi_v=0.02, disp=0.02, d=0.01)
        )
      }

      save(guesses, mf3,  file = "model8_filtering_pt1.rda")
# ------------------------------------------- filtering part2 ---------------------------------------------
      mf3iter2<- foreach(mf3item=iter(mf3),
                        .combine=c,
                        .packages = "pomp"
      
                  ) %dopar%
                  {
                    continue(mf3item, Nmif = 30)
                    }

      save(guesses, mf3iter2,  file = "model8_filtering_pt2.rda")

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
      # #staring a new csv file 
      lik_list %>%
      arrange(-loglik) %>%
      write.csv("model8_likelihoods.csv")
      # 
      # 
      # # #adding in more points to my csv of likelihoods
      # read.table("model4_likelihoods.csv",sep = ",", row.names = 1, header = TRUE) %>%
      # bind_rows(lik_list) %>%
      # arrange(-loglik) %>% write.csv("model4_likelihoods.csv")



