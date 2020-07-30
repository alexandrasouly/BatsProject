source("full_pomp_model.R", local = TRUE)

library(foreach)
library(iterators)
library(dplyr)
library(ggplot2)
library(ggpubr)

library(doParallel)
registerDoParallel(cores= 6)
library(doRNG)
registerDoRNG(625904618)
library(pomp)

sobolDesign(
  lower=model_4_params_lower,
  upper=model_4_params_upper,
  nseq=20
) -> guesses

Nmif = 1
Np = 2000

      mf3<- foreach(guess=iter(guesses,"row"),
                .combine=c,
                .packages = "pomp"

      ) %dopar%
      {
        mif2(
          pomp_object,
          params = guess,
          Nmif = Nmif,
          Np = Np,
          cooling.fraction.50=0.5,
          cooling.type="geometric",
          rw.sd=rw.sd(R0=0.02, rho_val=0.02, epsilon_val = 0.02, zeta=0.02,
                      c=0.02, s=0.02, phi=0.02, s_v=0.02, phi_v=0.02, disp=0.02, d=0.01)
        )
      }

save(Nmif, Np, guesses,  file = "Sili_filtering_1it.rda")

mf3iter20<- foreach(mf3item=iter(mf3),
                  .combine=c,
                  .packages = "pomp"
                  
            ) %dopar%
  {
              continue(mf3item, Nmif = 2)}

mf3 %>%
  traces() %>%
  melt() %>%
  filter(variable %in% c("loglik", "R0", "rho_val", "epsilon_val",  "zeta", "c",
         "s", "phi","s_v", "phi_v", "disp", "d") )%>%
  ggplot(aes(x=iteration,y=value,group=L1,color=L1))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  guides(color=FALSE)

