source("full_pomp_model.R")

library(foreach)
library(iterators)
library(dplyr)
library(ggplot2)
library(pomp)
library(ggpubr)
ggarrange(plotlist = plots, ncol = 4)
plprev

pf <- replicate(5, pfilter(pomp_object, params = params, Np=2000))

ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

plot(pf[[1]])


library(doParallel)
registerDoParallel(cores= 8)
library(doRNG)
registerDoRNG(625904618)
library(pomp)

system.time(
  foreach(i = 1, .combine=c, .packages = "pomp") %dopar%
    {print(i)
      mif2( pomp_object, params = model_4_params, Nmif=20, Np=500, 
            cooling.fraction.50=0.5,cooling.type="geometric",
            rw.sd=rw.sd(R0=0.02, rho_val=0.02, epsilon_val = 0.02, 
                        zeta=0.02, c=0.02,
                        s=0.02, phi=0.02,
                        s_v=0.02, phi_v=0.02, disp=0.02, d=0.01)
      ) }-> mf3
)

# evaluate the log likelihood here
foreach(start = iter(mf3), .combine=rbind) %dopar% 
  {replicate(5, start %>% pfilter() %>% logLik()) %>%
    logmeanexp(se=TRUE) -> ll} -> ll_list

rbind( coef(mf3[[1]]) , coef(mf3[[2]]), coef(mf3[[3]]), coef(mf3[[4]]), coef(mf3[[5]]),
coef(mf3[[6]]), coef(mf3[[7]]), coef(mf3[[8]]), coef(mf3[[9]]), coef(mf3[[10]])
)

# try:
#mf3_coef <- rbind(lapply(mf3, coef))

data.frame(mf3_coef,loglik=ll_list[,1],loglik.se=ll_list[,2]) -> estimates3

estimates3  %>%
  select(loglik, loglik.se, R0, rho_val, epsilon_val,  zeta, c,
         s, phi,s_v, phi_v, disp, d) %>%
  arrange(-loglik)

mf3 %>%
  traces() %>%
  melt() %>%
  filter( variable %in% c("loglik","R0", "gamma_val",  "zeta", "c",
                          "s", "phi", "disp", "d")) -> filtered_plots

ggplot(filtered_plots)+
  aes(x=iteration,y=value,group=L1,color=factor(L1))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")

