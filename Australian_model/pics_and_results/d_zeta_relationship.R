# lower <- c(zeta = 0.2, d = 1)
# upper <- c(zeta = 1, d = 20)

library(foreach)
library(iterators)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(doParallel)
registerDoParallel(cores= 4)
library(doRNG)
registerDoRNG(625904618)
library(pomp)


zeta_d_guesses <- sobolDesign(lower =lower, upper = upper, nseq = 100)
plot(zeta_d_guesses)

head(model8_likelihoods, 5) %>%
  select(-zeta, -d, -loglik, - loglik.se) -> other_params_guesses

guesses <- vector(mode = "list", length = 5)
for(i in 1:5)
{merge(zeta_d_guesses, other_params_guesses[i,]) -> guesses[[i]]}

#for now only doing best param set
source("pomp_model_fn.R", local = TRUE)
source("aarons_params.R", local = TRUE)
lik_list2 <- foreach(guess=iter(guesses[[2]],by = "row" ),
                    .combine=rbind,
                    .packages = "pomp"
                    
) %dopar%
  { pomp_model_fn("CLU", "all", model_8_before_equ_ini, as.list(guess)) -> pomp_object
    evals <- replicate(5,logLik(pfilter(pomp_object, params = guess, Np=250)))
    ll <- logmeanexp(evals,se=TRUE)
    c(guess,loglik=ll[1],loglik=ll[2])
  }

lik_list2 <- as.data.frame(lik_list2)


lik_list3 <- foreach(guess=iter(guesses[[3]],by = "row" ),
                     .combine=rbind,
                     .packages = "pomp"
                     
) %dopar%
  { pomp_model_fn("CLU", "all", model_8_before_equ_ini, as.list(guess)) -> pomp_object
    evals <- replicate(5,logLik(pfilter(pomp_object, params = guess, Np=250)))
    ll <- logmeanexp(evals,se=TRUE)
    c(guess,loglik=ll[1],loglik=ll[2])
  }

lik_list2 <- as.data.frame(lik_list3)


