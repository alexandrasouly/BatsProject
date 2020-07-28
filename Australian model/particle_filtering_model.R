source("full_pomp_model.R")

library(foreach)
library(iterators)
library(dplyr)
library(ggplot2)
library(pomp)

pf <- replicate(10, pfilter(pomp_object, params = params, Np=5000))

ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

plot(pf[[1]])

