library(pomp)

closed_sir_ode <- Csnippet("
  DS = -Beta*S*I/N;
  DI = Beta*S*I/N-gamma*I;
  DR = gamma*I;
")

init1 <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  ")

closed_sir <- pomp(data=data.frame(time=1:50,data=NA),
     times="time",t0=0,
     skeleton=vectorfield(closed_sir_ode),
     rinit=init1,
     statenames=c("S","I","R"),
     paramnames=c("Beta","gamma","N"))

params1 <- c(Beta=1,gamma=1/13,N=763)

x <- trajectory(closed_sir,params=params1,format="data.frame")

library(ggplot2)
ggplot(data=x,mapping=aes(x=time,y=I))+geom_line()