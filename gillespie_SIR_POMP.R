## Gillespie algorithm taken from a POMP tutorial
## SIR model with biths and deaths

library(pomp)

pomp(
  data=data.frame(
    time=seq(0.01,0.5,by=0.01),
    reports=NA
  ),
  times="time",t0=0,
  rprocess=gillespie_hl(
    birth=list("rate = mu*N;",c(N=1,X=1,Y=0,Z=0,cases=0)),
    deathS=list("rate=mu*X;",c(N=-1,X=-1,Y=0,Z=0,cases=0)),
    deathI=list("rate=mu*Y;",c(N=-1,X=0,Y=-1,Z=0,cases=0)),
    deathR=list("rate=mu*Z;",c(N=-1,X=0,Y=0,Z=-1,cases=0)),
    infection=list("rate=Beta*X*Y/N;",c(N=0,X=-1,Y=1,Z=0,cases=1)),
    recovery=list("rate=gamma*Y;",c(N=0,X=0,Y=-1,Z=1,cases=0)),
    hmax=0.01
  ),
  rmeasure=Csnippet("reports=rbinom(cases,rho);"),
  paramnames=c("rho","mu","Beta","gamma"),
  statenames=c("N","X","Y","Z","cases"),
  accumvars=c("cases"),
  params=c(X.0=392,Y.0=8,Z.0=0,N.0=400,cases.0=0,
           mu=0.02,Beta=60,gamma=365/13,rho=0.5)
) -> sir

simulate(sir,nsim=10,format="data.frame") -> sims

library(ggplot2)
library(reshape2)
ggplot(melt(sims,id=c("time",".id")),
       aes(x=time,y=value,group=.id,color=.id))+
  geom_line()+
  guides(color=FALSE)+
  facet_grid(variable~.,scales="free_y")