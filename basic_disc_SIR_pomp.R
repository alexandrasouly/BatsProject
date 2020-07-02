library(pomp)

# deterministic discrete SIR model in POMP
# no pop dyn.

closed_sir_ode <- Csnippet("
  DS = S - Beta*S*I/N;
  DI = I + Beta*S*I/N-gamma*I;
  DR = R + gamma*I;
")

init1 <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  ")

closed_sir <- pomp(data=data.frame(time=1:50,data=NA),
                   times="time",t0=0,
                   skeleton=map(closed_sir_ode),
                   rinit=init1,
                   statenames=c("S","I","R"),
                   paramnames=c("Beta","gamma","N"))

params1 <- c(Beta=1,gamma=(1/13),N=763)

x <- trajectory(closed_sir,params=params1,format="data.frame")

library(ggplot2)
ggplot(data=x,mapping=aes(x=time,y=I))+geom_line()


################################################################################

# deterministic discrete SIR model in POMP
# births and deaths

open_sir_ode <- Csnippet("
  DS = S-Beta*S*I/(N)+mu*(N-S);
  DI = I+Beta*S*I/(N)-gamma*I-mu*I;
  DR = R+gamma*I-mu*R;
")

init2 <- Csnippet("
  S = S_0;
  I = I_0;
  R = N-S_0-I_0;
")

open_sir <- pomp(data=data.frame(time=seq(0,50,by=1),cases=NA),
                 times="time",t0=-1,
                 skeleton=map(open_sir_ode),
                 rinit=init2,
                 statenames=c("S","I","R"),
                 paramnames=c("Beta","gamma","mu","S_0","I_0","N")
            )

params3 <- c(mu=1/100,Beta=1,gamma=1/13,
             N=763,S_0=762,I_0=1)
x <- trajectory(open_sir,params=params3,format="d")

library(ggplot2)
ggplot(data=x,mapping=aes(x=time,y=I))+geom_line()
ggplot(data=x,mapping=aes(x=S,y=I))+geom_path()
