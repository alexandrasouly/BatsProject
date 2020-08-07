library(pomp)
# stochastic discrete SIR model in POMP
# no pop dyn.

# stochastic stepping process
stochStep <- Csnippet("
  p_SI = 1 - exp(-Beta*I/N);
  p_IR = 1 - exp(-gamma);

  n_SI = rbinom(S,p_SI);
  n_IR = rbinom(I,p_IR);

  S = S - n_SI;
  I = I + n_SI - n_IR;
  R = R + n_IR;
  H = n_SI;
")

# deterministic skeleton
skel <- Csnippet("
  DS = S - Beta*S*I/N;
  DI = I + Beta*S*I/N-gamma*I;
  DR = R + gamma*I;
  DH = H + Beta*S*I/N;
")

# initial values of states
init1 <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  H = 0;
  ")

params1 <- c(Beta=1, gamma = 1/13, N= 763)


closed_sir <- pomp(data=data.frame(time=1:50,data=NA),
                   times="time",
                   t0=0,
                   rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                   skeleton = map(skel),
                   rinit=init1,
                   statenames=c("S","I","R", "p_SI", "p_IR", "n_SI", "n_IR", "H"),
                   paramnames=c("Beta","gamma","N"),
                   accumvars = "H")




sim <- simulate(closed_sir,params=params1,format = "data.frame", nsim = 100)

traj <- trajectory(closed_sir,params=params1,format="data.frame")

created_plot <-ggplot(NULL)+
geom_line(data=sim, alpha = 0.05,aes(x=time, y = I, group= factor(.id), colour = "Stochastic infected"))+
geom_line(data=sim, alpha = 0.05,aes(x=time, y = R, group= factor(.id), colour = "Stochastic recovered"))+
geom_line(data=sim, alpha = 0.05,aes(x=time, y = S, group= factor(.id), colour = "Stochastic susceptible"))+
geom_line(data=sim, alpha = 0.05,aes(x=time, y = H, group= factor(.id), colour = "Stochastic new infections"))+
  
scale_color_discrete(name="Legend")+
geom_line(data=traj, aes(x=time, y=I))+
geom_line(data=traj, aes(x=time, y=R))+
geom_line(data=traj, aes(x=time, y=S))+
geom_line(data=traj, aes(x=time, y=H))+
labs(title="Stochastic SIR model")

