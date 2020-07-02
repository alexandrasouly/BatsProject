library(pomp)
# stochastic discrete SIR model in POMP
# with birth and deaths included

# stochastic stepping process
stochStep <- Csnippet("
  p_SI = 1 - exp(-Beta*I/N);
  p_IR = 1 - exp(-gamma);
  p_b  = 1 - exp(-b);
  p_d  = 1 - exp(-d);

  n_SI = rbinom(S,p_SI);
  n_IR = rbinom(I,p_IR);
  
  b_N  = rbinom(N, p_b);
  d_I  = rbinom(I, p_d);
  d_R  = rbinom(R, p_d);
  d_S  = rbinom(S, p_d);
  
  S = S - n_SI        + b_N  - d_S;
  I = I + n_SI - n_IR        - d_I;
  R = R        + n_IR        - d_R;

")

# deterministic skeleton
skel <- Csnippet("
  DS = S - Beta*S*I/N           + b*N - d*S ;
  DI = I + Beta*S*I/N - gamma*I       - d*I;
  DR = R              + gamma*I       - d*R;
")

# initial values of states
init1 <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  ")

params1 <- c(Beta=1, gamma = 1/13, N= 763, b=0.05, d = 0.05)


closed_sir <- pomp(data=data.frame(time=1:50,data=NA),
                   times="time",
                   t0=0,
                   rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                   skeleton = map(skel),
                   rinit=init1,
                   statenames=c("S","I","R", "p_SI", "p_IR", "p_b", "p_d", "n_SI", "n_IR", "b_N", "d_I", "d_S", "d_R"),
                   paramnames=c("Beta","gamma","N", "b", "d"))


sim <- simulate(closed_sir,params=params1,format = "data.frame", nsim = 100)

traj <- trajectory(closed_sir,params=params1,format="data.frame")

created_plot <-ggplot(NULL)+
  geom_line(data=sim, alpha = 0.05,aes(x=time, y = I, group= factor(.id), colour = "Stochastic infected"))+
  geom_line(data=sim, alpha = 0.05,aes(x=time, y = R, group= factor(.id), colour = "Stochastic recovered"))+
  geom_line(data=sim, alpha = 0.05,aes(x=time, y = S, group= factor(.id), colour = "Stochastic susceptible"))+
  scale_color_discrete(name="Legend")+
  geom_line(data=traj, aes(x=time, y=I))+
  geom_line(data=traj, aes(x=time, y=R))+
  geom_line(data=traj, aes(x=time, y=S))+
  labs(title="Stochastic SIR model")

