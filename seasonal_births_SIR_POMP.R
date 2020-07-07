#stoch bd model with seasonal births
library(pomp)
stochStep <- Csnippet("
  double b = b0*(1+b1*cos(2*M_PI*t/365));
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
  double b = b0*(1+b1*cos(2*M_PI*t/365));
  DS = S - Beta*S*I/N           + b*N - d*S ;
  DI = I + Beta*S*I/N - gamma*I       - d*I;
  DR = R              + gamma*I       - d*R;
")

# initial values of states
init1 <- Csnippet("
  S = N-5;
  I = 5;
  R = 0;
  ")

params1 <- c(Beta=1, b0=0.1, b1=0.4, gamma = 1/13, N= 763, d = 0.1)


closed_sir <- pomp(data=data.frame(time=1:700,data=NA),
                   times="time",
                   t0=0,
                   rprocess=discrete_time(step.fun=stochStep,delta.t=1),
                   skeleton = map(skel),
                   rinit=init1,
                   statenames=c("S","I","R", "p_SI", "p_IR", "p_b", "p_d", "n_SI", "n_IR", "b_N", "d_I", "d_S", "d_R"),
                   paramnames=c("Beta","b0","b1","gamma","N", "d"))


sim <- simulate(closed_sir,params=params1,format = "data.frame", nsim = 100)

traj <- trajectory(closed_sir,params=params1,format="data.frame")

birth_plot <- ggplot(NULL)+ geom_line(data=traj,aes(x=time, y = I+S+R,))+ylim(0,NA)
birth_plot

cols <- c( "Stochastic infected"= "red", "Stochastic susceptible" = "blue", "Stochastic recovered" = "darkgreen")

created_plot <-ggplot(NULL)+
  geom_line(data=sim, alpha = 0.05,aes(x=time, y = I, group= factor(.id), colour = "Stochastic infected"))+
  geom_line(data=sim, alpha = 0.05,aes(x=time, y = R, group= factor(.id), colour = "Stochastic recovered"))+
  geom_line(data=sim, alpha = 0.05,aes(x=time, y = S, group= factor(.id), colour = "Stochastic susceptible"))+
  scale_colour_manual(values = cols)+
  geom_line(data=traj, aes(x=time, y=I))+
  geom_line(data=traj, aes(x=time, y=R))+
  geom_line(data=traj, aes(x=time, y=S))+
  labs(title="Stochastic SIR model")

