plots <- function(sim_plus_data, x, sim){
  
  # plotting prevalence from tests and actual from the model
  
  ggplot(NULL)+
    geom_line(data=sim_plus_data, alpha = 0.1,aes(x=time, y = sim_model_prev, group= factor(.id), colour = "Simulated model prevalence"
    ))+
    geom_line(data=sim_plus_data, alpha = 0.1,aes(x=time, y = sim_test_prev, group= factor(.id), colour = "Simulated tests prevalence"
    ))->plprev
  
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Im)) +
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Im, group= factor(.id),
                                        colour = "Stochastic infected"))->plIm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sm))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sm, group= factor(.id),
                                        colour = "Stochastic infected")) ->plSm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rm))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rm, group= factor(.id),
                                        colour = "Stochastic infected")) ->plRm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ij))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ij, group= factor(.id),
                                        colour = "Stochastic infected")) ->plIj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sj))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sj, group= factor(.id),
                                        colour = "Stochastic infected")) ->plSj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rj))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rj, group= factor(.id),
                                        colour = "Stochastic infected")) ->plRj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=In))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = In, group= factor(.id),
                                        colour = "Stochastic infected")) ->plIn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sn))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sn, group= factor(.id),
                                        colour = "Stochastic infected")) ->plSn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rn))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rn, group= factor(.id),
                                        colour = "Stochastic infected")) ->plRn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ma))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ma, group= factor(.id),
                                        colour = "Stochastic infected")) ->plMa
  
  plots <- list(plprev, 
                plIm, plRm, plSm, 
                plIj, plRj, plSj, 
                plIn, plRn, plSn,
                plMa)
  return(plots)
}
  
  