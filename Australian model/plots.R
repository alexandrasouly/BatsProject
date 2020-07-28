plots <- function(sim_plus_data, x, sim){

  
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Im)) +
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Im, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none")->plIm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sm))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sm, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plSm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rm))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rm, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plRm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Em))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Em, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plEm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ij))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ij, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plIj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sj))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sj, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plSj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rj))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rj, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plRj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ej))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ej, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plEj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=In))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = In, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plIn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sn))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sn, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plSn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rn))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rn, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plRn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=En))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = En, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plEn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ma))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ma, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") ->plMa
  
  plots <- list(plIm, plRm, plEm, plSm, 
                plIj, plRj, plEj, plSj, 
                plIn, plRn, plEn, plSn,
                plMa)
  return(plots)
}
  
  