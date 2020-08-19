plots <- function(sim_plus_data, x, sim){
  
# preparing plotting 
  x$Time <- factor(x$time)
  sim_plus_data$Time <- factor(sim_plus_data$time)
  read.csv("seasons.csv", sep = " ") -> seasons
  seasons[2:9,] -> seasons_CLU
  seasons_CLU[1,2] <-198
  seasons_CLU[8,3] <- 1061
  months_breaks <- c(c("0","151","243"),
                     as.character(c(0,151,243)+365),
                    as.character(c(0,151,243)+2*365))
  month_labels <- c("Jan 2017",  "Jun", "Sep", 
                        "Jan 2018",  "Jun", "Sep",
                        "Jan 2019",  "Jun", "Sep")
  
  
  read.table("data/catching_data.txt") -> catching_data
  catching_data %>% 
    filter(site_abbrev == "CLU")%>%
    select(min_date, hen_prevalence) -> clu_catching
  names(clu_catching)[names(clu_catching) == "min_date"] <- "time"
# plotting  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Im)) +
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Im, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none")+
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 3000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plIm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sm))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sm, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none")+
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 3000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plSm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rm))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rm, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 3000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plRm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Em))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Em, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none")+
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plEm
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ij))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ij, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plIj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sj))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sj, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plSj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rj))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rj, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plRj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ej))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ej, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plEj
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=In))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = In, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none")+
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plIn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Sn))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Sn, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plSn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Rn))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Rn, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plRn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=En))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = En, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plEn
  
  ggplot(data=x,mapping=aes(x=time))+
    geom_line(aes(x=time, y=Ma))+
    geom_line(data=sim, alpha = 0.1,aes(x=time, y = Ma, group= factor(.id),
                                        colour = "Stochastic infected"))+
    theme(legend.position = "none") +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 2000, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))->plMa
  
  
  ggplot(data=x, aes(x= Time), main = "New infection cases")+
    geom_line(aes( y=H, group = 1)) +
    geom_line(data=sim, alpha = 0.01,
              aes(x=time, y = H,group= factor(.id),
              colour = "Stochastic new infected"))+
    theme(legend.position = "none") +
    geom_vline(xintercept = c(0,365,2*365,3*365)) +
    geom_rect(data = seasons, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = max(sim$H), alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))+
    scale_x_discrete(breaks=months_breaks,
                     labels = month_labels)+
    ylab("Daily infections")+
    theme(axis.text.x=element_text(angle=90,hjust=1)) ->plH
  
  ggplot(data = sim_plus_data, aes(x= time))+
    geom_line( alpha = 0.1,aes( y = sim_model_prev, group= factor(.id), colour = "Simulated model prevalence"
    ))+
    geom_line( alpha = 0.1,aes( y = true_test_prev, group= factor(.id), colour = "Clunes underroost prevalence"
    ))+
    geom_line(alpha = 0.1,aes( y = sim_test_prev, group= factor(.id), colour = "Simulated underroost prevalence"
    ))+
    geom_vline(xintercept = c(365,2*365)) +
    geom_rect(data = seasons_CLU, 
              aes(x = NULL, xmin = start_day, xmax = end_day, fill = season),
              ymin = 0, ymax = 1, alpha = 0.1)+
    scale_fill_manual(values = c("green", "blue"),
                      labels = c("Not Winter", "Winter"))+
    ylab("Prevalence proportion")+
    theme(axis.text.x=element_text(angle=90,hjust=1))-> plprev
    plprev + geom_point(data = clu_catching,aes(x = time, 
                    y = hen_prevalence))->plprev
  
  
  plots <- list(plIm, plRm, plEm, plSm, 
                plIj, plRj, plEj, plSj, 
                plIn, plRn, plEn, plSn,
                plMa, plH, plprev)
  return(plots)
  
  
  
  
  
  
}
  
  