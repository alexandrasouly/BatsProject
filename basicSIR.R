library(deSolve)
library(ggplot2)

SIRmodel <- function(t,y, parms){
  with(as.list(c(y, parms)), {
    dS <- (-beta) * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I

    return(list(c(dS, dI, dR)))
    }
  )
}

init <- c(S = 1-1e-6, I = 1e-6, R=0.0)
params <- c(beta = 1.4247, gamma = 0.14286)
times <- seq(0, 70, by = 1)



results <- as.data.frame(ode(y = init, times = times, func = SIRmodel, parms = params))


res_plot<-ggplot(results,aes(x=time))+
  geom_line(aes(y=S,colour="Susceptible"))+
  geom_line(aes(y=I,colour="Infected"))+
  geom_line(aes(y=R,colour="Recovered"))+
  scale_colour_manual("Compartments",
                      breaks=c("Susceptible","Infected","Recovered"),
                      values=c("red","green","blue"))

res_plot
