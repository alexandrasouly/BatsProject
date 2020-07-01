# SIR deterministic models in deSolve package
# one model with no pop dyn, and one with birth and deaths added 

library(deSolve)
library(ggplot2)


SIRmodel <- function(t,y, parms){
  # no pop dyn
  # returns a list of grads for ode() input
  with(as.list(c(y, parms)), {
    dS <- (-beta) * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I

    return(list(c(dS, dI, dR)))
    }
  )
}

SIR_birth_death <- function(t,y, parms){
  # pop dyn added
  # returns a list of grads for ode() input
  with(as.list(c(y, parms)), {
    dS <- b - beta * S * I             - d * S
    dI <-     beta * S * I - gamma * I - d * I
    dR <-                    gamma * I - d * R
    
    return(list(c(dS, dI, dR)))
  }
  )
}


solve_and_plot <- function(init, params, times, func){
  # solves ODES specified in func
  
  results <- as.data.frame(ode(y = init, times = times, func = func, parms = params))
  res_plot<-ggplot(results,aes(x=time))+
    geom_line(aes(y=S,colour="Susceptible"))+
    geom_line(aes(y=I,colour="Infected"))+
    geom_line(aes(y=R,colour="Recovered"))+
    ylab(label="Proportion")+
    xlab(label="Day")+
    ggtitle("Basic deterministic SIR Model")+
    scale_colour_manual("Compartments",
                        breaks=c("Susceptible","Infected","Recovered"),
                        values=c("red","green","blue"))
  
    return(list(results, res_plot))
}

# example of using it

init <- c(S = 1- 1/763, I = 1/763, R=0.0)
params <- c(beta = 1, gamma = 1/13)
times <- seq(0, 50, by = 1)

out <- solve_and_plot(init, params, times, SIRmodel)

