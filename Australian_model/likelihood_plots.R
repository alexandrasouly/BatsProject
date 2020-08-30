pdf("model4_if_plots.pdf")
#plotting the traces
cat("Model 4 Iterated filtering on all Clunes data using 50 iter, 2000 particles")
mf3iter2%>%
  traces() %>%
  melt() %>%
  filter(variable %in% c("loglik", "R0", "gamma_val", "zeta", "c",
                         "s_v", "phi_v", "disp", "d") )%>%
  ggplot(aes(x=iteration,y=value,group=L1,color=L1))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  guides(color=FALSE) -> pl1

plot(pl1)

# stating values 
pairs(~R0+rho_val+epsilon_val+zeta+c+s+phi+s_v+phi_v+disp+d, data = guesses, pch = 20, main = "Model 4 fitting starting values", sub= " 5x Aaron's posterior values and 25 point from a Sobol grid of 95% CI")


# starting and end values
all <- ldply(list(guess=guesses,result=subset(lik_list,loglik>max(loglik)-20)))
pairs(~R0+rho_val+epsilon_val+zeta+
        c+s+phi+s_v+phi_v+disp+d,data=all,col=ifelse(all$.id=="guess",grey(0.5),"red"),pch=20,main ="Model 4 high likelihood areas" )

# end values loglik
pairs(~R0+rho_val+epsilon_val+zeta+
        c+s_v+phi_v+disp+d,data=model4_likelihoods_s_fixed,pch=20, main = "Model 4 loglikelihoods")
  

#######################################################

gather(model6_likelihoods_s_fixed, key, value, -loglik) %>%
  filter(loglik> -59) %>%
  filter(key %in% c( "R0", "gamma_val",  "zeta", "c",
                              "s_v", "phi_v", "disp", "d")) %>%
ggplot(aes(value, loglik), main = "Results of model 4 IF") + 
   geom_point() +
  facet_wrap(~ key , scales="free_x", ncol=4) -> pl4

plot(pl4)

dev.off()
  