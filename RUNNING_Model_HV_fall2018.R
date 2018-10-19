

rm(list=ls())
###################################################################################
library(jagsUI)

#####################
### SIMULATING SOME DATA - TO BE USED AS INITIAL VALUES HERE; 


source("Appendix 3 Simulated_Data.R")
Simulated_data <- SimPop(PHI3=0.9, mean_F=0.5, p1=0.5, p2=0.5, bias1=1, T=40)                                 # Run script to simmulate data 


#####################
### READING AND PREPARING DATA FROM HARDANGERVIDDA

d <- read.csv("data/Data_Hardangervidda_2018.csv", sep=";", header = T)

#d[30, c(11, 12, 13)] <- NA




#########################
# Bundle data
# Hardangervidda data;
T_max <- 34

bugs.data.HV <- list(SU = d[1:T_max,"SU"],
                        J = d[1:T_max,"J"],
                        H0f = d[1:T_max,"H0f"],
                        H0m = d[1:T_max,"H0m"],
                        H1f = d[1:T_max,"H1f"],
                        H1m = d[1:T_max,"H1m"],
                        Hadf = d[1:T_max,"Hadf"],
                        Hadm = d[1:T_max,"Hadm"],
                        T = T_max,
                        y_s=d[1:T_max, "N_s"],
                        y_w=d[2:T_max, "N_v"],
                        C0=d[1:T_max,"C0"],
                        Cm=d[1:T_max,"Cm"],
                        Cf=d[1:T_max,"Cf"])

#####
bugs.data.HV2 <- list(SU = d[1:T_max,"SU"],
                     J = d[1:T_max,"J"],
                     H0f = d[1:T_max,"H0f"],
                     H0m = d[1:T_max,"H0m"],
                     H1f = d[1:T_max,"H1f"],
                     H1m = d[1:T_max,"H1m"],
                     Hadf = d[1:T_max,"Hadf"],
                     Hadm = d[1:T_max,"Hadm"],
                     T = T_max,
                     y_w=d[2:T_max, "N_v"],
                     C0=d[1:T_max,"C0"],
                     Cm=d[1:T_max,"Cm"],
                     Cf=d[1:T_max,"Cf"])


#####
bugs.data.HV3 <- list(SU = d[1:T_max,"SU"],
                     J = d[1:T_max,"J"],
                     H0f = d[1:T_max,"H0f"],
                     H0m = d[1:T_max,"H0m"],
                     H1f = d[1:T_max,"H1f"],
                     H1m = d[1:T_max,"H1m"],
                     Hadf = d[1:T_max,"Hadf"],
                     Hadm = d[1:T_max,"Hadm"],
                     T = T_max,
                     y_s=d[1:T_max, "N_s"],
                     C0=d[1:T_max,"C0"],
                     Cm=d[1:T_max,"Cm"],
                     Cf=d[1:T_max,"Cf"])


###########################################################################################################
### PREPARING JAGS MODEL RUNS; 

# MCMC settings
niter <- 250000
nthin <- 3
nburn <- 50000
nchains <- 3


# Initial values
inits <- function() { list( phi3 = runif(1, 0.9, 0.99), 
                            mean.f = runif(1, 2, 3), 
                            mean.phi1 = runif(1, 2, 3), 
                            sd_phi1 = runif(1, 0, 1),
                            sd_f = runif(1, 0, 1),
                            p2 = runif(34, 0.6, 0.9),
                            p1 = runif(34, 0.6, 0.9),
                            N0f = Simulated_data$N[1, 4:37]*250, 
                            N0m = Simulated_data$N[2, 4:37]*250, 
                            N1f = Simulated_data$N[3, 4:37]*250,
                            N1m = Simulated_data$N[4, 4:37]*250,
                            Nadf = Simulated_data$N[5, 4:37]*250,
                            Nadm = Simulated_data$N[6, 4:37]*250) }  


HV_M1_var <- jags(data=bugs.data.HV, inits=inits, 
                 parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                      "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                 model.file="jagsMod_M1_HV_variant_fall2018.bug",n.chain=nchains, n.iter=niter, 
                 n.burnin=nburn, parallel=TRUE)

HV_M1_var2 <- jags(data=bugs.data.HV, inits=inits, 
                  parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                       "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                  model.file="jagsMod_M1_HV_variant2_fall2018.bug",n.chain=nchains, n.iter=niter, 
                  n.burnin=nburn, parallel=TRUE)

HV_M1_var3 <- jags(data=bugs.data.HV, inits=inits, 
                   parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                        "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                   model.file="jagsMod_M1_HV_variant3_fall2018.bug",n.chain=nchains, n.iter=niter, 
                   n.burnin=nburn, parallel=TRUE)






##########################################
### Summer abundance (Nt)

mean <- upper <- lower <- numeric()

for(i in 1:34){
  upper[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.975)
  lower[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.025)
  mean[i] <- mean(HV_M1_var$sims.list$Ntot[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2018, mean, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,15000), ylab="Population size (N[t])", 
     xlab="Year", xlim=c(1985, 2018))
lines(1985:2018, mean, col="dark red",  type="l")
arrows(x0=1985:2018, x1=1985:2018, y0=lower, y1=upper, angle=0)
arrows(x0=1985:2018, x1=1985:2018, y1=lower, y0=upper, angle=0)

points(x=1985:2018, y=d$N_s[1:34], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)

##########################################
### Winter abundance (Xt)

mean <- upper <- lower <- numeric()

for(i in 1:33){
  upper[i] <- quantile(HV_M1_var$sims.list$Xtot[,i], 0.975)
  lower[i] <- quantile(HV_M1_var$sims.list$Xtot[,i], 0.025)
  mean[i] <- mean(HV_M1_var$sims.list$Xtot[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1986:2018, mean, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,15000), ylab="Population size (X[t])", 
     xlab="Year", xlim=c(1986, 2018))
lines(1986:2018, mean, col="dark red",  type="l")
arrows(x0=1986:2018, x1=1986:2018, y0=lower, y1=upper, angle=0)
arrows(x0=1986:2018, x1=1986:2018, y1=lower, y0=upper, angle=0)

points(x=1986:2018, y=d$N_v[2:34], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)


##########################################
mean3 <- upper3 <- lower3 <- numeric()

for(i in 1:33){
  upper3[i] <- quantile(HV_M1_var$sims.list$f[,i], 0.975)
  lower3[i] <- quantile(HV_M1_var$sims.list$f[,i], 0.025)
  mean3[i] <- mean(HV_M1_var$sims.list$f[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2017, mean3, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,1), ylab="f (calves/female)", 
     xlab="Year", xlim=c(1985, 2017))
lines(1:33, mean, type="p", col="dark red", pch=16)
arrows(x0=1985:2017, x1=1985:2017, y0=lower3, y1=upper3, angle=0)
arrows(x0=1985:2017, x1=1985:2017, y1=lower3, y0=upper3, angle=0)

points(x=1985:2017, y=d$N[1:32], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)

##########################################
### Combined (Xt)

## Summer pop
mean1 <- upper1 <- lower1 <- numeric()

for(i in 1:34){
  upper1[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.975)
  lower1[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.025)
  mean1[i] <- mean(HV_M1_var$sims.list$Ntot[,i])
}


## Winter pop
mean <- upper <- lower <- numeric()

for(i in 1:33){
  upper[i] <- quantile(HV_M1_var$sims.list$Xtot[,i], 0.975)
  lower[i] <- quantile(HV_M1_var$sims.list$Xtot[,i], 0.025)
  mean[i] <- mean(HV_M1_var$sims.list$Xtot[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2018, mean1, type="n", col="dark red", lwd=1, lty=1, ylim=c(0,15000), ylab="Population size", 
     xlab="Year", xlim=c(1985, 2018))
lines(1986:2018, mean, col=1,  type="l", lwd=2)
lines(seq(1985.3, 2018.3, by=1), mean1, col="dark orange",  type="l", lwd=2, lty=3)

arrows(x0=1986:2018, x1=1986:2018, y0=lower, y1=upper, angle=0, col=1)
arrows(x0=1986:2018, x1=1986:2018, y1=lower, y0=upper, angle=0, col=1)

arrows(x0=seq(1985.3, 2018.3, by=1), x1=seq(1985.3, 2018.3, by=1), y0=lower1, y1=upper1, angle=0, col="dark orange")
arrows(x0=seq(1985.3, 2018.3, by=1), x1=seq(1985.3, 2018.3, by=1), y1=lower1, y0=upper1, angle=0, col="dark orange")

points(x=1986:2018, y=d$N_v[2:34], pch=16, cex=2, col=1)
points(x=seq(1985.3, 2018.3, by=1), y=d$N_s[1:34], pch=16, cex=2, col="dark orange")

mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)


###########################################
##########################################
### Comparing predctions based on models with 
### summer and winter counts, and with winter 
### counts only and summer counts only

mean <- upper <- lower <- numeric()

for(i in 1:34){
  upper[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.975)
  lower[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.025)
  mean[i] <- mean(HV_M1_var$sims.list$Ntot[,i])
}

mean2 <- upper2 <- lower2 <- numeric()

for(i in 1:34){
  upper2[i] <- quantile(HV_M1_var2$sims.list$Ntot[,i], 0.975)
  lower2[i] <- quantile(HV_M1_var2$sims.list$Ntot[,i], 0.025)
  mean2[i] <- mean(HV_M1_var2$sims.list$Ntot[,i])
}

mean3 <- upper3 <- lower3 <- numeric()

for(i in 1:34){
  upper3[i] <- quantile(HV_M1_var3$sims.list$Ntot[,i], 0.975)
  lower3[i] <- quantile(HV_M1_var3$sims.list$Ntot[,i], 0.025)
  mean3[i] <- mean(HV_M1_var3$sims.list$Ntot[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2018, mean, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,15000), ylab="Population size (N[t])", 
     xlab="Year", xlim=c(1985, 2018))
lines(1985:2018, mean, col="dark red",  type="l")
lines(1985:2018, mean2, col=1,  type="l", lty=3)
lines(1985:2018, mean3, col="green",  type="l", lty=2)


arrows(x0=1985:2018, x1=1985:2018, y0=lower, y1=upper, angle=0)
arrows(x0=1985:2018, x1=1985:2018, y1=lower, y0=upper, angle=0)

points(x=1985:2018, y=d$N_s[1:34], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)





