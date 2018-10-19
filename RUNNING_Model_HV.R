

rm(list=ls())
###################################################################################
library(jagsUI)

#####################
### SIMULATING SOME DATA - TO BE USED AS INITIAL VALUES HERE; 

setwd("C:/Users/erlend.nilsen/OneDrive - NINA/Prosjekter/Arbeidsbord/Villrein/Integrated Population Models/PrevalenceProject")

source("Appendix 3 Simulated_Data.R")
Simulated_data <- SimPop(PHI3=0.9, mean_F=0.5, p1=0.5, p2=0.5, bias1=1, T=40)                                 # Run script to simmulate data 


#####################
### READING AND PREPARING DATA FROM SNØHETTA AND
### KNUTSHØ

d <- read.csv("Data_Hardangervidda_2016.csv", sep="\t", header = T)

#d[30, c(11, 12, 13)] <- NA




#########################
# Bundle data
# Hardangervidda data; 
bugs.data.HV <- list(SU = d[1:33,"SU"],
                        J = d[1:33,"J"],
                        H0f = d[1:33,"H0f"],
                        H0m = d[1:33,"H0m"],
                        H1f = d[1:33,"H1f"],
                        H1m = d[1:33,"H1m"],
                        Hadf = d[1:33,"Hadf"],
                        Hadm = d[1:33,"Hadm"],
                        T = 32,
                        y=d[1:33, "N"],
                        C0=d[1:33,"C0"],
                        Cm=d[1:33,"Cm"],
                        Cf=d[1:33,"Cf"])

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
                            p2 = runif(32, 0.6, 0.9),
                            p1 = runif(32, 0.6, 0.9),
                            N0f = Simulated_data$N[1, 4:35]*250, 
                            N0m = Simulated_data$N[2, 4:35]*250, 
                            N1f = Simulated_data$N[3, 4:35]*250,
                            N1m = Simulated_data$N[4, 4:35]*250,
                            Nadf = Simulated_data$N[5, 4:35]*250,
                            Nadm = Simulated_data$N[6, 4:35]*250) }  


HV_M1_var <- jags(data=bugs.data.HV, inits=inits, 
                 parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                      "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                 model.file="jagsMod_M1_HV_variant.bug",n.chain=nchains, n.iter=niter, 
                 n.burnin=nburn, parallel=TRUE)

HV_M2_var <- jags(data=bugs.data.HV, inits=inits, 
                  parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                       "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                  model.file="jagsMod_M2_HV_variant.bug",n.chain=nchains, n.iter=niter, 
                  n.burnin=nburn, parallel=TRUE)



str(HV_M1_var$sims.list)
summary(HV_M1_var)

##########################################
mean <- upper <- lower <- numeric()

for(i in 1:32){
  upper[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.975)
  lower[i] <- quantile(HV_M1_var$sims.list$Ntot[,i], 0.025)
  mean[i] <- mean(HV_M1_var$sims.list$Ntot[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2016, mean, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,15000), ylab="Population size (N[t])", 
     xlab="Year", xlim=c(1985, 2016))
lines(1985:2016, mean, col="dark red",  type="l")
arrows(x0=1985:2016, x1=1985:2016, y0=lower, y1=upper, angle=0)
arrows(x0=1985:2016, x1=1985:2016, y1=lower, y0=upper, angle=0)

points(x=1985:2016, y=d$N[1:32], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)



##########################################
mean2 <- upper2 <- lower2 <- numeric()

for(i in 1:32){
  upper2[i] <- quantile(HV_M2_var$sims.list$Ntot[,i], 0.975)
  lower2[i] <- quantile(HV_M2_var$sims.list$Ntot[,i], 0.025)
  mean2[i] <- mean(HV_M2_var$sims.list$Ntot[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2016, mean2, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,15000), ylab="Population size (N[t])", 
     xlab="Year", xlim=c(1985, 2016))
lines(1:32, mean, type="p", col="dark red", pch=16)
arrows(x0=1985:2016, x1=1985:2016, y0=lower2, y1=upper2, angle=0)
arrows(x0=1985:2016, x1=1985:2016, y1=lower2, y0=upper2, angle=0)

points(x=1985:2016, y=d$N[1:32], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)



##########################################
mean3 <- upper3 <- lower3 <- numeric()

for(i in 1:31){
  upper3[i] <- quantile(HV_M2_var$sims.list$f[,i], 0.975)
  lower3[i] <- quantile(HV_M2_var$sims.list$f[,i], 0.025)
  mean3[i] <- mean(HV_M2_var$sims.list$f[,i])
}


par(bty="l", mar=c(4,4,3,1), cex=0.8)
plot(1985:2015, mean3, type="l", col="dark red", lwd=1, lty=1, ylim=c(0,1), ylab="f (calves/female)", 
     xlab="Year", xlim=c(1985, 2015))
lines(1:31, mean, type="p", col="dark red", pch=16)
arrows(x0=1985:2015, x1=1985:2015, y0=lower3, y1=upper3, angle=0)
arrows(x0=1985:2015, x1=1985:2015, y1=lower3, y0=upper3, angle=0)

points(x=1985:2016, y=d$N[1:32], pch=16, cex=2, col="dark orange")
mtext("HARDANGERVIDDA", outer=F, adj=0.01, line=1)


###########################################################################################################
###########################################################################################################

###########################################################################################################
### PREPARING JAGS MODEL RUNS; 

rm(list=ls())
###################################################################################
library(jagsUI)


setwd("C:/Users/erlend.nilsen/OneDrive - NINA/Prosjekter/Arbeidsbord/Villrein/Integrated Population Models/PrevalenceProject")

source("M:/My Documents/Prosjekter/Arbeidsbord/Villrein/Integrated Population Models/Models in PlosOne/PlosOneAppendixes/Appendix 3 Simulated_Data.R")
source("Appendix 3 Simulated_Data.R")

#####################
### READING AND PREPARING DATA FROM SNØHETTA AND
### KNUTSHØ

d <- read.csv("Data_Hardangervidda_2016.csv", sep="\t", header = T)



#########################
# Bundle data
# Hardangervidda data; 




start <- 16
stop <- 29

Time_s <- (stop-start)+1

bugs.data.HV <- list(SU = d[start:stop,"SU"],
                     J = d[start:stop,"J"],
                     H0f = d[start:stop,"H0f"],
                     H0m = d[start:stop,"H0m"],
                     H1f = d[start:stop,"H1f"],
                     H1m = d[start:stop,"H1m"],
                     Hadf = d[start:stop,"Hadf"],
                     Hadm = d[start:stop,"Hadm"],
                     T = Time_s,
                     y=d[start:stop, "N"],
                     C0=d[start:stop,"C0"],
                     Cm=d[start:stop,"Cm"],
                     Cf=d[start:stop,"Cf"])



# MCMC settings
niter <- 150000
nthin <- 3
nburn <- 50000
nchains <- 3


# Initial values
inits <- function() { list( phi3 = runif(1, 0.9, 0.99), 
                            mean.f = runif(1, 2, 3), 
                            mean.phi1 = runif(1, 2, 3), 
                            sd_phi1 = runif(1, 0, 1),
                            sd_f = runif(1, 0, 1),
                            p2 = runif(Time_s, 0.6, 0.9),
                            p1 = runif(Time_s, 0.6, 0.9),
                            N0f = Simulated_data$N[1, start:stop]*250, 
                            N0m = Simulated_data$N[2, start:stop]*250, 
                            N1f = Simulated_data$N[3, start:stop]*250,
                            N1m = Simulated_data$N[4, start:stop]*250,
                            Nadf = Simulated_data$N[5, start:stop]*250,
                            Nadm = Simulated_data$N[6, start:stop]*250) }  


HV_M1_var <- jags(data=bugs.data.HV, inits=inits, 
                  parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                       "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                  model.file="jagsMod_M1_HV_variant.bug",n.chain=nchains, n.iter=niter, 
                  n.burnin=nburn, parallel=TRUE)

HV_M2_var <- jags(data=bugs.data.HV, inits=inits, 
                  parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                       "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                  model.file="jagsMod_M2_HV_variant.bug",n.chain=nchains, n.iter=niter, 
                  n.burnin=nburn, parallel=TRUE)


HV_M1_temp <- jags(data=bugs.data.HV, inits=inits, 
                  parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                       "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                  model.file="jagsMod_M1_use_Preval_HV.bug",n.chain=nchains, n.iter=niter, 
                  n.burnin=nburn, parallel=TRUE)








