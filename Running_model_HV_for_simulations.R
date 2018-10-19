

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
### HARDANGERVIDDA

d <- read.csv("Data_Hardangervidda_2016.csv", sep="\t", header = T)



#########################
# Bundle data
# Hardangervidda data; 

Simulated_data <- SimPop(PHI3=0.9, mean_F=0.5, p1=0.5, p2=0.5, bias1=1, T=40)                                 # Run script to simmulate data 

############################################################################################################
###########################################################################################################
### PREPARING JAGS MODEL RUNS; 
###
## USING DATA FROM YEARS 2001-2017

start <- 17
stop <- 33

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


HV_M1_temp_structure <- jags(data=bugs.data.HV, inits=inits, 
                                          parameters.to.save=c("mean.phi1", "phi3", "phi1", "f", "p1", "p2", "Ntot", "Xtot", "Sjuv", "fert", "sd_f", "sd_phi1", "N0f", 
                                                               "N0m", "N1f", "N1m", "Nadf", "Nadm"), 
                                          model.file="jagsMod_M1_HV_variant.bug",n.chain=nchains, n.iter=niter, 
                                          n.burnin=nburn, parallel=TRUE)

summary(HV_M1_temp_structure)


#############################################################################



## Population structure 2017
res <- matrix(ncol=2, nrow=6)

res[1,1] <- HV_M1_temp_structure$mean$N0f[17]
res[2,1] <- HV_M1_temp_structure$mean$N0m[17]
res[3,1] <- HV_M1_temp_structure$mean$N1f[17]
res[4,1] <- HV_M1_temp_structure$mean$N1m[17]
res[5,1] <- HV_M1_temp_structure$mean$Nadf[17]
res[6,1] <- HV_M1_temp_structure$mean$Nadm[17]


res[1,2] <- HV_M1_temp_structure$sd$N0f[17]
res[2,2] <- HV_M1_temp_structure$sd$N0m[17]
res[3,2] <- HV_M1_temp_structure$sd$N1f[17]
res[4,2] <- HV_M1_temp_structure$sd$N1m[17]
res[5,2] <- HV_M1_temp_structure$sd$Nadf[17]
res[6,2] <- HV_M1_temp_structure$sd$Nadm[17]



## Population structure 2016
res <- matrix(ncol=2, nrow=6)

res[1,1] <- HV_M1_temp_structure$mean$N0f[16]
res[2,1] <- HV_M1_temp_structure$mean$N0m[16]
res[3,1] <- HV_M1_temp_structure$mean$N1f[16]
res[4,1] <- HV_M1_temp_structure$mean$N1m[16]
res[5,1] <- HV_M1_temp_structure$mean$Nadf[16]
res[6,1] <- HV_M1_temp_structure$mean$Nadm[16]


res[1,2] <- HV_M1_temp_structure$sd$N0f[16]
res[2,2] <- HV_M1_temp_structure$sd$N0m[16]
res[3,2] <- HV_M1_temp_structure$sd$N1f[16]
res[4,2] <- HV_M1_temp_structure$sd$N1m[16]
res[5,2] <- HV_M1_temp_structure$sd$Nadf[16]
res[6,2] <- HV_M1_temp_structure$sd$Nadm[16]










