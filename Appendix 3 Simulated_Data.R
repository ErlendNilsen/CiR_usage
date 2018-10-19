


#############################################################################################################################
#####
##### THIS SCRIPT USED TO SIMULATE DATA FOR ANALYSIS IN THE 
##### HIERARCHICAL INTEGRATED POP. MODEL PRESENTED IN NILSEN ET AL.; 


SimPop  <- function(m=1, phi1t="Nei", ft="Nei", p1t="Nei", p2t="Nei", T=31, PHI3 = 0.95, mean_F=0.85, f_sd=0.1, mean_phi1 = 0.85, phi1_sd=0.1, 
                     p1=0.6, p1_sd=0.1, p2=0.4, p2_sd=0.1, 
                     bias1 = 1, bias2 = 1,  bias3=1, error.dist="pois", p.obs=NULL){
  # m: Multiplier - for the initial poulation size
  # T: Number of time steps to simulate
  # PHI3: Adult annual survival probability
  # mean_F: Mean number of calves pr. female, 
  # PHI1: Juvenile summer survival probability
  # PHI1: JUvenile overwinter survival probability
  # p1: Observation probability for calf surveys in spring
  # p2: Observation probability for structur surveys post harvest  
  # bias1: Bias in observation data; Calf surveys in spring. bias1>1 implies higher detection probability 
  # for females/yearlings. bias1<1 implies higher detection probability 
  # for calves.
  # bias2: Bias in observation data; Structural surveys post harvest. bias2>1 implies higher detection probability for
  # females compared to calves and males. bias2<1 implies lower detection probability for calves than females and males
  # bias3: Bias in observation data; Structural surveys post harvest. bias3>1 implies higher detection probability for
  # males compared to calves and females
  # error.dist: Error distribution for count data
  # p.obs: Observation probaiblity for count data
  
  ## POPULATION VECTORS
  N <- matrix(ncol=T, nrow=6)               ## Pre harvest pop. vector. No monitoring
  X <- matrix(ncol=T, nrow=6)               ## Post harvest pop. vector. Monitored by str. surveys (calves, males, females) and potentially tot. surveys
  N_juv <- matrix(ncol=T, nrow=2)          ## Spring pop. vector. Monitored by recruitment surveys (calves vs yearlings and females)
  N_tot <- X_tot <- N_obs <- J <- SU <- C0 <- Cf <- Cm <- matrix(ncol=T, nrow=1)
  
  ## OBSERVATION PARAMETERS; PRE-HARVEST SURVEYS
  pa <- matrix(ncol=T, nrow=1) 
  
  if(p1t=="Ja"){
    for(i in 1:T){
    pa[i] <- betaval(p1, p1_sd)  
    }
  }
  
  if(p1t=="Nei"){
    for(i in 1:T){
      pa[i] <- p1  
    }
  }
  
  ## OBSERVATION PARAMETERS; POST-HARVEST SURVEYS
  pb <- matrix(ncol=T, nrow=1) 
  
  if(p2t=="Ja"){
    for(i in 1:T){
      pb[i] <- betaval(p1, p1_sd)  
    }
  }
  
  if(p2t=="Nei"){
    for(i in 1:T){
      pb[i] <- p2  
    }
  }
  
  
  ## DEMOGRAPHIC PARAMETERS
# fecundity;   
  f <- matrix(ncol=T, nrow=1) 
  
  if(ft=="Ja"){
    for(i in 1:T){
      f[i] <- betaval(mean_F, f_sd)  
    }
  }
  
  if(ft=="Nei"){
    for(i in 1:T){
      f[i] <- mean_F  
    }
  }

# Juvenile summer survival
  phi1 <- matrix(ncol=T, nrow=1) 
  
  if(phi1t=="Ja"){
    for(i in 1:T){
      phi1[i] <- betaval(mean_phi1, phi1_sd)  
    }
  }
  
  if(phi1t=="Nei"){
    for(i in 1:T){
      phi1[i] <- mean_phi1  
    }
  }
  
  ############################################################
  ## HARVEST OFF-TAKE
  
  h <- matrix(ncol=T, nrow=6)
  h[1,] <- sample(size=T, seq(0.2,0.3, by=0.01), replace=T)
  h[2,] <- sample(size=T, seq(0.15,0.30, by=0.01), replace=T)
  h[3,] <- sample(size=T, seq(0.1,0.15, by=0.01), replace=T)
  h[4,] <- sample(size=T, seq(0.2,0.3, by=0.01), replace=T)
  h[5,] <- sample(size=T, seq(0.05,0.1, by=0.01), replace=T)
  h[6,] <- sample(size=T, seq(0.1,0.40, by=0.01), replace=T)
  
  ############################################################
  ## DERIVED HARVEST NUMBERS
  
  H <- matrix(ncol=T, nrow=6)
  
  #############################################################
  ### INITIAL POPULATION SIZE AND STRUCTURE - FOR X
  
  N[1,1] <- 120*m
  N[2,1] <- 120*m
  N[3,1] <- 90*m
  N[4,1] <- 80*m
  N[5,1] <- 150*m
  N[6,1] <- 180*m
  
  
  
  
  for (t in 1:(T-1)){
    ###########################################################
    # STATE PROCESS;
    # PRE-HARVEST POPULATION VECTORS IN T+1  
    
    N[3,t+1] <- rbinom(1, round(N[1,t]*(1-h[1,t])), PHI3)
    N[4,t+1] <- rbinom(1, round(N[2,t]*(1-h[2,t])), PHI3)
    
    N[5,t+1] <- rbinom(1, round((N[3,t]*(1-h[3,t]))+(N[5,t]*(1-h[5,t]))), PHI3)
    N[6,t+1] <- rbinom(1, round((N[4,t]*(1-h[4,t]))+(N[6,t]*(1-h[6,t]))), PHI3)
    
    N[1,t+1] <- rbinom(1, N[5, t+1], phi1[t]*f[t]/2)   
    N[2,t+1] <- rbinom(1, N[5, t+1], phi1[t]*f[t]/2)   
    
    N_juv[1, t+1] <- round(N[5,t+1]*(f[t]/2))                    # Female calves in spring
    N_juv[2, t+1] <- round(N[5,t+1]*(f[t]/2))                    # Male calves in spring
    
    
    #############################################################
    # POST-HARVEST POPULATION VECTORS IN T+1
    X[,t] <- round(N[,t]*(1-h[,t]))
    
    #############################################################
    # DERIVED HARVEST NUMBERS  
    H[,t] <- round(N[,t]*h[,t])
    X_tot[t] <- sum(X[,t])
    N_tot[t] <- sum(N[,t])
    
    if(error.dist=="pois"){
    N_obs[t] <- rpois(1, N_tot[t])
    }
    
    if(error.dist=="bin"){
      N_obs[t] <- rbinom(1, N_tot[t], p.obs)
    }
    
    if(error.dist=="count"){
      N_obs[t] <- N_tot[t]
    }
    
    #############################################################
    # DERIVED OBSERVATIONS
    # Calf surveys in spring
    J[t] <- rbinom(1, round((N[1,t]+N[2, t])/phi1[t]), pa[t])
    SU[t] <- rbinom(1, (N[3, t]+N[4, t]+N[5, t]), pa[t]*bias1)
    
    ##############################################################
    # Structural surveys - post harvest
    C0[t] <- rbinom(1, (X[1,t]+X[2, t]), pb[t])
    Cf[t] <- rbinom(1, (X[3,t]+X[5, t]), pb[t]*bias2)
    Cm[t] <- rbinom(1, (X[4,t]+X[6, t]), pb[t]*bias3)  
    
  }
  
  
  out <- list(X, N, N_juv, h, H, X_tot, N_tot, N_obs, J, SU, C0, Cf, Cm, pa, pb, f, phi1, PHI3)
  names(out) <- c("X", "N", "N_juv", "h", "H", "X_tot", "N_tot","N_obs", "J", "SU", "C0", "Cf", "Cm", "p1", "p2", "f", "phi1", "phi3")
  out
  
}
