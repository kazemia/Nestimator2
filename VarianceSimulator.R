setwd("C:/Users/kazem/Dropbox/NESTIMATOR PhD/R Scripts/nestimator2/Nestimator2")
source("Functions.R")
library(doRNG)
library(ggplot2)
library(gridExtra)
library(latex2exp)
TLevels <- c("t1", "t2", "t3")
Z <- list("2010" = c("t1", "t2", "t3"),
          "2011" = c("t2", "t1", "t3"),
          "2012" = c("t2", "t3", "t1"),
          "2013" = c("t3", "t2", "t1")
)

nZ <- length(Z)
names(Z) <- 1:nZ
#Generate all the adherence sets
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)


#Specify arguments again for CatSimulator2
VT = as.matrix(data.frame(V1T = c(exp(3),0,0), V2T = c(0,exp(3),0)))
VP = c(0.5,0.5)
VY = c(0.1, 0.1)
TY = c(0.05,0.1,0.15)
YPintercept = 0.2
Vadj = c(TRUE, FALSE)

#Specify maximum number of observations in simulation, the step 
#and number of simulations pr step
sim_max <- 2000
sim_start <- 250
sim_step <- 250
n_sim <- 500
n_bs <- 500
TT_vec <- c("t1_t2", "t1_t3", "t2_t3")

#Simulate
myCluster <- makeCluster(15)
registerDoParallel(myCluster)
registerDoRNG(123)
sim_res <- list()
for(nP in floor(sim_start/sim_step):floor(sim_max/sim_step)){
  sim_res[[nP]] <- foreach(counter=idiv(n_sim, chunks=getDoParWorkers()), .combine = rbind,
                           .packages = c("dplyr", "tidyverse", "stringr","lpSolve")
  ) %dopar% {
    out <- data.frame(n = NA, LATE = NA, BS_LATE = NA, pseudo_LATE = NA)
    for(counter2 in 1:counter){
      #Simulate the data
      sim_data <- CatSimulator(nP*sim_step, Z, TLevels, VT,
                                VP, VY, TY, c(FALSE, FALSE))
      #Calculate the CIV estimator
      P_Z <- MakeP_Z(sim_data, "Z", "T")
      Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
      P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
      LATE <- LATEIdentifier(Q_Z, KB, b, P_Sigma)$t2_t3[1]
      pseudo_data <- PseudoPopulator("t2_t3", sim_data, KB, b, P_Sigma, "Z", "T", "Y") %>% mutate(Yw = Y*w)
      pseudo_LATE <- mean(pseudo_data$Yw[pseudo_data$T == "t2"]) - mean(pseudo_data$Yw[pseudo_data$T == "t3"])
      LATE_vec <- c()
      pseudo_LATE_vec <- c()
      for (counter3 in 1:n_bs) {
        sim_data_temp <- sim_data[sample(nrow(sim_data), nrow(sim_data), replace=T),]
        P_Z_temp <- MakeP_Z(sim_data_temp, "Z", "T")
        Q_Z_temp <- MakeQ_Z(sim_data_temp, "Z", "T", "Y")
        P_Sigma_temp <- P_SigmaIdentifier(P_Z_temp, KB, b)
        LATE_temp <- LATEIdentifier(Q_Z_temp, KB, b, P_Sigma_temp)$t2_t3[1]
        pseudo_data_temp <- pseudo_data[sample(nrow(pseudo_data), nrow(pseudo_data), replace=T),]
        pseudo_LATE_temp <- mean(pseudo_data_temp$Yw[pseudo_data_temp$T == "t2"]) - mean(pseudo_data_temp$Yw[pseudo_data_temp$T == "t3"])
        LATE_vec <- c(LATE_vec, LATE_temp)
        pseudo_LATE_vec <- c(pseudo_LATE_vec, pseudo_LATE_temp)
      }
      
      out <- rbind(out, data.frame(
        n = rep(nP*sim_step, n_bs),
        LATE = rep(LATE, n_bs),
        BS_LATE = LATE_vec,
        pseudo_LATE = pseudo_LATE_vec
      ))
    }
    out[2:nrow(out),]
  }
}
stopCluster(myCluster)

sim_results <- sim_res[[1]]
for (counter in 2:length(sim_res)) {
  sim_results <- rbind(sim_results, sim_res[[counter]])
}
save.image(file = "VarSim.Rdata")

remove(sim_res)

sim_results$BS_LATE[abs(sim_results$BS_LATE) > 1] <- NA
sim_results$pseudo_LATE[abs(sim_results$pseudo_LATE) > 1] <- NA


var_res <- sim_results %>% group_by(n, LATE) %>% summarise(BS_var = var(BS_LATE, na.rm = T),
                                                           pseudo_var = var(pseudo_LATE, na.rm = T),
                                                           BS_mean = mean(BS_LATE, na.rm = T),
                                                           pseudo_mean = mean(pseudo_LATE, na.rm = T)) %>%
  as.data.frame()
var_res <- var_res %>% group_by(n) %>% summarise(LATE_var = var(LATE, na.rm = T),
                                                 BS_var_mean = mean(BS_var, na.rm = T),
                                                 pseudo_var_mean = mean(pseudo_var, na.rm = T),
                                                 BS_mean_mean = mean(BS_mean, na.rm = T),
                                                 pseudo_mean_mean = mean(pseudo_mean, na.rm = T))
