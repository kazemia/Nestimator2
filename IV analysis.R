load("C:/Users/kazem/Documents/Data/tdToTSD/td/CompleteData.Rdata")
setwd("C:/Users/kazem/Dropbox/NESTIMATOR PhD/R Scripts/nestimator2/Nestimator2")
source("Functions.R")

Z <- list("2010" = c("Inf", "Gol", "Cer", "Eta", "Ada"),
          "2011" = c("Eta", "Inf", "Cer", "Gol", "Ada"),
          "2012" = c("Inf", "Eta", "Cer", "Gol", "Ada"),
          "2013" = c("Cer", "Inf", "Gol", "Eta", "Ada"),
          "2014" = c("Inf", "Cer", "Gol", "Ada", "Eta"),
          "2015" = c("Inf", "Cer", "Gol", "Eta", "Ada"),
          "2016" = c("Inf", "Cer", "Eta", "Gol", "Ada"),
          "2017" = c("Inf", "Eta", "Gol", "Cer", "Ada"),
          "2018" = c("Inf", "Eta", "Cer", "Ada", "Gol"),
          "2019" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2020" = c("Ada", "Inf", "Eta", "Cer", "Gol") ,
          "2021" = c("Ada", "Inf", "Eta", "Cer", "Gol") ,
          "2022" = c("Ada", "Inf", "Eta", "Cer", "Gol") ,
          "2023" = c("Ada", "Inf", "Eta", "Cer", "Gol")  
)
#Find the identifiables
nZ <- length(Z)

A <- list(
  AT_A = c("Ada"),
  AT_C = c("Cer"),
  AT_E = c("Eta"),
  AT_G = c("Gol"),
  AT_I = c("Inf"),
  C_AE = c("Ada", "Eta"),
  C_AG = c("Gol", "Ada"),
  C_AI = c("Inf", "Ada"),
  C_CG = c("Cer", "Gol"),
  NT_EI = c("Ada", "Cer", "Gol"),
  NT_CG = c("Ada", "Eta", "Inf"),
  NT_CE = c("Ada", "Gol", "Inf"),
  NT_I = c("Ada", "Cer", "Gol", "Eta"),
  C = c("Inf", "Eta", "Ada", "Cer", "Gol")
)

TLevels <- levels(CompleteData$trtgrp)

R <- MakeR(A, Z, T_decider)
KB. <- MakeKB(R, TLevels, 4)
b. <- KbSolver(KB., 3)
Pis <- PiIdentifier(b.)
CompleteData$Z_value <- as.numeric(CompleteData$Z_value)
P_Z <- MakeP_Z(CompleteData, "Z_value", "trtgrp")
P_Sigma <- P_SigmaIdentifier(P_Z, KB., b.)
ReliableLatesList. <- P_Sigma$WA %>% 
  lapply(function(x) max(x) > 0.1) %>% 
  unlist() %>% which() %>% names()

ResMaker <- function(CompleteData, KB = KB., b = b., ReliableLatesList = ReliableLatesList., onlyIV = FALSE, ContrastVec = NULL){
  P_Z <- MakeP_Z(CompleteData, "Z_value", "trtgrp")
  Q_Z <- MakeQ_Z(CompleteData, "Z_value", "trtgrp", "Rem")
  P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
  LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)
  
  BSCIs <- BSCICalculator(10000, CompleteData, "Z_value", "trtgrp", "Rem", 
                          KB, b, Cap = TRUE)
  largestSet <- lapply(b, function(x){
    if (length(x) <= ncol(KB[[1]]$K))
      return(1)
    return(which(rowSums(x) == max(rowSums(x))))
    }) %>% unlist()
  ReliableLates <- lapply(1:length(LATEs), function(x) {
    if(length(LATEs[[x]]) == 1)
      return(LATEs[[x]])
    return(LATEs[[x]][largestSet[x]])
    })
  names(ReliableLates) <- names(LATEs)
  ReliableLates <- ReliableLates[ReliableLatesList] %>% 
    unlist() %>% as.data.frame()
  ReliableLates$contrast = rownames(ReliableLates)
  colnames(ReliableLates) <- c("IV_estimate", "contrast")
  rownames(ReliableLates) <- 1:nrow(ReliableLates)
  ReliableCIs <- BSCIs$CIs
  ReliableCIs$contrast <- rownames(ReliableCIs)
  whichCIs <- sapply(ReliableLatesList, function(x){
    if(x %in% ReliableCIs$contrast)
      return(x)
    return(paste0(x, largestSet[x]))
  })
  ReliableCIs <- ReliableCIs %>% filter(contrast %in% whichCIs) %>%
    mutate(contrast = gsub('[0-9]+', '', contrast))
  
  colnames(ReliableCIs) <- c("IV_conf.low", "IV_conf.high", "contrast")
  rownames(ReliableCIs) <- 1:nrow(ReliableCIs)
  IV_comparisons <- ReliableLates %>% left_join(ReliableCIs, by = "contrast")
  
  if(onlyIV){
    for (row in 1:nrow(IV_comparisons)) {
      t1 <- IV_comparisons$contrast[row] %>% substr(1,3)
      t2 <- IV_comparisons$contrast[row] %>% substr(5,7)
      if(paste(t2, "-", t1) %in% ContrastVec){
        IV_comparisons$contrast[row] <- paste(t2, "-", t1)
        IV_comparisons$IV_estimate[row] <- -IV_comparisons$IV_estimate[row]
        temp <- IV_comparisons$IV_conf.low[row]
        IV_comparisons$IV_conf.low[row] <- -IV_comparisons$IV_conf.high[row]
        IV_comparisons$IV_conf.high[row] <- -temp
      }else{
        IV_comparisons$contrast[row] <- paste(t1, "-", t2)
      }
    }
    my_comparisons <- IV_comparisons %>% 
      mutate_if(is.numeric, round, digits = 2) %>% 
      mutate(IV_CI = paste0("(", IV_conf.low, ",", IV_conf.high, ")"))%>%
      select(contrast, IV_estimate, IV_CI)
    return(my_comparisons)
  }else{
    my_base_model <- glm(Rem ~ ., data = CompleteData %>% select(-Z_value), family = "binomial")
    my_base_comparisons <- my_base_model %>% marginaleffects::avg_comparisons(variables = list(trtgrp = "pairwise"), type = "response", newdata = "marginalmeans")
    
    
    for (row in 1:nrow(IV_comparisons)) {
      t1 <- IV_comparisons$contrast[row] %>% substr(1,3)
      t2 <- IV_comparisons$contrast[row] %>% substr(5,7)
      if(paste(t2, "-", t1) %in% my_base_comparisons$contrast){
        IV_comparisons$contrast[row] <- paste(t2, "-", t1)
        IV_comparisons$IV_estimate[row] <- -IV_comparisons$IV_estimate[row]
        temp <- IV_comparisons$IV_conf.low[row]
        IV_comparisons$IV_conf.low[row] <- -IV_comparisons$IV_conf.high[row]
        IV_comparisons$IV_conf.high[row] <- -temp
      }else{
        IV_comparisons$contrast[row] <- paste(t1, "-", t2)
      }
    }
    
    my_comparisons <- my_base_comparisons %>% 
      select(contrast, estimate, conf.low, conf.high) %>% 
      left_join(IV_comparisons, by = "contrast") %>% 
      mutate_if(is.numeric, round, digits = 2) %>%
      mutate(R_CI = paste0("(", conf.low, ",", conf.high, ")"),
             IV_CI = ifelse(!is.na(IV_estimate), paste0("(", IV_conf.low, ",", IV_conf.high, ")"), "NA")) %>% 
      select(contrast, estimate, R_CI, IV_estimate, IV_CI)
    
    return(my_comparisons)
  }
}

med_age <- median(CompleteData$age)
#Make results table for primary and startified analyses
my_comparisons <- CompleteData %>% ResMaker()
my_old_comparisons <- CompleteData %>% filter(age > med_age) %>% ResMaker()
my_young_comparisons <- CompleteData %>% filter(age <= med_age) %>% ResMaker()
my_seropositive_comparisons <- CompleteData %>% filter(anticcp) %>% ResMaker()
my_seronegative_comparisons <- CompleteData %>% filter(!anticcp) %>% ResMaker()

colnames(my_old_comparisons)[2:5] <- paste(colnames(my_old_comparisons)[2:5], "old", sep = "_")
colnames(my_young_comparisons)[2:5] <- paste(colnames(my_young_comparisons)[2:5], "young", sep = "_")
colnames(my_seropositive_comparisons)[2:5] <- paste(colnames(my_seropositive_comparisons)[2:5], "seropositive", sep = "_")
colnames(my_seronegative_comparisons)[2:5] <- paste(colnames(my_seronegative_comparisons)[2:5], "seronegative", sep = "_")

knitr::kable(my_comparisons, "latex")
knitr::kable(cbind(my_young_comparisons, my_old_comparisons %>% select(-contrast)), "latex")
knitr::kable(cbind(my_seronegative_comparisons, my_seropositive_comparisons %>% select(-contrast)), "latex")

#Make tables for secondary and sensitivity analysis for different Pis
foundR <- MakeR(foundA, Z, T_decider)
foundKB <- MakeKB(foundR, TLevels, 4)
foundb <- KbSolver(foundKB, 3)
foundPis <- PiIdentifier(foundb)
foundP_Sigma <- P_SigmaIdentifier(P_Z, foundKB, foundb)
foundReliableLatesList <- foundP_Sigma$WA %>% 
  lapply(function(x) max(x) > 0.1) %>% 
  unlist() %>% which() %>% names()
my_found_comparisons <- CompleteData %>% ResMaker(foundKB, foundb, foundReliableLatesList, onlyIV = TRUE, ContrastVec = my_comparisons$contrast)
knitr::kable(my_found_comparisons, "latex")

fullA <- GenerateA(TLevels)
fullR <- MakeR(fullA, Z, T_decider)
fullKB <- MakeKB(fullR, TLevels, 4)
fullb <- KbSolver(fullKB, 3)
fullPis <- PiIdentifier(fullb)
fullP_Sigma <- P_SigmaIdentifier(P_Z, fullKB, fullb)
fullReliableLatesList <- fullP_Sigma$WA %>% 
  lapply(function(x) max(x) > 0.1) %>% 
  unlist() %>% which() %>% names()
my_full_comparisons <- CompleteData %>% ResMaker(fullKB, fullb, fullReliableLatesList, onlyIV = TRUE, ContrastVec = my_comparisons$contrast)
knitr::kable(my_full_comparisons, "latex")

Pz <- (summary(as.factor(CompleteData$Z_value))/nrow(CompleteData))

#Make table 1 for target trials

SumTbMaker <- function(tt, SumData = CompleteData){
  my_index <- ifelse(tt == "Cer_Gol", 3, 1)
  P_Z2 <- MakeP_Z(SumData, "Z_value", "trtgrp")
  P_Sigma2 <- P_SigmaIdentifier(P_Z2, KB, b)
  pseudo_data <- PseudoPopulator(tt, SumData, KB, b, P_Sigma2, 
                                 "Z_value", "trtgrp", "Rem", my_index)
  my_weights <- pseudo_data$w
  pseudo_data <- pseudo_data %>% 
    mutate(BL_DAS28 = (0.56*sqrt(BL_tjc28)) + 
             (0.28*sqrt(BL_sjc28)) + 
             (0.014*BL_pga) + 
             (0.36*log(BL_crp+1)) + 0.96)
  
  pseudo_data <- pseudo_data[, c(2,3,5:13,30,32,55)]
  
  temp <- pseudo_data$trtgrp
  Summary_table_mat <- pseudo_data %>% 
    mutate(BL_otherDmards = BL_otherDmards | BL_Prednisolone | BL_Sulfasalazine) %>%
    select(-BL_Prednisolone, -BL_Sulfasalazine) %>%
    mutate_if(is.logical, as.integer) %>%
    select(-trtgrp) %>%
    fastDummies::dummy_cols(remove_most_frequent_dummy = TRUE,
                            remove_selected_columns = TRUE) %>% 
    mutate_if(is.integer, function(x) 100*x) %>% as.data.frame()
  
  Summary_table_mat$trtgrp <- temp
  
  Summary_table <- Summary_table_mat %>% 
    mutate_if(is.numeric, function(x, w) x*w, my_weights) %>% 
    group_by(trtgrp) %>% summarise_all(mean) %>% mutate_if(is.numeric, round)
  
  Summary_table[, !(colnames(Summary_table) %in% 
                      c("trtgrp", "age", "BL_diagdur", "BL_DAS28"))] <- 
    Summary_table[, !(colnames(Summary_table) %in% 
                        c("trtgrp", "age", "BL_diagdur", "BL_DAS28"))] %>% 
    mutate_all(function(x) ifelse(x > 100, NA, ifelse(x < 0, NA, x)))
  
  Summary_table <- Summary_table[, c(1,3,11:14,4:7,9:10,8)]
  
  Summary_table$trtgrp <- as.character(Summary_table$trtgrp)
  
  colnames(Summary_table) <- c("Treatment", 
                               "Age (years)", "Male (%)", "Smoker current (%)", 
                               "Smoker occationally (%)", "Smoker previously (%)", 
                               "RF positive (%)", "Anti-CCP positive (%)", 
                               "Previous MTX user (%)",  "Methotrexate at BL (%)", 
                               "Other DMARDs  at BL (%)", "DAS28 CRP", 
                               "Diagnosis duration (years)")
  return(Summary_table)
}



myCluster <- makeCluster(14)
clusterExport(myCluster, c("SumTbMaker", "ReliableLatesList", "CompleteData",
                           "PseudoPopulator", "MakeP_Z", "P_SigmaIdentifier",
                           "KB", "b"))
registerDoParallel(myCluster)
registerDoRNG(1234)
my_sum_list <- foreach(bb=idiv(10000, chunks=getDoParWorkers()),
                  .packages = c("dplyr", "stringr")) %dopar% {
                    l_e <- list()
                    for (counter in 1:bb) {
                      my_sum2 <- lapply(ReliableLatesList, SumTbMaker,
                                        CompleteData[sample(nrow(CompleteData), 
                                                            nrow(CompleteData), 
                                                            replace=T),])
                      names(my_sum2) <- ReliableLatesList
                      l_e[[counter]] <- my_sum2
                    }
                    l_e
                    }
stopCluster(myCluster)

my_sum_biglist <- do.call(c, my_sum_list)
my_sum_list <- list()
my_low_list <- list()
my_high_list <- list()
for(tt in ReliableLatesList){
  my_sum_smallList <- lapply(my_sum_biglist, function(x) x[[tt]])
  my_sum_df <- do.call(rbind,my_sum_smallList)
  my_sum_list[[tt]] <- my_sum_df %>% group_by(Treatment) %>% summarise_all(median, na.rm = TRUE)
  my_low_list[[tt]] <- my_sum_df %>% group_by(Treatment) %>% summarise_all(quantile, na.rm = TRUE, probs = 0.025)
  my_high_list[[tt]] <- my_sum_df %>% group_by(Treatment) %>% summarise_all(quantile, na.rm = TRUE, probs = 0.975)
}
my_sum <- do.call(rbind,my_sum_list)
my_low <- do.call(rbind,my_low_list)
my_high <- do.call(rbind,my_high_list)

Summary_table <- my_sum
for (counter in 2:ncol(my_sum)) {
  Summary_table[[counter]] <- paste0(round(my_sum[[counter]]), " (", round(my_low[[counter]]), ", ", round(my_high[[counter]]), ")")
}
Summary_table <- t(Summary_table)
knitr::kable(Summary_table[,7:10], "latex")

remove(my_sum, my_sum_biglist, my_sum_list, my_sum_smallList, l_e)

#
treeData <- CompleteData %>% select(-Z_value, -Rem)
colnames(treeData)[31:51] <- paste0("comorb", 1:21)
set.seed(2)
myTree <- party::ctree((trtgrp=="Gol") ~ ., data = treeData %>%
                         select(-center, -all_of(paste0("comorb", 1:21))),
             control = party::ctree_control(mincriterion=0.95, 
                                            minsplit=20, 
                                            minbucket=10))
plot(myTree)




adj_cols <- c("center", "sex", "age", "smoker", "rhmfact", "prevmtx",
              "BL_Methotrexate", "BL_tjc28", "BL_sjc28", "BL_pga", "BL_crp", "BL_esr",
              "BL_diagdur", "BL_sympdur")
P_Z_adj <- MakeP_Z(CompleteData, "Z_value", "trtgrp", adj_cols,T)
Q_Z_adj <- MakeQ_Z(CompleteData, "Z_value", "trtgrp", "Rem", adj_cols, T, "binomial", P_Z_adj)
LATEs_adj <- LATEIdentifier(Q_Z_adj, KB, b, P_SigmaIdentifier(P_Z_adj, KB, b))
BSCIs_adj <- BSCICalculator(1000, CompleteData, "Z_value", "trtgrp", "Rem", 
                        KB, b, Cap = TRUE, C_columns = adj_cols, parametric = T,
                        family = "binomial")
save(BSCIs_adj, file = "BSCIs_adj.Rdata")





