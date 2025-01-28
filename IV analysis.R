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
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)
CompleteData$Z_value <- as.numeric(CompleteData$Z_value)
P_Z <- MakeP_Z(CompleteData, "Z_value", "trtgrp")
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
ReliableLatesList <- P_Sigma$WA %>% 
  lapply(function(x) max(x) > 0.1) %>% 
  unlist() %>% which() %>% names()

ResMaker <- function(CompleteData, KB = KB, b = b, ReliableLatesList = ReliableLatesList, onlyIV = FALSE){
  P_Z <- MakeP_Z(CompleteData, "Z_value", "trtgrp")
  Q_Z <- MakeQ_Z(CompleteData, "Z_value", "trtgrp", "Rem")
  P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
  LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)
  
  BSCIs <- BSCICalculator(10000, CompleteData, "Z_value", "trtgrp", "Rem", 
                          KB, b, Cap = TRUE)
  largestSet <- lapply(b, function(x){
    if (is.null(dim(x)))
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
    IV_comparisons %>% 
      mutate_if(is.numeric, round, digits = 2) %>% 
      mutate(IV_CI = paste0("(", IV_conf.low, ",", IV_conf.high, ")"))%>%
      select(contrast, IV_estimate, IV_CI) %>% return()
  }
  
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

med_age <- median(CompleteData$age)
#Make results table for primary and sub analyses
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


foundR <- MakeR(foundA, Z, T_decider)
foundKB <- MakeKB(foundR, TLevels, 4)
foundb <- KbSolver(foundKB, 3)
foundPis <- PiIdentifier(foundb)
foundP_Sigma <- P_SigmaIdentifier(P_Z, foundKB, foundb)
foundReliableLatesList <- foundP_Sigma$WA %>% 
  lapply(function(x) max(x) > 0.1) %>% 
  unlist() %>% which() %>% names()

LATEs <- LATEIdentifier(Q_Z, foundKB, foundb, foundP_Sigma)

BSCIs <- BSCICalculator(10000, CompleteData, "Z_value", "trtgrp", "Rem", 
                        foundKB, foundb, Cap = TRUE)


my_found_comparisons <- CompleteData %>% ResMaker(foundKB, foundb)

Pz <- (summary(as.factor(CompleteData$Z_value))/nrow(CompleteData))

#Make table 1 for target trials
myCluster <- makeCluster(length(ReliableLatesList))
registerDoParallel(myCluster)
registerDoRNG(1234)
my_sum <- foreach(tt=ReliableLatesList, .combine = "rbind",
                  .packages = c("dplyr", "stringr")) %dopar% {
  my_index <- ifelse(tt == "Cer_Gol", 3, 1)
  pseudo_data <- PseudoPopulator(tt, CompleteData, KB, b, P_Sigma, 
                                 "Z_value", "trtgrp", "Rem", my_index)
  #t1 <- substr(tt, 1,3)
  #t2 <- substr(tt, 5,7)
  #if(tt == "Cer_Gol"){
  #  headcount1 <- round(P_Sigma$E1[[tt]][3] * min(1, Pz %*% KB[[t1]]$B_t %*% b[[tt]][3,]) * nrow(CompleteData))[1]
  #  headcount2 <- round(P_Sigma$E2[[tt]][3] * min(1, Pz %*% KB[[t2]]$B_t %*% b[[tt]][3,]) * nrow(CompleteData))[1]
  #}else if(tt %in% c("Cer_Inf", "Eta_Gol")){
  #  headcount1 <- round(P_Sigma$E1[[tt]] * min(1, Pz %*% KB[[t1]]$B_t %*% t(b[[tt]])) * nrow(CompleteData))[1]
  #  headcount2 <- round(P_Sigma$E2[[tt]] * min(1, Pz %*% KB[[t2]]$B_t %*% t(b[[tt]])) * nrow(CompleteData))[1]
  #}else{
  #  headcount1 <- round(P_Sigma$E1[[tt]] * min(1, Pz %*% KB[[t1]]$B_t %*% b[[tt]]) * nrow(CompleteData))[1]
  #  headcount2 <- round(P_Sigma$E2[[tt]] * min(1, Pz %*% KB[[t2]]$B_t %*% b[[tt]]) * nrow(CompleteData))[1]
  #}
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
    mutate_all(function(x) ifelse(x > 100, 100, ifelse(x < 0, 0, x)))
  
  Summary_table <- Summary_table[, c(1,3,11:14,4:7,9:10,8)]
  
  Summary_table$trtgrp <- as.character(Summary_table$trtgrp)
  
  #Summary_table$trtgrp[Summary_table$trtgrp == t1] <- 
  #  paste0(Summary_table$trtgrp[Summary_table$trtgrp == t1], " (", headcount1, ")")
  #Summary_table$trtgrp[Summary_table$trtgrp == t2] <- 
  #  paste0(Summary_table$trtgrp[Summary_table$trtgrp == t2], " (", headcount2, ")")
  
  colnames(Summary_table) <- c("Treatment", 
                               "Age (years)", "Male (%)", "Smoker current (%)", 
                               "Smoker occationally (%)", "Smoker previously (%)", 
                               "RF positive (%)", "Anti-CCP positive (%)", 
                               "Previous MTX user (%)",  "Methotrexate at BL (%)", 
                               "Other DMARDs  at BL (%)", "DAS28 CRP", 
                               "Diagnosis duration (years)")
  Summary_table
                  }
stopCluster(myCluster)
Summary_table <- t(my_sum)
knitr::kable(Summary_table, "latex")

CompleteData2 <- CompleteData
colnames(CompleteData2)[33:53] <- paste0("comorb", 1:21)
adj_cols <- colnames(CompleteData2)[4:ncol(CompleteData2)]
P_Z_adj <- MakeP_Z(CompleteData2, "Z_value", "trtgrp", adj_cols,T)
Q_Z_adj <- MakeQ_Z(CompleteData2, "Z_value", "trtgrp", "Rem", adj_cols, T, "binomial", P_Z_adj)
LATEs_adj <- LATEIdentifier(Q_Z_adj, KB, b, P_SigmaIdentifier(P_Z_adj, KB, b))
BSCIs_adj <- BSCICalculator(10000, CompleteData2, "Z_value", "trtgrp", "Rem", 
                        KB, b, Cap = TRUE, C_columns = adj_cols, parametric = T,
                        family = "binomial")





