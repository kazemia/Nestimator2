#This file includes all the functions used in NESTIMATOR IV analysis.
#Run the whole file to load the functions before doing any experiments.

library(tidyverse)
library(combinat)
library(stringr)
library(lpSolve)
library(doParallel)
library(dplyr)
library(doRNG)

MakeP_Z <- function(data, Z_column, T_column, C_columns = NULL, 
                    parametric = FALSE){
  #Function: MakeP_Z
  #Purpose: Calculate P(Treatment|Instrument)
  #Input:
  #       data: is the original data set as data frame
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       C_columns: is a vector of the names of the columns in the data to 
  #                  adjust for
  #       parametric: is an indicator saying whether or not the estimation 
  #                   should be done parametrically with multinomial regression
  #Returns: A data frame of probabilities with a column for every value of the
  #         treatment and a row for every value of the instrument.
  #Note: Adjustment for C_columns is only available with parametric modelling
  if(parametric){
    my_formula <- paste(T_column, "~", Z_column)
    if(!is.null(C_columns)){
      for (c in C_columns) {
        my_formula <- paste(my_formula, "+", c)
      }
    }
    data[[Z_column]] <- as.factor(data[[Z_column]])
    wide_effects <- my_formula %>% as.formula() %>% 
      nnet::multinom(data = data) %>% 
      marginaleffects::avg_predictions(variables = Z_column) %>%
      select(all_of(c("group", "estimate", Z_column))) %>%
      pivot_wider(names_from = group, values_from = estimate) %>% 
      arrange(by = eval(parse(text = Z_column)))
    wide_effects[[Z_column]] <- as.numeric(wide_effects[[Z_column]])
    return(wide_effects)
    
  }else{
    P_Z <- data %>% group_by(eval(parse(text = Z_column)), eval(parse(text = T_column))) %>%
      summarise(P = n(), .groups = "keep") %>%
      as.data.frame()
    colnames(P_Z) <- c(Z_column, T_column, "P")
    P_Z <- P_Z %>% reshape(direction = "wide", idvar = Z_column, timevar = T_column)
    P_Z[is.na(P_Z)] <- 0
    P_Z[,2:ncol(P_Z)] <- P_Z[,2:ncol(P_Z)] %>% apply(1, function(x) x/sum(x)) %>% t()
    colnames(P_Z)[2:ncol(P_Z)] <- substring(colnames(P_Z)[2:ncol(P_Z)], 3)
    rownames(P_Z) <- 1:nrow(P_Z)
    return(P_Z)
  }
}

MakeQ_Z <- function(data, Z_column, T_column, Y_column, C_columns = NULL,
                    parametric = FALSE, family = NULL, P_Z = NULL){
  #Function: MakeQ_Z
  #Purpose: Calculate E(Response*1[Treatment = t]|Instrument)
  #Input:
  #       data: is the original data set as data frame
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       Y_column: is the name of the column in data for the response
  #       C_columns: is a vector of the names of the columns in the data to 
  #                  adjust for
  #       parametric: is an indicator saying whether or not the estimation 
  #                   should be done parametrically with glm
  #       family: family input in the glm regression model
  #       P_Z: The P_Z to be used to multiply with the marginal effects.
  #            Must be passed for parametric modelling.
  #Returns: A data frame of expectations with a column for every value of the
  #         treatment and a row for every value of the instrument.
  #Note: Adjustment for C_columns is only available with parametric modelling
  if(parametric){
    my_formula <- paste0(Y_column, "~", Z_column, "*",T_column)
    if(!is.null(C_columns)){
      for (c in C_columns) {
        my_formula <- paste(my_formula, "+", c)
      }
    }
    data[[Z_column]] <- as.factor(data[[Z_column]])
    long_effects <- my_formula %>% as.formula() %>% 
      glm(data = data, family = family) %>% 
      marginaleffects::avg_predictions(variables = c(Z_column, T_column)) %>%
      select(all_of(c(T_column, "estimate", Z_column)))
    colnames(long_effects) <- c("group", "estimate", Z_column)
    wide_effects <- long_effects %>% 
      pivot_wider(names_from = group, values_from = estimate) %>%
      select(all_of(colnames(P_Z))) %>% 
      arrange(by = eval(parse(text = Z_column)))
    wide_effects[[Z_column]] <- as.numeric(wide_effects[[Z_column]])
    Z_col <- wide_effects[[Z_column]]
    result <- (wide_effects %>% select(-all_of(Z_column)))*(P_Z %>% select(-all_of(Z_column)))
    result[[Z_column]] <- Z_col
    return(result %>% select(all_of(colnames(P_Z))))
    
  }else{
    Q_Z <- data %>% group_by(eval(parse(text = Z_column)), eval(parse(text = T_column))) %>%
      summarise(Y = mean(eval(parse(text = Y_column))),
                Weight = n(), .groups = "keep") %>% 
      as.data.frame()
    colnames(Q_Z) <- c(Z_column, T_column, "Y", "Weight")
    Q_Z <- Q_Z %>% group_by(eval(parse(text = Z_column))) %>%
      mutate(Denominator = sum(Weight)) %>%
      mutate(Y = Y*Weight/Denominator) %>% as.data.frame()
    Q_Z <- Q_Z %>% select(all_of(c(Z_column, T_column, "Y")))
    Q_Z <- Q_Z %>% reshape(direction = "wide", idvar = Z_column, timevar = T_column)
    colnames(Q_Z)[2:ncol(Q_Z)] <- substring(colnames(Q_Z)[2:ncol(Q_Z)], 3)
    rownames(Q_Z) <- 1:nrow(Q_Z)
    Q_Z[is.na(Q_Z)] <- 0
    return(Q_Z)
  }
}

MakeV_Z <- function(data, Z_column, T_column, Y_column, C_columns = NULL,
                    parametric = FALSE, family = NULL, P_Z = NULL){
  #Function: MakeV_Z
  #Purpose: Calculate Var(Response*1[Treatment = t]|Instrument)
  #Input:
  #       data: is the original data set as data frame
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       Y_column: is the name of the column in data for the response
  #       C_columns: is a vector of the names of the columns in the data to 
  #                  adjust for
  #       parametric: is an indicator saying whether or not the estimation 
  #                   should be done parametrically with glm
  #       family: family input in the glm regression model
  #       P_Z: The P_Z to be used to multiply with the marginal effects.
  #            Must be passed for parametric modelling.
  #Returns: A data frame of expectations with a column for every value of the
  #         treatment and a row for every value of the instrument.
  #Note: Adjustment for C_columns is only available with parametric modelling
  if(parametric){
    return("Error: parametric variance estimation is not developed yet")
    
  }else{
    V_Z <- data %>% group_by(eval(parse(text = Z_column)), eval(parse(text = T_column))) %>%
      summarise(Y = var(eval(parse(text = Y_column))),
                Weight = n(), .groups = "keep") %>% 
      as.data.frame()
    colnames(V_Z) <- c(Z_column, T_column, "Y", "Weight")
    V_Z <- V_Z %>% group_by(eval(parse(text = Z_column))) %>%
      mutate(Denominator = sum(Weight)) %>%
      mutate(Y = Y*((Weight/Denominator)^2)) %>% as.data.frame()
    V_Z <- V_Z %>% select(all_of(c(Z_column, T_column, "Y")))
    V_Z <- V_Z %>% reshape(direction = "wide", idvar = Z_column, timevar = T_column)
    colnames(V_Z)[2:ncol(V_Z)] <- substring(colnames(V_Z)[2:ncol(V_Z)], 3)
    rownames(V_Z) <- 1:nrow(V_Z)
    V_Z[is.na(V_Z)] <- 0
    return(V_Z)
  }
}

GenerateA <- function(TLevels){
  #Function: GenerateA
  #Purpose: Generate a list of all possible adherence sets
  #Input:
  #       TLevels: A vector of possible values for the treatment
  #Returns: A list of all possible adherence sets
  A <- unlist(lapply(1:length(TLevels),
                     combn,
                     x = TLevels,
                     simplify = FALSE),
              recursive = FALSE)
  names(A) <- paste0("A", 1:length(A))
  return(A)
}

T_decider <- function(Z_value, A_value){
  #Function: T_decider
  #Purpose: The choice function, determining Treatment|Instrument,Adherence set
  #Input:
  #       Z_value: A vector of all treatments sorted 
  #                 in decreasing order of encouragement
  #       A_value: A vector of the treatments that the agent adheres to
  #Returns: The treatment in A_value that is the most encouraged by Z_value
  return(A_value[which.min(match(A_value, Z_value))])
}

MakeR <- function(A, Z, T_decider = T_decider){
  #Function: MakeR
  #Purpose: Making a response matrix
  #Input:
  #       A: A list of all possible adherence sets, i.e. output of GenerateA
  #       Z: A list of all available value of the instrument
  #Returns: A data frame containing Treatment|Instrument, Adherence set,
  #         with a row for every value of the instrument and
  #         a column for every value of the adherence set.
  outer(Z,A, Vectorize(T_decider)) %>% as.data.frame() %>% return()
}

MakeKB <- function(R, TLevels, tolerance){
  #Function: MakeKB
  #Purpose: Making binary indicator matrices and K
  #Input:
  #       R: Output of MakeR
  #       TLevels: A vector of possible values for the treatment
  #       tolerance: The number of decimals of accuracy.
  #Returns: A list of lists of matrices. For every treatment, it calculates the 
  #         indicator matrix B_t, its psudo-inverse B_t^+ and K_t = I - B_tB_t^+
  KB_t <- function(t, R, tolerance){
    B_t <- 1*(R==t)
    B_t_i <- MASS::ginv(B_t)
    temp <- B_t_i %*% B_t
    K = round(diag(ncol(R)) - temp, tolerance)
    list(K = K,
         B_t = B_t,
         B_t_i = B_t_i)
  }
  KB <- lapply(TLevels, KB_t, R, tolerance)
  names(KB) <- TLevels
  return(KB)
}

KbSolver <- function(KB, tolerance){
  #Function: KbSolver
  #Purpose: Solving the equation bK_t = bK_t' = 0 for every pair of t and t'
  #Input:
  #       KB: Output of MakeKB
  #       tolerance: The number of decimals of accuracy.
  #                  It must be smaller than tolerance for makeKB.
  #Returns: A list of matrices. For every pair of treatments, it calculates
  #         all binary non-zero b's such that bK_t = bK_t' = 0 and r-binds them
  K <- lapply(KB, function(x) x[["K"]])
  KK <- list()
  for(t1_ind in 1:(length(K)-1)){
    for(t2_ind in (t1_ind+1):length(K)){
      KK[[paste(names(K)[t1_ind], names(K)[t2_ind], sep = "_")]] <- K[c(t1_ind, t2_ind)]
    }
  }
  names(K) <- names(KB) #Should be removed?
  K_tbSolver <- function(KK_t){
    K_t <- rbind(KK_t[[1]], KK_t[[2]])
    VN <- colnames(K_t)
    numVars <- ncol(K_t)
    freeVars <- which(colSums(K_t == 0) == (2*numVars))
    fixedVars <- which(rowSums(K_t != 0) == 1)
    fixedVars[fixedVars > numVars] <- fixedVars[fixedVars > numVars] - numVars
    fixedVars <- unique(fixedVars)
    ntVars <- (1:numVars)[-c(freeVars, fixedVars)]
    K_t_backUp <- K_t
    K_t <- rbind(KK_t[[1]][-c(freeVars, fixedVars), -c(freeVars, fixedVars)], 
                 KK_t[[2]][-c(freeVars, fixedVars), -c(freeVars, fixedVars)])
    numcols <- ncol(K_t)
    mylp <- lpSolve::lp('max', rep(0,numcols), rbind(K_t, K_t),
                        c(rep("<", numcols), rep(">", numcols)),
                        c(rep(10^(-tolerance), numcols), rep(-10^(-tolerance), numcols)),
                        all.bin = TRUE, num.bin.solns=((2^numcols) + 1), use.rw = FALSE)
    numsols <- mylp$num.bin.solns
    solutions <- matrix(head(mylp$solution, numcols*numsols), nrow=numsols, byrow=TRUE)
    allSolutions <- matrix(0, nrow = nrow(solutions), ncol = numVars)
    allSolutions[, ntVars] <- solutions
    for (fv in freeVars) {
      temp <- allSolutions
      temp[, fv] <- 1
      allSolutions <- rbind(allSolutions, temp)
    }
    solutions <- allSolutions
    solutions <- solutions[(rowSums(solutions)>0) & 
                             (colSums(abs(K_t_backUp %*% t(solutions)) < 10^(-tolerance)) == (2*numVars)),]
    if(is.null(nrow(solutions))){
      names(solutions) <- VN
    }else{
      colnames(solutions) <- VN
      solutions <- unique(solutions)
    }
    return(solutions)
  }
  lapply(KK, K_tbSolver) %>% return()
}

PiIdentifier <- function(b){
  #Function: PiIdentifier
  #Purpose: Identifying all sets of adherence sets for which a causal effect is
  #         identifiable
  #Input:
  #       b: Output of KbSolver
  #Returns: A list of lists of vector. For every pair of treatments, it returns
  #         a list of all sets of adherence sets for which LATE is identifiable
  Pi_decider <- function(b_t){
    if(is.null(nrow(b_t))){
      list(names(which(b_t==1)))
    }else{
      apply(b_t, 1, function(x) names(which(x==1)))
    }
  }
  lapply(b, Pi_decider) %>% return()
  
}

P_SigmaIdentifier <- function(P_Z, KB, b){
  #Function: P_SigmaIdentifier
  #Purpose: Identifying the probability of a agent belonging to an adherence set
  #         in every Pi identified by PiIdentifier
  #Input:
  #       P_Z: Output of MakeP_Z
  #       KB: Output of makeKB
  #       b: Output of KbSolver
  #Returns: 3 lists of vectors of probabilities. First list, named E1 is
  #         calculated based on the first treatment in the pair t-t'.
  #         Second list, named E2 is
  #         calculated based on the second treatment in the pair t-t'.
  #         Third list, named WA is the weighted average of E1 and E2.
  #         Each list contains an element for every pair of treatments t-t', 
  #         and every element is a vector of probabilities of an agent
  #         belonging to an adherence set in every Pi.
  PrTIdentifier <- function(b_t, B_t_i, P_Z_t) (b_t %*% (B_t_i %*% P_Z_t))
  
  BP_S1 <- lapply(names(b), 
                  function(t) 
                    PrTIdentifier(b[[t]], KB[[str_split(t, "_")[[1]][1]]]$B_t_i,
                                  P_Z[[str_split(t, "_")[[1]][1]]]))
  names(BP_S1) <- names(b)
  BP_S2 <- lapply(names(b),
                  function(t) 
                    PrTIdentifier(b[[t]], KB[[str_split(t, "_")[[1]][2]]]$B_t_i,
                                  P_Z[[str_split(t, "_")[[1]][2]]]))
  names(BP_S2) <- names(b)
  w_t <- colSums(P_Z)[-1]
  names(w_t) <- colnames(P_Z)[-1]
  BP_S <- lapply(names(b), 
                 function(t) 
                   ((BP_S1[[t]]*w_t[str_split(t, "_")[[1]][1]]) +
                      (BP_S2[[t]]*w_t[str_split(t, "_")[[1]][2]]))/
                   (w_t[str_split(t, "_")[[1]][1]] +
                      w_t[str_split(t, "_")[[1]][2]]))
  names(BP_S) <- names(b)
  return(list(E1 = BP_S1, E2 = BP_S2, WA = BP_S))
}

MaxAIdentifier <- function(TLevels, Z, P_Z, minA, likelyA, 
                           minConnections = 2, 
                           ReliabilityCutOff = 0.1, n_cores = 14, upwardsStep = 4){
  likelyA <- lapply(likelyA, sort)
  minA <- lapply(minA, sort)
  likelyA <- likelyA[!(likelyA %in% minA)]
  nonAT <- likelyA[sapply(likelyA, length) > 1]
  AT <- likelyA[sapply(likelyA, length) == 1]
  maxSize <- length(nonAT)
  Found <- FALSE
  AllSets <- expand.grid(rep(list(c(FALSE, TRUE)), maxSize))
  SetSizes <- rowSums(AllSets)
  Found_func <- function(A_temp_index){
    A_temp <- nonAT[unlist(possibleSets[A_temp_index,])]
    A_temp <- c(A_temp, minA, AT)
    R <- MakeR(A_temp, Z, T_decider)
    KB <- MakeKB(R, TLevels, 4)
    b <- KbSolver(KB, 3)
    P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
    ReliableLatesList <- P_Sigma$WA %>% 
      lapply(function(x) max(x) > ReliabilityCutOff) %>% 
      unlist() %>% which() %>% names()
    degrees <- sapply(TLevels, function(x) sum(grepl(x, ReliableLatesList)))
    return(min(degrees))
  }
  myCluster <- makeCluster(n_cores)
  clusterExport(myCluster, c("MakeR", "T_decider", "MakeKB", "KbSolver", 
                             "P_SigmaIdentifier"))
  clusterExport(myCluster, c("Z", "TLevels", "P_Z", "minA", "AT",
                             "ReliabilityCutOff", "nonAT"), envir = environment())
  invisible(clusterEvalQ(myCluster, {
    library(dplyr)
    library(stringr)
  }))
  while ((maxSize >= 0) & (!Found)) {
    possibleSets <- AllSets[SetSizes == maxSize,]
    clusterExport(myCluster, c("possibleSets"), envir = environment())
    degrees_vec <- parSapply(myCluster, 1:nrow(possibleSets), Found_func)
    Found_vec <- (degrees_vec >= minConnections)
    Found <- (sum(Found_vec) > 0)
    maxSize <- maxSize - 1
  }
  #ADD A FOR LOOP HERE FOR MORE THAN 1 FOUND SET
  minA <- c(minA, nonAT[unlist(possibleSets[Found_vec,])])
  
  fullA <- GenerateA(TLevels)
  fullA <- lapply(fullA, sort)
  minA <- lapply(minA, sort)
  fullA <- fullA[!(fullA %in% minA)]
  nonAT <- fullA[sapply(fullA, length) > 1]
  maxSize <- max(upwardsStep, length(nonAT))
  AllSets <- expand.grid(rep(list(c(FALSE, TRUE)), length(nonAT)))
  SetSizes <- rowSums(AllSets)
  clusterExport(myCluster, c("minA", "nonAT"), envir = environment())
  for (s in maxSize:1) {
    possibleSets <- AllSets[SetSizes == s,]
    clusterExport(myCluster, c("possibleSets"), envir = environment())
    degrees_vec <- parSapply(myCluster, 1:nrow(possibleSets), Found_func)
    Found_vec <- (degrees_vec >= minConnections)
    Found <- (sum(Found_vec) > 0)
    if (Found){      
      minA <- c(minA, nonAT[unlist(possibleSets[Found_vec,])])
      nonAT <- nonAT[!(nonAT %in% minA)]
      AllSets <- expand.grid(rep(list(c(FALSE, TRUE)), length(nonAT)))
      SetSizes <- rowSums(AllSets)
      clusterExport(myCluster, c("minA", "nonAT"), envir = environment())
    }
  }
  s <- 16
  while((!Found) & (s>0)){
    possibleSets <- AllSets[SetSizes == s,]
    clusterExport(myCluster, c("possibleSets"), envir = environment())
    degrees_vec <- parSapply(myCluster, 1:nrow(possibleSets), Found_func)
    Found_vec <- (degrees_vec >= minConnections)
    Found <- (sum(Found_vec) > 0)
    s <- s-1
  }
  
  stopCluster(myCluster)
  
  if(Found){
    c(minA, AT) %>% return()
  }else{
    return("We found no sets.")
  }
}





LATEIdentifier <- function(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = FALSE){
  #Function: LATEIdentifier
  #Purpose: Identifying the local average treatment effect for every Pi
  #Input:
  #       P_Z: Output of MakeP_Z
  #       KB: Output of makeKB
  #       b: Output of KbSolver
  #       P_Sigma: Output of P_SigmaIdentifier
  #       RR: Whether or not a risk ratio should be returned 
  #           (if FALSE the treatment effect will be returned)
  #       AverageProb: Wheter or not the average prob should be used
  #Returns: A lists of vectors of expectations. The list contains an element
  #         for every pair of treatments t-t', and every element is a vector
  #         of expected difference (or riskt ratio) in treatment effect for 
  #         agents belonging to an adherence set in every Pi.
  PrTIdentifier <- function(b_t, B_1_i, B_2_i, Q_Z_1, Q_Z_2,
                            P_Sigma_t_1, P_Sigma_t_2, P_Sigma_WA,
                            RR. = RR, AverageProb. = AverageProb) {
    if(RR. & AverageProb.){
      Q_Sigma_t <- ((b_t %*% (B_1_i %*% Q_Z_1))/
                      (b_t %*% (B_2_i %*% Q_Z_2)))
    }else if(RR.){
      Q_Sigma_t <- (((b_t %*% (B_1_i %*% Q_Z_1))/P_Sigma_t_1)/
                      ((b_t %*% (B_2_i %*% Q_Z_2))/P_Sigma_t_2))
    }else if(AverageProb.){
      Q_Sigma_t <- (((b_t %*% (B_1_i %*% Q_Z_1))/P_Sigma_WA) - 
                      ((b_t %*% (B_2_i %*% Q_Z_2))/P_Sigma_WA))
    }else{
      Q_Sigma_t <- (((b_t %*% (B_1_i %*% Q_Z_1))/P_Sigma_t_1) - 
                      ((b_t %*% (B_2_i %*% Q_Z_2))/P_Sigma_t_2))
    }
    return(Q_Sigma_t)
  }
  
  Q_Sigma <- lapply(names(b),
                    function(t)
                      PrTIdentifier(b[[t]], 
                                    KB[[str_split(t, "_")[[1]][1]]]$B_t_i, 
                                    KB[[str_split(t, "_")[[1]][2]]]$B_t_i,
                                    Q_Z[[str_split(t, "_")[[1]][1]]],
                                    Q_Z[[str_split(t, "_")[[1]][2]]],
                                    P_Sigma$E1[[t]], P_Sigma$E2[[t]], P_Sigma$WA[[t]]))
  
  names(Q_Sigma) <- names(b)
  return(Q_Sigma)
}

CatSimulator <- function(n, Z, TLevels, VT, VP, VY, TY, VReturn, ZT = "exp"){
  #Function: CatSimulator
  #Purpose: Simulating data with >2 treatments
  #Input:
  #       n: is the total number of patients
  #       Z: A list of all available value of the instrument
  #       TLevels: A vector of possible values for the treatment
  #       VT: is the effect of the confounders on the treatments.
  #           This is a matrix, where the value in row i and column j is the
  #           effect of the j-th confounder on the i-th treatment
  #       VP: is the incidence rates of the confounders. This is a vector
  #           of the length the number of confounders where the value are
  #           between 0 and 1.
  #       VY: is the effect of the confounder on the response
  #       TY: is the treatment effect on the response
  #       VReturn: Is an indicator vector for wheter or not the confounders
  #                are observed.
  #       ZT: is the shape of the effect of Z on T
  #Returns: A simulated data set in data frame format with columns
  #         Z for the instrument, T for the treatment, Y for the response and
  #         Vs for the observed confounders.
  Z_obs_ind <- sample.int(length(Z), size = n, replace = TRUE)
  Z_obs <- Z[Z_obs_ind]
  names(Z_obs) <- 1:n
  nT <- length(TLevels)
  nV <- length(VP)
  V_obs <- matrix(NA, nrow = n, ncol = length(VP))
  for(v_counter in 1:nV)
    V_obs[,v_counter] <- rbinom(n,1,VP[v_counter])
  VTeffect <- (V_obs %*% t(VT))
  T_probs <- data.frame(matrix(NA, nrow = n, ncol = 0))
  for(t in TLevels)
    T_probs[[t]] <- sapply(Z_obs, function(z) (nT - which(z == t)))
  if(ZT == "exp")
    T_probs <- exp(T_probs)
  T_probs <- T_probs + VTeffect
  T_probs <- t(apply(T_probs, 1, function(r) r/sum(r)))
  T_obs <- apply(T_probs, 1, function(x) sample(TLevels, size = 1, prob = x))
  YP <- V_obs %*% VY
  YP <- TY[match(T_obs, TLevels)] + (V_obs %*% VY)
  YP[YP < 0] = 0
  YP[YP > 1] = 1
  Y_obs <- rbinom(rep(n,n),1,YP)
  out_obj <- data.frame(Z = Z_obs_ind, T = T_obs, Y = Y_obs)
  for(v_counter in 1:nV){
    if(VReturn[v_counter])
      out_obj[[paste0("V", v_counter)]] <- V_obs[,v_counter]
  }
  return(out_obj)
}

BinarySimulator <- function(ZT, VT, VY, TY, n, OR = FALSE){
  #Function: BinarySimulator
  #Purpose: Simulating data with 2 treatments
  #Input:
  #       ZT: is the effect of the instrument of the treatment
  #       VT: is the effect of the confounder on the treatment
  #       VY: is the effect of the confounder on the response
  #       TY: is the treatment effect on the response
  #       n: is the total number of patients
  #       OR: Is an indicator of the effects being exponential or linear.
  #Returns: A simulated data set in data frame format with columns
  #         Z for the instrument, T for the treatment and Y for the response.
  V <- rbinom(n,1,0.5)
  Z <- rbinom(n,1,0.5)
  if(OR){
    T <- rbinom(rep(n,n),1,exp((VT*V) + (ZT*Z))/(1+exp((VT*V) + (ZT*Z))))
    Y <- rbinom(rep(n,n),1,exp((VY*V) + (TY*T))/(1+exp((VY*V) + (TY*T))))
  }else{
    T <- rbinom(rep(n,n),1,(VT*V) + (ZT*Z))
    Y <- rbinom(rep(n,n),1,(VY*V) + (TY*T))
  }
  return(data.frame(Z = Z, T = T, Y = Y))
}

BSCICalculator <- function(n, data, Z_column, T_column, Y_column, KB, b,
                           alpha = 0.05, n_cores = 14, Cap = FALSE, 
                           C_columns = NULL, parametric = FALSE, family = NULL,
                           MakeP_Z. = MakeP_Z, MakeQ_Z. = MakeQ_Z,
                           P_SigmaIdentifier. = P_SigmaIdentifier,
                           LATEIdentifier. = LATEIdentifier){
  #Function: BSCICalculator
  #Purpose: Calculating confidence intervals with bootstrapping
  #Input:
  #       n: is the number of bootstraps
  #       data: is the original data set
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       Y_column: is the name of the column in data for the response
  #       TLevels: A vector of possible values for the treatment
  #       Alpha: P-value for two sided test
  #       n_cores: The number of CPU cores used for bootstrapping
  #       Cap: Whether or not the effects should be capped with max-min
  #Returns: A list with three elements:
  #         Alpha: P-value for two sided test
  #         BSData: The result of bootstrapping
  #         CIs: data frame of confidence intervals
  StatisticCalculator <- function(counter, data, Z_column, T_column, Y_column,
                                  KB, b, C_columns, 
                                  parametric, family,
                                  MakeP_Z.. = MakeP_Z.,
                                  MakeQ_Z.. = MakeQ_Z.,
                                  P_SigmaIdentifier.. = P_SigmaIdentifier.,
                                  LATEIdentifier.. = LATEIdentifier.){
    data <- data[sample(nrow(data), nrow(data), replace=T),]
    P_Z <- MakeP_Z..(data, Z_column, T_column, C_columns, parametric)
    Q_Z <- MakeQ_Z..(data, Z_column, T_column, Y_column, C_columns, parametric, 
                     family, P_Z)
    P_Sigma <- P_SigmaIdentifier..(P_Z, KB, b)
    Contrasts <- LATEIdentifier..(Q_Z, KB, b, P_Sigma)
    return(unlist(Contrasts))
  }
  myCluster <- makeCluster(n_cores)
  registerDoParallel(myCluster)
  registerDoRNG(1234)
  boot_b <- foreach(bb=idiv(n, chunks=getDoParWorkers()), .combine = "cbind",
                    .packages = c("dplyr", "tidyverse", "stringr")) %dopar% {
                                    sapply(1:bb, StatisticCalculator, data,
                                           Z_column, T_column,Y_column, 
                                           KB, b, C_columns, parametric, family)
                    }
  stopCluster(myCluster)
  bootstrap_data <- t(boot_b)
  bootstrap_data <- as.data.frame(bootstrap_data)
  CIs <- data.frame()
  Cap_value <- max(data[[Y_column]]) - min(data[[Y_column]])
  if(Cap){
    for(c in colnames(bootstrap_data)){
      CIs <- rbind(CIs,
                   quantile(bootstrap_data[[c]][
                     (bootstrap_data[[c]] <= Cap_value) &
                       (bootstrap_data[[c]] >= (-Cap_value)) &
                       (!is.na(bootstrap_data[[c]]))],
                            c(alpha/2,(1-(alpha/2)))))
    }
  }else{
    for(c in colnames(bootstrap_data)){
      CIs <- rbind(CIs,
                   quantile(bootstrap_data[[c]],
                            c(alpha/2,(1-(alpha/2))), na.rm = T))
    }
  }
  colnames(CIs) <- c("Lower_bound", "Upper_bound")
  rownames(CIs) <- colnames(bootstrap_data)
  return(list(alpha = alpha, BSData = bootstrap_data, CIs = CIs))
}

NaiveIV <- function(tt, data, Z, T_column, Z_column, Y_column){
  #Function: NaiveIV
  #Purpose: Providing a nive IV estimator
  #Input:
  #       tt: pair of treatments of interest in the form of "t1_t2"
  #       data: is the original data set as data frame
  #       Z: a list of all available values of the instrument
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       Y_column: is the name of the column in data for the response
  #Returns: The LATE ignoring the other treatments than the ones of interest
  t1 <- str_split(tt, "_")[[1]][1]
  t2 <- str_split(tt, "_")[[1]][2]
  new_z <- sapply(Z, function(z) (which(z == t1) < which(z == t2)))
  new_data <- data[data[[T_column]] %in% c(t1,t2),]
  new_data$new_Z <- new_z[new_data[[Z_column]]]
  return((mean(new_data[[Y_column]][new_data$new_Z]) - 
            mean(new_data[[Y_column]][!new_data$new_Z]))/
           (mean(new_data[[T_column]][new_data$new_Z] == t1) - 
              mean(new_data[[T_column]][!new_data$new_Z] == t1)))
}

PseudoPopulator <- function(tt, data, KB, b, P_Sigma,
                            Z_column, T_column, Y_column, Pi_index = 1){
  #Function: PseudoPopulator
  #Purpose: Generate a pseudo-population with weights corresponding to the
  #         non-parametric method
  #Input:
  #       tt: pair of treatments of interest in the form of "t1_t2"
  #       data: is the original data set as data frame  
  #       KB: Output of makeKB
  #       b: Output of KbSolver
  #       P_Sigma: Output of P_SigmaIdentifier
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       Y_column: is the name of the column in data for the response
  #       Pi_index: The population for which LATE should be calculated.
  #                 It should be an integer between 1 and length(Pis$tt), where
  #                 Pis is the output of PiIdentifier
  
  #Returns: A dataframe with weights where regression analysis would yield 
  #         the same result as LATEIdentifier
  
  
  t1 <- str_split(tt, "_")[[1]][1]
  t2 <- str_split(tt, "_")[[1]][2]
  w_z <- nrow(data)/(data %>% group_by(eval(parse(text = Z_column))) %>% 
                         summarise(number = n()))$number
  w_t1 <- sum(data[[T_column]] == t1)/nrow(data)
  w_t2 <- sum(data[[T_column]] == t2)/nrow(data)
  w_beta1 <- (b[[tt]] %*% KB[[t1]]$B_t_i)[Pi_index,]
  w_beta2 <- (b[[tt]] %*% KB[[t2]]$B_t_i)[Pi_index,]
  w_sigma1 <- 1/P_Sigma$E1[[tt]][Pi_index,1]
  w_sigma2 <- 1/P_Sigma$E2[[tt]][Pi_index,1]
  pseudo <- data[data[[T_column]] %in% c(t1, t2),]
  pseudo$w <- w_z[pseudo$Z]
  pseudo$w[pseudo[[T_column]] == t1] <- pseudo$w[pseudo[[T_column]] == t1]*
    w_t1*w_beta1[pseudo[[Z_column]][pseudo[[T_column]] == t1]]*w_sigma1
  pseudo$w[pseudo[[T_column]] == t2] <- pseudo$w[pseudo[[T_column]] == t2]*
    w_t2*w_beta2[pseudo[[Z_column]][pseudo[[T_column]] == t2]]*w_sigma2
  
  return(pseudo)
  
}
