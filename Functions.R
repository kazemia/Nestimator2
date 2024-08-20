#This file includes all the functions used in NESTIMATOR IV analysis.
#Run the whole file to load the functions before doing any experiments.

library(tidyverse)
library(combinat)
library(stringr)
library(lpSolve)
library(doParallel)
library(dplyr)

MakeP_Z <- function(data, Z_column, T_column){
  #Function: MakeP_Z
  #Purpose: Calculate P(Treatment|Instrument)
  #Input:
  #       data: is the original data set as data frame
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #Returns: A data frame of probabilities with a column for every value of the
  #         treatment and a row for every value of the instrument.
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

MakeQ_Z <- function(data, Z_column, T_column, Y_column){
  #Function: MakeQ_Z
  #Purpose: Calculate E(Response*1[Treatment = t]|Instrument)
  #Input:
  #       data: is the original data set as data frame
  #       Z_column: is the name of the column in data for the instrument
  #       T_column: is the name of the column in data for the treatment
  #       Y_column: is the name of the column in data for the response
  #Returns: A data frame of expectations with a column for every value of the
  #         treatment and a row for every value of the instrument.
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

MakeKB <- function(R, TLevels, tolerance, balanced = FALSE, weights = NULL){
  #Function: MakeKB
  #Purpose: Making binary indicator matrices and K
  #Input:
  #       R: Output of MakeR
  #       TLevels: A vector of possible values for the treatment
  #       tolerance: The number of decimals of accuracy.
  #Returns: A list of lists of matrices. For every treatment, it calculates the 
  #         indicator matrix B_t, its psudo-inverse B_t^+ and K_t = I - B_tB_t^+
  KB_t <- function(t, R, tolerance, balanced, weights){
    B_t <- 1*(R==t)
    if(balanced){
      for(c in colnames(B_t)){
        B_t[,c] <- B_t[,c] * weights
      }
    }
    B_t_i <- MASS::ginv(B_t)
    temp <- B_t_i %*% B_t
    K = round(diag(ncol(R)) - temp, tolerance)
    list(K = K,
         B_t = B_t,
         B_t_i = B_t_i)
  }
  KB <- lapply(TLevels, KB_t, R, tolerance, balanced, weights)
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
    mylp <- lp('max', rep(0,numcols), rbind(K_t, K_t),
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

CatSimulator <- function(n, Z, A, TLevels, VLevels, VP, VA, VY, TY, T_decider){
  #Function: CatSimulator
  #Purpose: Simulating data with more than 2 treatments
  #Input:
  #       n: is the total number of patients
  #       Z: is a list of ordered vectors of treatments
  #       A: a list of adherence sets
  #       TLevels: a vector of all possible treatments
  #       VLevels: a vector of levels of confounders
  #       VP: a vector of probability of being in each confounder level 
  #           (sums to 1)
  #       VA: a list of vecotr of probability of A for each level of V 
  #           (each element sums to 1)
  #       VY: a vector of effect of each level of V on Y
  #       TY: a vector of effects of T on Y
  #       T_decider: the choice function
  #OPS: max(VY) + max(TY) <= 1
  #Returns: A simulated data set in data frame format with columns
  #         Z for the instrument, T for the treatment and Y for the response.
  V_obs <- sample(VLevels, n, replace = TRUE, prob = VP)
  A_obs <- A[unlist(lapply(VA[as.character(V_obs)],
                           function(x) 
                             sample.int(length(A), size = 1, prob = x)))]
  
  names(A_obs) <- 1:n
  Z_obs_ind <- sample.int(length(Z), size = n, replace = TRUE)
  Z_obs <- Z[Z_obs_ind]
  names(Z_obs) <- 1:n
  T_obs <- sapply(1:n, function(x) T_decider(Z_obs[[x]], A_obs[[x]]))
  YP <- TY[match(T_obs, TLevels)] + VY[match(V_obs,VLevels)]
  Y_obs <- rbinom(rep(n,n),1,YP)
  return(data.frame(Z = Z_obs_ind, T = T_obs, Y = Y_obs))
}

CatSimulator2 <- function(n, Z, TLevels, VT, VP, VY, TY, YPintercept, VReturn, ZT = "exp"){
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
  YP <- TY[match(T_obs, TLevels)] + (V_obs %*% VY) + YPintercept
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

BSCICalculator <- function(n, data, Z_column, T_column, Y_column, TLevels, Z,
                           alpha = 0.05, n_cores = 14, Cap = FALSE, balanced = FALSE,
                           MakeP_Z. = MakeP_Z, MakeQ_Z. = MakeQ_Z,
                           GenerateA. = GenerateA, T_decider. = T_decider, 
                           MakeR. = MakeR, MakeKB. = MakeKB,
                           KbSolver. = KbSolver,
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
                                  TLevels, Z, balanced, MakeP_Z.. = MakeP_Z.,
                                  MakeQ_Z.. = MakeQ_Z.,
                                  GenerateA.. = GenerateA.,
                                  T_decider.. = T_decider., 
                                  MakeR.. = MakeR., MakeKB.. = MakeKB.,
                                  KbSolver.. = KbSolver.,
                                  P_SigmaIdentifier.. = P_SigmaIdentifier.,
                                  LATEIdentifier.. = LATEIdentifier.){
    data <- data[sample(nrow(data), nrow(data), replace=T),]
    Z_weights <- data %>% group_by(Z_value) %>% summarise(N = n())
    Z_weights$Prob <- Z_weights$N/sum(Z_weights$N)
    P_Z <- MakeP_Z..(data, Z_column, T_column)
    Q_Z <- MakeQ_Z..(data, Z_column, T_column, Y_column)
    if(balanced){
      for(t in TLevels){
        P_Z[[t]] <- P_Z[[t]] * Z_weights$Prob
        Q_Z[[t]] <- Q_Z[[t]] * Z_weights$Prob
      }
    }
    A <- GenerateA..(TLevels)
    R <- MakeR..(A, Z, T_decider..)
    KB <- MakeKB..(R, TLevels, 5, balanced, Z_weights$Prob)
    b <- KbSolver..(KB,4)
    P_Sigma <- P_SigmaIdentifier..(P_Z, KB, b)
    Contrasts <- LATEIdentifier..(Q_Z, KB, b, P_Sigma)
    return(unlist(Contrasts))
  }
  myCluster <- makeCluster(n_cores)
  registerDoParallel(myCluster)
  registerDoRNG(1234)
  boot_b <- foreach(b=idiv(n, chunks=getDoParWorkers()), .combine = "cbind",
                    .packages = c("dplyr", "tidyverse", "stringr",
                                  "lpSolve")) %dopar% {
                                    sapply(1:b, StatisticCalculator, data,
                                           Z_column, T_column,Y_column,
                                           TLevels, Z, balanced)
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
                            Z_column, T_column, Y_column,
                            tolerance, Pi_index = 1){
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
  #       tolerance: The number of decimals of accuracy
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
