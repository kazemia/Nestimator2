

#Non-linear effect
#Model continous outcome
Denominator_model <- glm(T ~ Z + X, family = "binomial")

Numarator_model <- lm(Y2~Z + X + Z*X + I(X^2))
Confounded_model <- e1071::svm(Y2~ T + X)

#IV and confounded effect estimation
Correction_df <- data.frame(X_r = X_range,
                            tau = sapply(X_range, function(x.) ((alpha*x.*x.) + beta*x. + gamma)),
                            C_E = sapply(X_range, Confounded_estimator),
                            IV_E = sapply(X_range, IV_estimator),
                            Conditional_C_prob = sapply(X_range, C_prob_estimator))

Correction_df <- Correction_df %>% 
  mutate(Marginal_C_prob = dnorm(X_r, mean = mean(X), sd = sd(X)) * 
           Conditional_C_prob,
         high_C = Marginal_C_prob > C_cutoff)


#Corrected model
Corrected_model <- WeightSVM::wsvm(IV_E ~ X_r, 
                      data = Correction_df,
                      weight = Correction_df$Marginal_C_prob)
Correction_df$LIV_E <- predict(Corrected_model, newdata = data.frame(X_r = Correction_df$X_r), type = "response")


ggplot(Correction_df) + 
  geom_line(aes(x = X_r, y = tau,
                color = "tau")) + 
  geom_line(aes(x = X_r, y = IV_E,
                color = "IV estimate")) + 
  geom_line(aes(x = X_r, y = C_E,
                color = "Confounded estimate")) + 
  geom_line(aes(x = X_r, y = LIV_E,
                color = "Corrected estimate")) + xlab("X") + ylab("CATE")

#Effect estimation in the complier region
ggplot(Correction_df %>% filter(high_C)) + 
  geom_line(aes(x = X_r, y = tau,
                color = "tau")) + 
  geom_line(aes(x = X_r, y = IV_E,
                color = "IV estimate")) + 
  geom_line(aes(x = X_r, y = C_E,
                color = "Confounded estimate")) + 
  geom_line(aes(x = X_r, y = LIV_E,
                color = "Corrected estimate")) + xlab("X") + ylab("CATE")

Correction_df <- Correction_df %>% mutate(B = C_E - IV_E)
Bias_model <- lm(B ~ X_r, 
                 data = Correction_df %>% filter(high_C),
                 weights = Marginal_C_prob)
Correction_df$My_E <- sapply(Correction_df$X_r, Corrected_estimator)



ggplot(Correction_df) + 
  geom_line(aes(x = X_r, y = tau,
                color = "tau")) + 
  geom_line(aes(x = X_r, y = My_E,
                color = "Combined estimate")) + 
  geom_line(aes(x = X_r, y = C_E,
                color = "Confounded estimate")) + 
  geom_line(aes(x = X_r, y = LIV_E,
                color = "Corrected IV estimate")) + xlab("X") + ylab("CATE")















#Categorical treatment
library(Hmisc)

TLevels <- c("t1", "t2", "t3")
Z <- list("2010" = c("t1", "t2", "t3"),
          "2011" = c("t2", "t1", "t3"),
          "2012" = c("t2", "t3", "t1"),
          "2013" = c("t3", "t2", "t1")
)

A <- GenerateA(TLevels = TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

UT1 <- -3
UT2 <- -4
XT1 <- 2
XT2 <- 3
ZT <- 4

alpha1 <- 1
beta1 <- 1
alpha2 <- 2
beta2 <- 2

#Simulate
U <- rbinom(n,1,UP)
X <- rnorm(n, mean = X_mean, sd = X_sd)

Z_obs_ind <- sample.int(length(Z), size = n, replace = TRUE)

Rt_1 <- (lapply(Z, function(z) (3 - which(z == "t1"))) %>% unlist())[Z_obs_ind]
Rt_2 <- (lapply(Z, function(z) (3 - which(z == "t2"))) %>% unlist())[Z_obs_ind]
Rt_3 <- (lapply(Z, function(z) (3 - which(z == "t3"))) %>% unlist())[Z_obs_ind]


# covariate matrix
mX = matrix(c(Rt_1, Rt_2, Rt_3, X, U), n, 5)
# coefficients for each choice
TCoef1 = c(ZT, 0, 0, XT1, UT1)
TCoef2 = c(0, ZT, 0, XT2, UT2)
TCoef3 = c(0, 0, ZT, 0, 0)

# vector of probabilities
TProb = cbind(exp(mX%*%TCoef1), exp(mX%*%TCoef2), exp(mX%*%TCoef3))

# multinomial draws
T = TLevels[apply(t(apply(TProb, 1, rmultinom, n = 1, size = 1)), 1,
                  function(x) which(x==1))]

delta <- rnorm(n, 0, delta_sd)
epsilon <- rnorm(n, 0, epsilon_sd)
Y <- alpha1*X*(T == "t1") + beta1*(T == "t1") + 
  alpha2*X*(T == "t2") + beta2*(T == "t2") + UY*U + XY*X + epsilon


Cat_df <- data.frame(Z_value = as.factor(Z_obs_ind), X_value = X, T_value = as.factor(T), Y_value = Y)


Denominator_model <- nnet::multinom(T_value ~ Z_value + X_value, data = Cat_df)

Numarator_model <- lm(Y_value~Z_value + X_value + Z_value*X_value + I(X_value^2), data = Cat_df)

Confounded_model <- lm(Y_value~ T_value + X_value + T_value*X_value, data = Cat_df)

Pi_probability <- function(x_v){
  mX_0 <- matrix(c(lapply(Z, function(z) (3 - which(z == "t1"))) %>% unlist(), 
                   lapply(Z, function(z) (3 - which(z == "t2"))) %>% unlist(), 
                   lapply(Z, function(z) (3 - which(z == "t3"))) %>% unlist(), 
                   rep(x_v, length(Z)), rep(0, length(Z))), 4, 5)
  mX_1 <- matrix(c(lapply(Z, function(z) (3 - which(z == "t1"))) %>% unlist(), 
                   lapply(Z, function(z) (3 - which(z == "t2"))) %>% unlist(), 
                   lapply(Z, function(z) (3 - which(z == "t3"))) %>% unlist(), 
                   rep(x_v, length(Z)), rep(1, length(Z))), 4, 5)
  TP0 <- cbind(exp(mX_0%*%TCoef1), exp(mX_0%*%TCoef2), exp(mX_0%*%TCoef3))
  TP1 <- cbind(exp(mX_1%*%TCoef1), exp(mX_1%*%TCoef2), exp(mX_1%*%TCoef3))
  
  TP0 <- t(apply(TP0, 1, function(r) r/sum(r)))
  TP1 <- t(apply(TP1, 1, function(r) r/sum(r)))
  
  P_Z_T <- cbind(data.frame(Z_value = 1:length(Z)), 
                 as.data.frame(UP*TP1 + (1-UP)*TP0))
  colnames(P_Z_T)[2:4] <- TLevels
  P_Sigma <- P_SigmaIdentifier(P_Z_T, KB, b)
  P_Sigma$WA$t2_t3 %>% return()
}

Pi_prob_estimator <- function(x_v, D_model = Denominator_model){
  P_Z_E <- predict(D_model, 
                   newdata = data.frame(Z_value = c("1", "2", "3", "4"), 
                                        X_value = rep(x_v, length(Z)), 
                                        T_value = c("t1","t2", "t3", "t1")), 
                   type = "probs")
  P_Z_E <- cbind(data.frame(Z_value = 1:length(Z)), P_Z_E)
  colnames(P_Z_E)[2:4] <- TLevels
  P_Sigma <- P_SigmaIdentifier(P_Z_E, KB, b)
  P_Sigma$WA$t2_t3 %>% return()
}

IV_estimator <- function(x_v, N_model = Numarator_model,
                         D_model = Denominator_model, Y_bin = FALSE){
  #Function to estimate the CATE with IV
  P_Z_E <- predict(D_model, 
                   newdata = data.frame(Z_value = c("1", "2", "3", "4"), 
                                        X_value = rep(x_v, length(Z)), 
                                        T_value = c("t1","t2", "t3", "t1")), 
                   type = "probs")
  Q_Z_E <- predict(N_model, 
                   newdata = data.frame(Z_value = c("1", "2", "3", "4"), 
                                        X_value = rep(x_v, length(Z)), 
                                        T_value = c("t1","t2", "t3", "t1")), 
                   type = "response")
  Q_Z_E <- Q_Z_E * P_Z_E
  
  P_Z_E <- cbind(data.frame(Z_value = 1:length(Z)), P_Z_E)
  colnames(P_Z_E)[2:4] <- TLevels
  Q_Z_E <- cbind(data.frame(Z_value = 1:length(Z)), Q_Z_E)
  colnames(Q_Z_E)[2:4] <- TLevels
  
  P_Sigma_E <- P_SigmaIdentifier(P_Z_E, KB, b)
  LATEs <- LATEIdentifier(Q_Z_E, KB, b, P_Sigma_E)
  
  return(LATEs$t2_t3)
}
Confounded_estimator <- function(x_v, C_model = Confounded_model){
  #Function to estimate the CATE with confounded model
  out <- predict(C_model, newdata = data.frame(X_value = c(x_v,x_v), T_value = c("t3","t2")),
                 type = "response")
  return(out[2] - out[1])
}

Corrected_estimator <- function(x_v, C_model = Confounded_model,
                                B_model = Bias_model){
  #Function to estimate the CATE combining the confounded and the IV model
  C_out <- predict(C_model, newdata = data.frame(X_value = c(x_v,x_v), T_value = c("t3","t2")),
                   type = "response")
  B_out <- predict(B_model, newdata = data.frame(X_r = x_v),
                   type = "response")
  return(C_out[2] - C_out[1] - B_out)
}

#IV and confounded effect estimation
Correction_df <- data.frame(X_r = X_range,
                            tau = sapply(X_range, function(x.) ((alpha2*x.) + beta2)),
                            C_E = sapply(X_range, Confounded_estimator),
                            IV_E = sapply(X_range, IV_estimator),
                            Conditional_C_prob = sapply(X_range, Pi_prob_estimator))

Correction_df <- Correction_df %>% 
  mutate(Marginal_C_prob = dnorm(X_r, mean = mean(X), sd = sd(X)) * 
           Conditional_C_prob,
         high_C = Marginal_C_prob > (C_cutoff))

#Compliance estimation
ggplot(Correction_df) + 
  geom_line(aes(x = X_r, y = Conditional_C_prob, 
                color = "Estimated P(S = Pi|X)")) + 
  geom_line(aes(x = X_r, y = sapply(X_range, Pi_probability), 
                color = "True P(S = Pi|X)")) + ylab("Density") + xlab("X")
ggplot(Correction_df) +
  geom_line(aes(x = X_r, y = Marginal_C_prob)) + 
  ylab("Estimated P(S = Pi & X)") + xlab("X")



grid.arrange(ggplot(Correction_df) + 
               geom_line(aes(x = X_r, 
                             y = Marginal_C_prob)) + 
               geom_hline(yintercept = C_cutoff) +
               xlab("X") + ylab("Estimated P(S = C & X)"),
             ggplot(Correction_df) + 
               geom_line(aes(x = X_r, y = tau,
                             color = "tau")) + 
               geom_line(aes(x = X_r, y = IV_E,
                             color = "IV estimate")) + 
               geom_line(aes(x = X_r, y = C_E,
                             color = "Confounded estimate")) + 
               xlab("X") + ylab("CATE") + 
               theme(legend.position="none"),
             ncol = 1, nrow = 2)


#Corrected model
Corrected_model <- lm(IV_E ~ X_r, 
                      data = Correction_df,
                      weights = Marginal_C_prob)
Correction_df$LIV_E <- predict(Corrected_model, newdata = data.frame(X_r = Correction_df$X_r), type = "response")



#Effect estimation in the complier region
ggplot(Correction_df %>% filter(high_C)) + 
  geom_line(aes(x = X_r, y = tau,
                color = "tau")) + 
  geom_line(aes(x = X_r, y = IV_E,
                color = "IV estimate")) + 
  geom_line(aes(x = X_r, y = C_E,
                color = "Confounded estimate")) + xlab("X") + ylab("CATE")

Correction_df <- Correction_df %>% mutate(B = C_E - IV_E)
Bias_model <- lm(B ~ X_r, 
                 data = Correction_df %>% filter(high_C),
                 weights = Marginal_C_prob)
Correction_df$My_E <- sapply(Correction_df$X_r, Corrected_estimator)



ggplot(Correction_df) + 
  geom_line(aes(x = X_r, y = tau,
                color = "tau")) + 
  geom_line(aes(x = X_r, y = My_E,
                color = "Combined estimate")) + 
  geom_line(aes(x = X_r, y = LIV_E,
                color = "Corrected estimate")) + 
  geom_line(aes(x = X_r, y = C_E,
                color = "Confounded estimate")) + xlab("X") + ylab("CATE")




load("C:/Users/kazem/Documents/Data/tdToTSD/td/CompleteData3.Rdata")

CompleteData <- CompleteData %>% 
  mutate(BL_DAS28 = (0.56*sqrt(BL_tjc28)) + 
           (0.28*sqrt(BL_sjc28)) + 
           (0.014*BL_pga) + 
           (0.36*log(BL_crp+1)) + 0.96,
         T_DAS28 = (0.56*sqrt(T_tjc28)) + 
           (0.28*sqrt(T_sjc28)) + 
           (0.014*T_pga) + 
           (0.36*log(T_crp+1)) + 0.96) %>% select(Z_value, trtgrp, T_DAS28, BL_DAS28, age)
CompleteData$trtgrp <- as.character(CompleteData$trtgrp)
CompleteData$trtgrp <- ifelse(CompleteData$trtgrp %in% c("Ada", "Cer", "Gol"), "Hum", CompleteData$trtgrp)
CompleteData$trtgrp<- as.factor(CompleteData$trtgrp)

CC_model <- e1071::svm(T_DAS28 ~ trtgrp + age, data = CompleteData)
age_range <- (150:900)/10
Confounded_estimator <- function(age_v, C_model = CC_model){
  #Function to estimate the CATE with confounded model
  out <- predict(C_model, newdata = data.frame(age = c(age_v,age_v), trtgrp = c("Hum","Inf")),
                 type = "response")
  return(out[2] - out[1])
}

Correction_df <- data.frame(age_r = age_range,
                            C_E = sapply(age_range, Confounded_estimator))

ggplot(Correction_df) + 
  geom_line(aes(x = age_r, y = C_E,
                color = "Confounded estimate")) + xlab("Age") + ylab("CATE")



























#Binary Outcome
Binary_effect <- function(x, alpha. = alpha, beta. = beta, XY. = XY, UY. = UY,
                          UP. = UP){
  #Function to calculate the true conditional risk difference for binary outcome
  return((UP.*(plogis((alpha.*x) + beta. + UY. + (XY.*x)) - 
                 plogis(UY.+(XY.*x)))) + 
           ((1-UP.)*(plogis((alpha.*x) + beta. + (XY.*x)) - plogis(XY.*x))))
}

#Generate new Y
Y <- rbinom(rep(n,n),1,plogis(alpha*X*T + beta*T + UY*U + XY*X + epsilon)) %>% 
  as.logical()


#Plot true risk difference for binary Y
ggplot() + 
  geom_line(aes(x = X_range, y = sapply(X_range, Binary_effect))) + 
  ylab("True risk difference") + xlab("X")

#Make new models for Y

Numarator_model <- glm(Y~Z + X + Z*X, family = "binomial")

Confounded_model <- glm(Y~ T + X + T*X, family = "binomial")


#IV and confounded effect estimation
Correction_df <- data.frame(X_r = X_range,
                            tau = sapply(X_range, Binary_effect),
                            C_E = sapply(X_range, Confounded_estimator),
                            IV_E = sapply(X_range, IV_estimator, Y_bin = TRUE),
                            Conditional_C_prob = sapply(X_range, C_prob_estimator))

Correction_df <- Correction_df %>% 
  mutate(Marginal_C_prob = dnorm(X_r, mean = mean(X), sd = sd(X)) * 
           Conditional_C_prob,
         high_C = Marginal_C_prob > C_cutoff)

grid.arrange(ggplot(Correction_df) + 
               geom_line(aes(x = X_r, 
                             y = Marginal_C_prob)) + 
               geom_hline(yintercept = C_cutoff) +
               xlab("X") + ylab("Estimated P(S = C & X)"),
             ggplot(Correction_df) + 
               geom_line(aes(x = X_r, y = tau,
                             color = "tau")) + 
               geom_line(aes(x = X_r, y = IV_E,
                             color = "IV estimate")) + 
               geom_line(aes(x = X_r, y = C_E,
                             color = "Confounded estimate")) + 
               xlab("X") + ylab("CATE") + 
               theme(legend.position="none"),
             ncol = 1, nrow = 2)
#Effect estimation in the complier region
ggplot(Correction_df %>% filter(high_C)) + 
  geom_line(aes(x = X_r, y = tau,
                color = "True effect")) + 
  geom_line(aes(x = X_r, y = IV_E,
                color = "IV estimate")) + 
  geom_line(aes(x = X_r, y = C_E,
                color = "Confounded estimate")) + xlab("X") + ylab("CATE")

Correction_df <- Correction_df %>% mutate(B = C_E - IV_E)
library(splines)
Bias_model <- lm(B ~ bs(X_r, df = 5), 
                 data = Correction_df %>% filter(high_C),
                 weights = Marginal_C_prob)
Correction_df$My_E <- sapply(Correction_df$X_r, Corrected_estimator)


ggplot(Correction_df) + 
  geom_line(aes(x = X_r, y = tau,
                color = "True effect")) + 
  geom_line(aes(x = X_range, y = My_E,
                color = "Corrected estimate")) + 
  geom_line(aes(x = X_range, y = C_E,
                color = "Confounded estimate")) + xlab("X") + ylab("CATE")






#Good code ends here
iv_model <- ivreg::ivreg(Y~T + X|X+Z)
avg_comparisons(iv_model, newdata = data.frame(X = c(4,4), T = c(0,1), Z = c(0,0)))
























HTE_Simulator <- function(counter, start, stop, step, UP, XP, ZP, UT, XT, ZT, UY, XY, alpha, beta){
  out <- data.frame(n = NA, IV = NA, Confounded = NA, Adj = NA, Interact = NA, HIV = NA)
  for (n in ((start/step):(stop/step))*step) {
    U <- rbinom(n,1,UP)
    X <- rnorm(n, mean = X_mean, sd = X_sd)
    Z <- rbinom(n, 1, ZP)
    delta <- rnorm(n, 0, delta_sd)
    epsilon <- rnorm(n, 0, epsilon_sd)
    TP <- plogis(ZT*Z + XT*X + UT*U + delta)
    T <- rbinom(rep(n,n),1,TP)
    Y <- alpha*X*T + beta*T + UY*U + XY*X + epsilon
    Naive_adjusted_c <- lm(Y ~ T + X)$coefficients
    Naive_interact_c <- lm(Y ~ T + X + X*T)$coefficients
    X_complier_mean1 <- (mean(X[T & Z])*mean(T[Z]) - 
                           mean(X[T & (!Z)])*mean(T[!Z]))/
      (mean(T[Z])-mean(T[!Z]))
    X_complier_mean2 <- (mean(X[(!T) & (!Z)])*mean(!T[!Z]) - 
                           mean(X[(!T) & Z])*mean(!T[Z]))/
      (mean(!T[!Z])-mean(!T[Z]))
    X_complier_mean <- mean(c(X_complier_mean1, X_complier_mean2)) #CONTINUE HERE
    X2 <- (X > X_complier_mean)
    X2 <- X
    Fc_numarator_model <- lm(Y~Z + X2 + Z*X2)
    Fc1_numarator <- 
      avg_comparisons(Fc_numarator_model, 
                      newdata = data.frame(X2 = c(1,1), Z = c(0,1)))
    Fc1_numarator <- Fc1_numarator$estimate[which(Fc1_numarator$term == "Z")]
    
    Fc0_numarator <- 
      avg_comparisons(Fc_numarator_model, 
                      newdata = data.frame(X2 = c(0,0), Z = c(0,1)))
    Fc0_numarator <- Fc0_numarator$estimate[which(Fc0_numarator$term == "Z")]
    
    Fc_denominator_model <- glm(T ~ Z + X2, family = "binomial")
    Fc1_denominator <- 
      avg_comparisons(Fc_denominator_model, 
                      newdata = data.frame(X2 = c(1,1), Z = c(0,1)))
    Fc1_denominator <- Fc1_denominator$estimate[which(Fc1_denominator$term == "Z")]
    
    Fc0_denominator <- 
      avg_comparisons(Fc_denominator_model, 
                      newdata = data.frame(X2 = c(0,0), Z = c(0,1)))
    Fc0_denominator <- Fc0_denominator$estimate[which(Fc0_denominator$term == "Z")]
    
    Fc1 <- Fc1_numarator/Fc1_denominator
    Fc0 <- Fc0_numarator/Fc0_denominator
    out <- rbind(out, c(
      n,
      (mean(Y[Z]) - mean(Y[!Z]))/(mean(T[Z]) - mean(T[!Z])), #Regular IV
      mean(Y[T])-mean(Y[!T]), #Confounded
      Naive_adjusted_c[2], #adjusted
      Naive_interact_c[2] + mean(X)*Naive_interact_c[4], #interact
      Fc0 + (Fc1-Fc0)*mean(X2)
    ))
  }
  return(out[2:nrow(out),])
}

myCluster <- makeCluster(10)
registerDoParallel(myCluster)
registerDoRNG(1234)
sim_set <- foreach(b=idiv(100, chunks=getDoParWorkers()), .combine = "cbind", .packages = c("marginaleffects")) %dopar% {
                                  sapply(1:b, HTE_Simulator, start, stop, step, UP, XP, ZP, UT, XT, ZT, UY, XY, alpha, beta)
                                }
stopCluster(myCluster)

sim_set <- t(sim_set) %>% as.data.frame() %>% tidyr::unnest(cols = c(n, IV, Confounded, Adj, Interact, HIV))

mean_set <- sim_set %>% group_by(n) %>% summarise_all(mean)
sd_set <- sim_set %>% group_by(n) %>% summarise_all(sd)
bias_set <- mean_set - (alpha*X_mean + beta)
bias_set$n <- mean_set$n
# The ATE is alpha*0.5 + beta = 0.375