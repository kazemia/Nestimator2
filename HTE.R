library(doParallel)
library(doRNG)
library(dplyr)

UP <- 0.5
X_mean <- 0.5
X_sd <- 1
ZP <- 0.5
delta_sd <- 1
epsilon_sd <- 1

UT <- 2
XT <- 2
ZT <- 3

UY <- 2
XY <- 2
alpha <- 2.5
beta <- 2.5

start <- 500
stop <- 5000
step <- 500


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
    X <- X + X_complier_mean - mean(X)
    HIV_numarator_model <- lm(Y ~ Z + X + X*Z)
    HIV_denominator_model <- glm(T ~ Z + X + Z*X, family = "binomial")
    predict()
    Fc1 <- (mean(Y[Z & X]) - mean(Y[(!Z) & X]))/
      (mean(T[Z & X]) - mean(T[(!Z) & X]))
    Fc0 <- (mean(Y[Z & (!X)]) - mean(Y[(!Z) & (!X)]))/
      (mean(T[Z & (!X)]) - mean(T[(!Z) & (!X)]))
    out <- rbind(out, c(
      n,
      (mean(Y[Z]) - mean(Y[!Z]))/(mean(T[Z]) - mean(T[!Z])), #Regular IV
      mean(Y[T])-mean(Y[!T]), #Confounded
      Naive_adjusted_c[2], #adjusted
      Naive_interact_c[2] + mean(X)*Naive_interact_c[4], #interact
      Fc0 + (Fc1-Fc0)*mean(X)
    ))
  }
  return(out[2:nrow(out),])
}

myCluster <- makeCluster(10)
registerDoParallel(myCluster)
registerDoRNG(1234)
sim_set <- foreach(b=idiv(1000, chunks=getDoParWorkers()), .combine = "cbind") %dopar% {
                                  sapply(1:b, HTE_Simulator, start, stop, step, UP, XP, ZP, UT, XT, ZT, UY, XY, alpha, beta)
                                }
stopCluster(myCluster)

sim_set <- t(sim_set) %>% as.data.frame() %>% tidyr::unnest()

mean_set <- sim_set %>% group_by(n) %>% summarise_all(mean)
sd_set <- sim_set %>% group_by(n) %>% summarise_all(sd)
bias_set <- mean_set - (alpha*X_mean + beta)
bias_set$n <- mean_set$n
# The ATE is alpha*0.5 + beta = 0.375