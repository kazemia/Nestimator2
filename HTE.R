library(doParallel)
library(doRNG)
library(dplyr)

UP <- 0.5
XP <- 0.5
ZP <- 0.5

UT <- 0.2
XT <- 0.2
ZT <- 0.3

UY <- 0.2
XY <- 0.2
alpha <- 0.25
beta <- 0.25

start <- 500
stop <- 5000
step <- 500


HTE_Simulator <- function(counter, start, stop, step, UP, XP, ZP, UT, XT, ZT, UY, XY, alpha, beta){
  out <- data.frame(n = NA, IV = NA, Confounded = NA, Adj = NA, Interact = NA, HIV = NA)
  for (n in ((start/step):(stop/step))*step) {
    U <- rbinom(n,1,UP)
    X <- rbinom(n,1,XP)
    Z <- rbinom(n, 1, ZP)
    TP <- ZT*Z + XT*X + UT*U
    T <- rbinom(rep(n,n),1,TP)
    YP <- alpha*X*T + beta*T + UY*U + XY*X
    Y <- rbinom(rep(n,n),1,YP)
    Naive_interact_c <- lm(Y ~ T + X + X*T)$coefficients
    Fc1 <- (mean(Y[Z & X]) - mean(Y[(!Z) & X]))/
      (mean(T[Z & X]) - mean(T[(!Z) & X]))
    Fc0 <- (mean(Y[Z & (!X)]) - mean(Y[(!Z) & (!X)]))/
      (mean(T[Z & (!X)]) - mean(T[(!Z) & (!X)]))
    out <- rbind(out, c(
      n,
      (mean(Y[Z]) - mean(Y[!Z]))/(mean(T[Z]) - mean(T[!Z])),
      mean(Y[T])-mean(Y[!T]),
      (mean(Y[T & X])-mean(Y[(!T) & X]))*mean(X) + 
        (mean(Y[T & (!X)])-mean(Y[(!T)  & (!X)]))*(1-mean(X)),
      Naive_interact_c[2] + mean(X)*Naive_interact_c[4],
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
bias_set <- mean_set - (alpha*XP + beta)
bias_set$n <- mean_set$n
# The ATE is alpha*0.5 + beta = 0.375