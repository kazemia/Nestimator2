library(CIMTx)
set.seed(1)
data <- data_sim(
  sample_size = 500, n_trt = 3,
  x = c("rnorm(0, 0.5)",  # x1
        "rbeta(2, .4)",   # x2
        "runif(0, 0.5)",  # x3
        "rweibull(1, 2)", # x4
        "rbinom(1, .4)"),   # x5
  # linear terms in parallel response surfaces
  lp_y = rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3), 
  # nonlinear terms in parallel response surfaces
  nlp_y  = rep(".7*x1*x1  - .1*x2*x3", 3), 
  align = F,# different predictors used in treatment and outcome models
  # linear terms in treatment assignment model
  lp_w = c(".4*x1 + .1*x2  - .1*x4 + .1*x5",   # w = 1
           ".2*x1 + .2*x2  - .2*x4 - .3*x5"),  # w = 2
  # nonlinear terms in treatment assignment model
  nlp_w = c("-.5*x1*x4  - .1*x2*x5", # w = 1
            "-.3*x1*x4 + .2*x2*x5"), # w = 2
  tau = c(-1.5, 0, 1.5), delta = c(0.5, 0.5), psi = 1)

ra_ate_res <- ce_estimate(y = data$y, x = data$covariates, w = data$w, 
                          method = "RA", estimand = "ATE", ndpost = 100)
summary(ra_ate_res)

data2 <- cbind(data$w, data$y, data$covariates)
colnames(data2)[1:2] <- c("w", "y")
data2$w <- as.factor(data2$w)
data2$y <- as.logical(data2$y)

my_model <- glm(y ~ ., data = data2, family = "binomial")
my_model %>% marginaleffects::avg_comparisons(variables = list(w = "pairwise"), type = "response", newdata = "marginalmeans")


iptw_sl_trim_ate_res <- ce_estimate(y = data$y, x = data$covariates , w = data$w, 
                                    method = "IPTW-SL", estimand = "ATE", 
                                    sl_library =  c("SL.glm", "SL.glmnet", "SL.rpart"),
                                    trim_perc = c(0.05,0.95), boot = FALSE, verbose_boot = F)

summary(iptw_sl_trim_ate_res)

bart_res <- ce_estimate(y = data$y, x = data$covariates, w = data$w, method = "BART", 
                        estimand = "ATE", ndpost=100, reference_trt = 1, discard = TRUE)

summary(bart_res)

tmle_res_boot <- ce_estimate(y = data$y, x = data$covariates, w = data$w, nboots = 100, 
                             method = "TMLE", estimand = "ATE", boot = TRUE, 
                             sl_library = c("SL.glm", "SL.glmnet", "SL.rpart"))
summary(tmle_res)
