library(dplyr)
library(tidyr)
library(kernlab)
library(furrr)
library(drf)
library(readr)


simulate <- function(N, d, e, m, tau, seed = 1, Xtest = NULL) {
  set.seed(seed)
  if(is.null(Xtest)) {
    X <- matrix(runif(N * d), ncol = d)
  }
  else {
    X <- matrix(rep(Xtest, N), byrow = TRUE, ncol = d)
  }
  colnames(X) <- paste0("X", 1:d)
  W <- rbinom(N, 1, e(X))
  Y0 <- rnorm(N, mean = m(X) - 0.5 * tau(X), sd = 1)
  Y1 <- rnorm(N, mean = m(X) + 0.5 * tau(X), sd = 1)
  Y <- rnorm(N, mean = m(X) + (W - 0.5) * tau(X), sd = 1)

  as_tibble(X) %>%
    mutate(W = W, Y0 = Y0, Y1 = Y1, Y = Y, cate = tau(X))
}

simulate_both <- function(N, d, seed = 1, Xtest = NULL) {
  e <- function(X) 1/4 * (1 + dbeta(X[,3], 2, 4))
  m <- function(X) 2 * X[,3] - 1
  chi <- function(x) 1 + 1/(1 + exp(-20 * (x - 1/3)))
  tau <- function(X) chi(X[,1]) * chi(X[,2])
  simulate(N, d, e, m, tau, seed = seed, Xtest = Xtest)
}

simulate_effect <- function(N, d, seed = 1, Xtest = NULL) {
  e <- function(X) 0.5
  m <- function(X) 2 * X[,3] - 1
  chi <- function(x) 1 + 1/(1 + exp(-20 * (x - 1/3)))
  tau <- function(X) chi(X[,1]) * chi(X[,2])
  simulate(N, d, e, m, tau, seed = seed, Xtest = Xtest)
}

# Second simulation setup
# no treatment effect, but confounding
simulate_confounding <- function(N, d, seed = 1, Xtest = NULL) {
  e <- function(X) 1/4 * (1 + dbeta(X[,3], 2, 4))
  m <- function(X) 2 * X[,3] - 1
  tau <- function(X) 0
  simulate(N, d, e, m, tau, seed = seed, Xtest = Xtest)
}


simulate_nothing <- function(N, d, seed = 1, Xtest = NULL) {
  e <- function(X) 0.5
  m <- function(X) 2 * X[,3] - 1
  tau <- function(X) 0
  simulate(N, d, e, m, tau, seed = seed, Xtest = Xtest)
}

get_truth <- function(simulator, N, d, y, Xtest = NULL) {
  data <- simulator(N, d, seed = 23553, Xtest = Xtest)
  bandwidth_Y <- drf:::medianHeuristic(data$Y)
  k_Y <- kernlab::rbfdot(sigma = bandwidth_Y)

  K1_mean <-
    kernelMatrix(k_Y, data$Y1, y = y) %>%
    colMeans()
  K0_mean <-
    kernelMatrix(k_Y, data$Y0, y = y) %>%
    colMeans()

  approxfun(y, K1_mean - K0_mean)
}

Xtest <- matrix(c(0.7, 0.3, 0.5, 0.68, 0.43), nrow = 1)
y <- seq(-8, 8, 0.02)
nothing_truth <- function(x) 0 * x
confounding_truth <- function(x) 0 * x
effect_truth <- get_truth(simulate_effect, 8e3, 5, y = y, Xtest = Xtest)
both_truth <- get_truth(simulate_both, 8e3, 5, y = y, Xtest = Xtest)

run_simulation <- function(N, seed, d, dgp = "nothing", num_trees = 1e3, ci_group_size = 5, witness_type = "original") {
  write(glue::glue("{date()} N = {N}, seed = {seed}, dgp = {dgp}"), file = "log.txt", append = TRUE)

  if(dgp == "nothing") {
    dat <- simulate_nothing(N, d, seed = seed)
    truth <- nothing_truth
  }
  else if(dgp == "both") {
    dat <- simulate_both(N, d, seed = seed)
    truth <- both_truth
  }
  else if(dgp == "effect") {
    dat <- simulate_effect(N, d, seed = seed)
    truth <- effect_truth
  }
  else {
    dat <- simulate_confounding(N, d, seed = seed)
    truth <- confounding_truth
  }

  X <- as.matrix(select(dat, starts_with("X")))
  Y <- matrix(dat$Y, ncol = 1)
  W <- as.matrix(dat$W)
  fit <- drf_causal(X, Y, W, splitting.rule = "CausalEffectFourierMMD", num.trees = num_trees, ci.group.size = ci_group_size, num.threads = 5, response.scaling = FALSE)
  print(fit$bandwidth)

  #ghat_fit <- glm(W ~ X, family = binomial(link = "logit"))
  #ghat <- predict(ghat_fit, type = "response")
  #ghat <- ifelse(dat$W == 1, 1 / ghat, 1 / (1 - ghat))

  w0 <- get_causal_sample_weights(fit, newdata = X, newtreatment = matrix(rep(0, N), ncol = 1), g = rep(1, nrow(dat)))
  w1 <- get_causal_sample_weights(fit, newdata = X, newtreatment = matrix(rep(1, N), ncol = 1), g = rep(1, nrow(dat)))

  #w0_ghat <- get_causal_sample_weights(fit, newdata = X, newtreatment = matrix(rep(0, N), ncol = 1), g = ghat)
  #w1_ghat <- get_causal_sample_weights(fit, newdata = X, newtreatment = matrix(rep(1, N), ncol = 1), g = ghat)

  dat$Yhat0 <- (as.matrix(w0) %*% Y)[,1]
  dat$Yhat1 <- (as.matrix(w1) %*% Y)[,1]
  dat$cate_hat <- dat$Yhat1 - dat$Yhat0

  #dat$Yhat0_ghat <- (as.matrix(w0_ghat) %*% Y)[,1]
  #dat$Yhat1_ghat <- (as.matrix(w1_ghat) %*% Y)[,1]
  #dat$cate_ghat <- dat$Yhat1_ghat - dat$Yhat0_ghat

  #witness <- predict_witness(fit)
  if(witness_type == "original") {
    Xtest <- matrix(c(0.7, 0.3, 0.5, 0.68, 0.43), nrow = 1)
    witness <- predict_witness_orig(fit, alpha = 0.05, newdata = Xtest, newtreatment = matrix(1), g = rep(1, nrow(dat)))
    #witness_ghat <- predict_witness_orig(fit, alpha = 0.05, g = ghat)

    dat$witness <- witness[1,]
    dat$lower <- witness[2, ]
    dat$upper <- witness[3, ]

    #dat$witness_ghat <- witness_ghat[1,]
    #dat$lower_ghat <- witness_ghat[2, ]
    #dat$upper_ghat <- witness_ghat[3, ]
  }
  else {
    witness <- predict_witness(fit, alpha = 0.05, g = rep(1, nrow(dat)))
    #witness_ghat <- predict_witness(fit, alpha = 0.05, g = ghat)

    dat$witness <- witness[1,]
    dat$lower <- witness[1,] + qnorm(0.025) * sqrt(witness[2, ])
    dat$upper <- witness[1,] + qnorm(0.975) * sqrt(witness[2, ])

    #dat$witness_ghat <- witness_ghat[1,]
    #dat$lower_ghat <- witness_ghat[1,] + qnorm(0.025) * sqrt(witness_ghat[2,])
    #dat$upper_ghat <- witness_ghat[1,] + qnorm(0.975) * sqrt(witness_ghat[2,])
  }
  
  dat$truth <- truth(dat$Y)

  dat
}


