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

source("https://raw.githubusercontent.com/JeffNaef/drfinference/874156b191d77b212c835cdf6f2ad9718a1ed9d2/drf-foo.R")

drf_separate <- function(data, num_trees, Xtest, ci_group_size) {
  B <- 100
  X <- as.matrix(select(data, starts_with("X")))
  Y <- matrix(data$Y, ncol = 1)
  W <- as.matrix(data$W)

  DRF0 <-
    drfCI(X = X[W == 0, , drop = FALSE],
          Y = Y[W == 0, , drop = FALSE],
          B = B, num.trees = ci_group_size)
  DRF1 <-
    drfCI(X = X[W == 1, , drop = FALSE],
          Y = Y[W == 1, , drop = FALSE],
          B = B, num.trees = ci_group_size)
  # predict the DRFs on testdata x
  DRFpred0 <- predictdrf(DRF0, x = Xtest)
  DRFpred1 <- predictdrf(DRF1, x = Xtest)

  # kernel
  bandwidth_Y <- drf:::medianHeuristic(data$Y)
  k_Y <- rbfdot(sigma = bandwidth_Y)

  DRFpred0_all <- predictdrf(DRF0, x = X)
  DRFpred1_all <- predictdrf(DRF1, x = X)

  K1 <- kernelMatrix(k_Y, Y[W == 1], y = Y[W == 1])
  K0 <- kernelMatrix(k_Y, Y[W == 0], y = Y[W == 0])
  K <- kernelMatrix(k_Y, Y[W == 0], y = Y[W == 1])

  data$Yhat0 <- (as.matrix(DRFpred0_all$weights) %*% Y[W == 0, 1])[,1]
  data$Yhat1 <- (as.matrix(DRFpred1_all$weights) %*% Y[W == 1, 1])[,1]
  data$cate_hat <- data$Yhat1 - data$Yhat0

  # simulated null distribution
  nulldist <- sapply(seq_len(length(DRFpred1$weightsb)), function(j) {
    diag(as.matrix(DRFpred0$weightsb[[j]] - DRFpred0$weights) %*% tcrossprod(K0, as.matrix(DRFpred0$weightsb[[j]] - DRFpred0$weights))) +
      diag(as.matrix(DRFpred1$weightsb[[j]] - DRFpred1$weights) %*% tcrossprod(K1, as.matrix(DRFpred1$weightsb[[j]] - DRFpred1$weights))) -
      2 * diag( as.matrix(DRFpred0$weightsb[[j]] - DRFpred0$weights) %*% tcrossprod(K, as.matrix(DRFpred1$weightsb[[j]] - DRFpred1$weights))) %>%
      as.numeric()
  })

  # Choose the right quantile
  alpha <- 0.05
  right_quantile <- quantile(nulldist, 1 - alpha)

  # witness func
  hatmun <- function(y) {
    Ky <- t(kernelMatrix(k_Y, Y, y = y))

    K1y <- t(kernelMatrix(k_Y, Y[W == 1], y = y))
    K0y <- t(kernelMatrix(k_Y, Y[W == 0], y = y))

    return(tcrossprod(K1y, as.matrix(DRFpred1$weights)) - tcrossprod(K0y, as.matrix(DRFpred0$weights)))
  }
  data$witness <- hatmun(Y) %>% as.numeric()

  # confidence bands for witness func
  data$lower <- data$witness - sqrt(right_quantile)
  data$upper <- data$witness + sqrt(right_quantile)

  data %>%
    select(Yhat0, Yhat1, cate_hat, witness, lower, upper)
}

drf_combined <- function(data, num_trees, Xtest, ci_group_size, witness_type) {
  X <- as.matrix(select(data, starts_with("X")))
  Y <- matrix(data$Y, ncol = 1)
  W <- as.matrix(data$W)
  N <- nrow(data)

  # Combined fit
  fit <- drf_causal(X, Y, W, splitting.rule = "CausalEffectFourierMMD", num.trees = num_trees, ci.group.size = ci_group_size, num.threads = 5, response.scaling = FALSE)
  print(fit$bandwidth)

  w0 <- get_causal_sample_weights(fit, newdata = X, newtreatment = matrix(rep(0, N), ncol = 1), g = rep(1, N))
  w1 <- get_causal_sample_weights(fit, newdata = X, newtreatment = matrix(rep(1, N), ncol = 1), g = rep(1, N))

  data$Yhat0 <- (as.matrix(w0) %*% Y)[,1]
  data$Yhat1 <- (as.matrix(w1) %*% Y)[,1]
  data$cate_hat <- data$Yhat1 - data$Yhat0

  if(witness_type == "original") {
    Xtest <- matrix(c(0.7, 0.3, 0.5, 0.68, 0.43), nrow = 1)
    witness <- predict_witness_orig(fit, alpha = 0.05, newdata = Xtest, newtreatment = matrix(1), g = rep(1, N))

    data$witness <- witness[1,]
    data$lower   <- witness[2, ]
    data$upper   <- witness[3, ]
  }
  else {
    witness <- predict_witness(fit, alpha = 0.05, g = rep(1, nrow(dat)))

    data$witness <- witness[1,]
    data$lower   <- witness[1,] + qnorm(0.025) * sqrt(witness[2, ])
    data$upper   <- witness[1,] + qnorm(0.975) * sqrt(witness[2, ])
  }

  data %>%
    select(Yhat0, Yhat1, cate_hat, witness, lower, upper)
}

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

  dat$index <- 1:nrow(dat)

  Xtest <- matrix(c(0.7, 0.3, 0.5, 0.68, 0.43), nrow = 1)

  # Our method
  combined <- drf_combined(dat, num_trees, Xtest, ci_group_size, witness_type)
  combined$truth <- truth(dat$Y)

  # Separate fits
  separate <- drf_separate(dat, num_trees, Xtest, ci_group_size)
  separate$truth <- truth(dat$Y)

  bind_rows(
    cbind(dat, combined) %>% mutate(method = "Combined"),
    cbind(dat, separate) %>% mutate(method = "Separate")
  )
}

