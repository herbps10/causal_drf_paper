library(dplyr)
library(tidyr)
library(kernlab)
library(furrr)
library(drf)
library(readr)

source("R/simulation_study_setup.R")

plan(multisession, workers = 20)

options = furrr_options(globals = c("run_simulation", "nothing_truth", "confounding_truth", "effect_truth", "both_truth", "simulate", "simulate_nothing", "simulate_confounding", "simulate_effect", "simulate_both"), packages = c("kernlab", "drf", "tidyr", "dplyr"), seed = NULL)

results <- expand_grid(
  N = c(250, 500, 1e3),
  dgp = c("nothing", "confounding", "effect", "both"),
  seed = 1:250,
  d = c(5, 10),
  #num_trees = c(2e3, 1e3 * 100)
  num_trees = c(2e3, 1e3 * 100)
) %>%
  mutate(ci_group_size = ifelse(num_trees < 5e3, 5, 1e3),
         witness_type = ifelse(num_trees < 5e3, "adaptive", "original")) %>%
  mutate(coverage = future_pmap(list(N, seed, d, dgp, num_trees, ci_group_size, witness_type), run_simulation, .options = options))

write_rds(results, "results_uniform_intervals.rds")
