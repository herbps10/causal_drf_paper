library(dplyr)
library(tidyr)
library(kernlab)
library(furrr)
library(drf)
library(readr)
library(purrr)

source("R/simulation_study_setup.R")

plan(multisession, workers = 10)

options = furrr_options(globals = c("run_simulation", "drf_separate", "drf_combined", "nothing_truth", "confounding_truth", "effect_truth", "both_truth", "simulate", "simulate_nothing", "simulate_confounding", "simulate_effect", "simulate_both", "drfCI", "predictdrf", "drfown"), packages = c("kernlab", "drf", "tidyr", "dplyr", "Matrix"), seed = NULL)

for(N in c(5e3)) {
  N_simulations <- 500
  results <- expand_grid(
    simulation = 1:N_simulations,
    N = N,
    dgp = c("nothing", "confounding", "effect", "both"),
    #d = c(5, 10),
    d = c(5),
    #num_trees = c(2e3, 1e3 * 100)
    num_trees = c(50 * 50)
  ) %>%
    mutate(seed = 1:n()) %>%
    #mutate(ci_group_size = ifelse(num_trees < 5e3, 5, 1e3),
    #       witness_type = ifelse(num_trees < 5e3, "adaptive", "original")) %>%
    mutate(ci_group_size = round(num_trees / 50), witness_type = "original") %>%
    mutate(coverage = future_pmap(list(N, seed, d, dgp, num_trees, ci_group_size, witness_type), run_simulation, .options = options))
  
  write_rds(results, glue::glue("/gpfs/scratch/susmah01/causal_drf_paper/results_{N}.rds"))
}
