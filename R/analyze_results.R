library(tidyverse)

#results <- read_rds("/gpfs/home/susmah01/causal_drf_paper/results/results.rds")
results <- bind_rows(
  read_rds("/gpfs/home/susmah01/causal_drf_paper/results/results_250.rds"),
  read_rds("/gpfs/home/susmah01/causal_drf_paper/results/results_500.rds"),
  read_rds("/gpfs/home/susmah01/causal_drf_paper/results/results_1000.rds"),
  read_rds("/gpfs/home/susmah01/causal_drf_paper/results/results_5000.rds")
)

analysis <- results %>%
  unnest(coverage) %>%
  mutate(covered = lower <= truth & upper >= truth, error = truth - witness) %>%
  group_by(N, dgp, d, method, num_trees, seed) %>%
  mutate(covered = all(covered)) %>%
  group_by(N, dgp, d, method, num_trees) %>%
  summarize(coverage = mean(covered), mae = mean(abs(error)), me = mean(error))

analysis_wide <- analysis %>%
  mutate(method = ifelse(method == "Combined", "comb", "sep")) %>%
  pivot_wider(names_from = c("method"), values_from = c("coverage", "mae", "me"))
