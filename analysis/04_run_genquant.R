# This is the standalone block to run the generated quantities to compute
# infection probabilities by time period of interest.
# 
# In this version probabilities are computed for the period prior to baseline (
# 6 months prior to first sample), and the period between each survey round.
# 

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(optparse)
library(here)
library(posterior)

source(here("analysis/utils.R"))

# User-supplied options
option_list <- list(
  make_option(c("-d", "--stan_data"), 
              default = here("generated_data/stan_data_vibriocidal_agegrpall_J100_ogawa.rds"), 
              action ="store", type = "character", help = "Stan input data"),
  make_option(c("-g", "--genquant"), 
              default = here("analysis/stan/seroincidence_estimation_additive_full_censored_rates_foi_envlambda_generate.stan"), 
              action ="store", type = "character", help = "Stan generated quantities"),
  make_option(c("-p", "--output_prefix"), 
              default = "final_envlambda_agegrpall_J100", 
              action ="store", type = "character", help = "Stan model")
)

opt <- parse_args(OptionParser(option_list = option_list)) 

cat("---- Running with output prefix: ", opt$output_prefix, "/n")

source("analysis/utils.R")


# Genquant block ----------------------------------------------------------

genquant_model <- cmdstan_model(opt$genquant)

# Get dates for output creation
dates_wide <- readRDS(here("data/dates_full_wide_agegrpall.rds"))
# Indices of round 1B
ind_round1B <- dates_wide$R0 > "2021-05-01"

# Create time bounds
output_periods_v1 <- tibble(
  tl = c(min(dates_wide$R0) - 180, 
         min(dates_wide$R0[ind_round1B]) - 180, 
         min(dates_wide$R0[!ind_round1B]),
         min(dates_wide$R0[ind_round1B]),
         min(dates_wide$R1),
         min(dates_wide$R0[!ind_round1B]),
         as.Date("2021-03-27"),
         as.Date("2021-03-27"),
         as.Date("2021-06-30")
  ),
  tr = c(min(dates_wide$R0), 
         min(dates_wide$R0[ind_round1B]), 
         min(dates_wide$R1), 
         min(dates_wide$R1),
         max(dates_wide$R2),
         max(dates_wide$R2),
         as.Date("2022-02-13"),
         as.Date("2021-06-29"),
         as.Date("2022-02-13")
  )
)

saveRDS(output_periods_v1, "generated_data/output_periods_prob_inf.rds")

# Create time bounds
output_periods_overall <- getOutputPeriodOverall()

# Sequential output periods for figure
output_periods <- getWeeklyOutputPeriods()
# For cumulative incidence
output_periods_cumulative <- output_periods %>% 
  mutate(tl = tl[1])

# minimum date before baseline = 6 months prior to first baseline sample
t_min <- getTmin()    # minimum date before baseline = 6 months prior to first baseline sample
t_max <- getTmax()

# Full vector of dates for incidence inference
full_dates <- seq.Date(t_min, t_max, by = "1 days")

# Make data for stan
stan_data_output <- readRDS(opt$stan_data)

stan_data_output <- makeGenOutputPeriodData(stan_data_output = stan_data_output,
                                            output_periods = output_periods,
                                            ref_dates = full_dates)

stan_data_output_cumulative <- makeGenOutputPeriodData(stan_data_output = stan_data_output,
                                                       output_periods = output_periods_cumulative,
                                                       ref_dates = full_dates)

stan_data_output_overall <- makeGenOutputPeriodData(stan_data_output = stan_data_output,
                                                    output_periods = output_periods_overall,
                                                    ref_dates = full_dates)

# Run generated quantities
genquant <- genquant_model$generate_quantities(
  fitted_params = getStanDraws(opt = opt),
  data = stan_data_output,
  parallel_chains = 4)

genquant_cumulative <- genquant_model$generate_quantities(
  fitted_params = getStanDraws(opt = opt),
  data = stan_data_output_cumulative,
  parallel_chains = 4)

genquant_overall <- genquant_model$generate_quantities(
  fitted_params = getStanDraws(opt = opt),
  data = stan_data_output_overall,
  parallel_chains = 4)

# Save output
genquant$save_object(makeGenquantOutputName(opt = opt))
genquant_cumulative$save_object(str_replace(makeGenquantOutputName(opt = opt), "\\.rds", "_cumulative.rds"))
genquant_overall$save_object(str_replace(makeGenquantOutputName(opt = opt), "\\.rds", "_overall.rds"))

# Extract statistics of interest ------------------------------------------
# 
# period_names <- c("pre_baseline", "pre_baseline", "baseline-followup1",
#                   "baseline-followup1", "followup1-followup2", "baseline-followup2", 
#                   "whole_period", "start_rdt1", "rdt2_end")
# 
# groups <- c("1A", "1B", "1A", "1B", "1A", "1A", "1A", "1A", "1A")
# 
# prob_output_stats <- genquant$summary(variables = "mu_output") %>% 
#   mutate(age_cat = str_extract(variable, "(?<=\\,)[0-9]"),
#          period_id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
#          period = period_names[period_id],
#          period = factor(period, levels = unique(period)), 
#          group = groups[period_id],
#          tl = output_periods$tl[period_id],
#          tr = output_periods$tr[period_id]) %>% 
#   {
#     x <- .
#     bind_rows(x, slice_tail(x, n = 5) %>% mutate(group = "1B"))
#   } %>%
#   rowwise() %>% 
#   mutate(tmid = mean.Date(c(tl, tr)))
# 
# 
# prob_output_stats %>% 
#   ggplot(aes(x = period, y  = mean, color = age_cat)) +
#   geom_point()
# 
# genquant_output <- list(
#   prob_output_stats = prob_output_stats,
#   output_periods = output_periods
# )
# 
# saveRDS(genquant_output, makeGenquantOutputStatsName(opt = opt))

