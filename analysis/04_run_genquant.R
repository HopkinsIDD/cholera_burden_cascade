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
              default = "20230510_sitakunda_joint_fullpooling_agegrpall_J100", 
              action ="store", type = "character", help = "Stan model")
)

opt <- parse_args(OptionParser(option_list = option_list)) 

cat("---- Running with output prefix: ", opt$output_prefix, "/n")

source("analysis/utils.R")


# Genquant block ----------------------------------------------------------

genquant_model <- cmdstan_model(opt$genquant)

# Get dates for output creation
dates_wide <- readRDS(here("generated_data/dates_full_wide.rds"))
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



# Scraps ------------------------------------------------------------------


# post_stratify 
# sero_dat <- readRDS("generated_data/vc_ogawa_titers_full_agegrpall.rds")
# 
# age_cat_prop <- sero_dat %>% 
#   distinct(id, age_grp) %>% 
#   count(age_grp) %>% 
#   mutate(prop = n/sum(n))
# 
# 
# postrat <- genquant %>% 
#   tidybayes::spread_draws(mu_output[i, j]) %>% 
#   rename(age_grp = j,
#          period_id = i) %>% 
#   inner_join(age_cat_prop) %>% 
#   group_by(period_id, .draw) %>% 
#   summarise(val = weighted.mean(mu_output, prop)) %>% 
#   group_by(period_id) %>% 
#   summarise(mean = mean(val),
#             q5 = quantile(val, .05),
#             q95 = quantile(val, .95))
# 
# postrat %>% 
#   mutate(tl = output_periods$tl[period_id],
#          tr = output_periods$tr[period_id]) %>% 
#   rowwise() %>% 
#   mutate(tmid = mean.Date(c(tl, tr))) %>% 
#   filter(tmid > min(dates_wide$R0)) %>% 
#   ggplot(aes(x = tmid, y  = mean)) +
#   geom_ribbon(aes(ymin = q5, ymax = q95), alpha = .3) +
#   geom_line() +
#   theme_bw()
# 

# period_data <- readRDS("generated_data/period_data.rds")
# 
# lab_data <- read_csv("data/clinical_tests_df_10282022.csv") %>% 
#   select(date = dt_interview, rdt_res = rdt_result, 
#          pcr_res = ind_pcr_result, cult_res = ind_cul_result) %>% 
#   mutate(week = lubridate::epiweek(date), #- lubridate::epiweek(date) %% 2,
#          year = lubridate::epiyear(date),
#          epiweek = str_c(year, week, sep = "-")) %>% 
#   mutate(across(contains("res"), function(x) x == "Positive"))
# 
# wide_period_data <- lab_data %>% 
#   group_by(epiweek) %>% 
#   group_modify(computeCounts)
# 
# y_A_obs <- wide_period_data %>% 
#   select(epiweek, count = y_A)
# 
# y_B_obs <- wide_period_data %>% 
#   select(epiweek, y_B) %>% 
#   mutate(counts_1 = map_dbl(y_B, ~ .[1]),
#          counts_2 = map_dbl(y_B, ~ .[2])) %>% 
#   select(-y_B) %>% 
#   pivot_longer(cols = contains("counts"),
#                names_to = "case",
#                values_to = "count")
# 
# 
# p_rdt <- y_B_obs %>% 
#   pivot_wider(names_from = case, 
#               values_from = count) %>% 
#   rowwise() %>% 
#   mutate(tot_rdt_neg = counts_1+counts_2,
#          frac_false_neg = counts_1/tot_rdt_neg,
#          lo = Hmisc::binconf(counts_1, tot_rdt_neg)[2],
#          hi = Hmisc::binconf(counts_1, tot_rdt_neg)[3]) %>% 
#   inner_join(period_data %>% select(epiweek, tl)) %>% 
#   ggplot(aes(x = tl, y = frac_false_neg)) +
#   geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
#   geom_point(aes(size = tot_rdt_neg)) +
#   theme_bw() +
#   labs(x = "date", y = "fraction RDT- & PCR+/total RDT-")
# 
# ggsave(p_rdt, filename = "figures/rdt_neg_frac_pos.png", width = 8, height = 4)
# 
# y_C_obs <- wide_period_data %>% 
#   select(epiweek, y_C) %>% 
#   mutate(counts_1 = map_dbl(y_C, ~ .[1]),
#          counts_2 = map_dbl(y_C, ~ .[2]),
#          counts_3 = map_dbl(y_C, ~ .[3]),
#          counts_4 = map_dbl(y_C, ~ .[4])) %>% 
#   select(-y_C) %>% 
#   pivot_longer(cols = contains("counts"),
#                names_to = "case",
#                values_to = "count")
# 
# y_A_gen <- extractCountStats(genquant = genquant,
#                              period_data = period_data,
#                              variable = "y_A_gen") %>% 
#   inner_join(y_A_obs) %>% 
#   mutate(case = "<RDT-, PCR:NA, CUL:NA>")
# 
# ggplot(y_A_gen, aes(x = count, y = mean)) +
#   geom_abline(aes(intercept = 0, slope = 1), col = "red") +
#   geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
#   geom_point() +
#   theme_bw()
# 
# ggplot(y_A_gen, aes(x = date, y = mean)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
#   geom_point() +
#   geom_point(aes(y = count), col = "red") +
#   theme_bw()
# 
# 
# y_B_gen <- extractCountStats(genquant = genquant,
#                              period_data = period_data,
#                              variable = "y_B_gen") %>% 
#   inner_join(y_B_obs) %>% 
#   mutate(case = case_when(case == "counts_1" ~ "<RDT-, PCR+, CUL:NA>",
#                           T ~ "<RDT-, PCR-, CUL:NA>"))
# 
# ggplot(y_B_gen, aes(x = date, y = mean)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
#   geom_point() +
#   geom_point(aes(y = count), col = "red") +
#   theme_bw() +
#   facet_grid(case ~.)
# 
# y_C_gen <- extractCountStats(genquant = genquant,
#                              period_data = period_data,
#                              variable = "y_C_gen") %>% 
#   inner_join(y_C_obs) %>% 
#   mutate(case = case_when(case == "counts_1" ~ "<RDT+, PCR+, CUL+>",
#                           case == "counts_2" ~ "<RDT+, PCR+, CUL->",
#                           case == "counts_3" ~ "<RDT+, PCR-, CUL+>",
#                           case == "counts_4" ~ "<RDT+, PCR-, CUL->",
#                           T ~ NA_character_))
# 
# 
# bind_rows(
#   y_A_gen,
#   y_B_gen,
#   y_C_gen
# ) %>% 
#   ggplot(aes(x = date, y = mean)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95), width = 0, alpha = .5) +
#   geom_point(alpha = .7) +
#   geom_point(aes(y = count), col = "red", alpha = .7) +
#   theme_bw() +
#   facet_wrap(case ~., scales = "free_y", ncol = 2)
# 
# 
# y_gen <- genquant$summary(variables = "p", .cores = 4) %>% 
#   mutate(time = str_extract(variable, "[0-9]+(?=,)") %>% as.numeric(),
#          case = str_extract(variable, "(?<=,)[0-9]+"),
#          date = period_data$tl[time],
#          var = str_extract(variable, "(.)*(?=\\[)"))
# 
# y_gen %>% 
#   ggplot(aes(x = date, y = mean)) +
#   geom_ribbon(aes(ymin = q5, ymax = q95, fill = case), alpha = .3) +
#   geom_line(aes(col = case)) +
#   theme_bw() +
#   facet_wrap(~ case)
# 
# gen_frac <- genquant$summary(variables = c("phi_avg_output"), .cores = 4,
#                              mean, ~posterior::quantile2(.x, 
#                                                          probs = c(0.025, 0.05, 0.95, 0.975))) %>%
#   bind_cols(period_data)
# 
# gen_frac %>%
#   ggplot(aes(x = tr)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95), col = "red", alpha = .5) +
#   geom_point(aes(y = mean), col = "red", alpha = .5, size = 2) +
#   geom_point(aes(y = frac_pos), alpha = .5) +
#   geom_errorbar(aes(ymin = lo, ymax = hi), alpha = .5) +
#   theme_bw() # +
# # geom_rect(data = round_periods, inherit.aes = F,
# #           aes(xmin = tl, xmax = tr, ymin = 0, ymax = Inf,
# #               fill = round), alpha = .4)
# 
# case_data_subset_full <- readRDS("generated_data/case_data_subset_full.rds")
# 
# gen_awd <- genquant$summary(variables = c("n_awd_gen"), .cores = 4) %>%
#   bind_cols(case_data_subset_full[1:nrow(.), ])
# 
# gen_awd %>%
#   ggplot(aes(x = date)) +
#   geom_ribbon(aes(ymin = q5, ymax = q95), fill = "red", alpha = .2) +
#   geom_line(aes(y = mean), col = "red", alpha = .5, size = 2) +
#   geom_line(aes(y = n_tot)) +
#   theme_bw()

