# This script aims at testing the framework to infer incidence of cholera
# and total AWD.

# Results from this inference will be used as priors for the infection times
# in the full model.
# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(future)

source("analysis/utils.R")

future::plan("multisession", workers = 6)

# Age cuts
age_cuts <- c(0, 5, 65, Inf)
do_plots <- T

# Load and prepare data ------------------------------------------------------------

lab_data<-readRDS("data/lab_data.rds")

tot_counts <- lab_data %>% 
  filter(date %in% makeFullDates()) %>% 
  group_by(date, age_cat) %>% 
  summarise(n_tot = n()) %>% 
  group_by(age_cat) %>% 
  complete(date = makeFullDates()) %>% 
  replace_na(list(n_tot = 0)) %>% 
  ungroup() %>% 
  arrange(date, age_cat) %>% 
  filter(date >= min(lab_data$date))

surveillance_dates <- sort(unique(tot_counts$date))

# Pass weekly numbers for positivity
period_data <- lab_data %>% 
  filter(date %in% surveillance_dates) %>% 
  mutate(week = lubridate::epiweek(date),
         year = lubridate::epiyear(date),
         epiweek = str_c(year, week, sep = "-")) %>% 
  group_by(epiweek) %>% 
  mutate(tl = min(date),
         tr = max(date)) %>% 
  group_by(epiweek, age_cat) %>% 
  summarise(n_tot = n(),
            n_pos = sum(rdt_res),
            tl = tl[1],
            tr = tr[1],
            frac_pos = n_pos/n_tot,
            lo = Hmisc::binconf(n_pos, n_tot)[2],
            hi = Hmisc::binconf(n_pos, n_tot)[3]) %>% 
  arrange(tl) %>% 
  ungroup() %>% 
  group_by(epiweek) %>% 
  # Complete missing data
  complete(age_cat = levels(lab_data$age_cat)) %>% 
  mutate(tl = ifelse(is.na(tl), 
                     min(tl, na.rm = T), tl) %>% 
           as.Date(origin = "1970-01-01"),
         tr = ifelse(is.na(tr), 
                     max(tr, na.rm = T), tr) %>% 
           as.Date(origin = "1970-01-01")) %>% 
  replace_na(list(n_tot = 0, n_pos = 0)) %>% 
  arrange(epiweek, age_cat, tl) %>% 
  ungroup()

period_data_dates <- period_data %>% 
  ungroup() %>% 
  select(epiweek, tl, tr) %>% 
  group_by(epiweek) %>% 
  slice_min(tl, n = 1, with_ties = F)

if (do_plots) {
  p_incid_ages <- period_data %>% 
    select(-hi, -lo, -frac_pos) %>% 
    pivot_longer(cols = contains("n_"),
                 values_to = "count",
                 names_to = "what") %>% 
    ggplot(aes(x = tl, y = count, color = age_cat)) +
    geom_point() +
    geom_line() +
    facet_grid(what ~ ., space = "free", scales = "free") +
    theme_bw()
  
  ggsave(p_incid_ages, filename = "figures/incidence_by_age.png",
         width = 10, height = 6)
  
  period_data %>% 
    ggplot(aes(x = tr, y = frac_pos, color = age_cat)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), alpha = .5) +
    theme_bw() +
    facet_grid(age_cat ~ .)
}


frac_rdtneg_pcrpos <- lab_data %>% 
  filter(!is.na(pcr_res), is.na(cult_res)) %>% 
  mutate(month = format(date, "%Y-%m")) %>% 
  group_by(epiweek, age_cat) %>% 
  summarise(tl = min(date),
            tmax = max(date),
            n_tot = n(),
            n_pcr_pos = sum(pcr_res),
            frac_pos = n_pcr_pos/n_tot,
            lo = Hmisc::binconf(n_pcr_pos, n_tot)[2],
            hi = Hmisc::binconf(n_pcr_pos, n_tot)[3])

# Aggregate counts of each type of testing protocol by age category
wide_count_data <- furrr::future_map_dfr(
  unique(lab_data$date),
  function(x) {
    lab_data %>% 
      filter(date == x) %>% 
      group_by(date, age_cat) %>% 
      group_modify(computeCounts)
  })

# Complete missing data
wide_count_data2 <- wide_count_data %>% 
  filter(date %in% surveillance_dates) %>% 
  ungroup() %>% 
  group_by(age_cat) %>% 
  complete(date = surveillance_dates) %>% 
  mutate(y_B = case_when(is.na(y_A) ~ list(rep(0, 2)),
                         T ~ y_B),
         y_C = case_when(is.na(y_A) ~ list(rep(0, 4)),
                         T ~ y_C),
         y_A = ifelse(is.na(y_A), 0, y_A)) %>% 
  arrange(date)

# Variables for stan
full_dates <- makeFullDates()
N_incid <- length(full_dates)
N_obs <- length(surveillance_dates)

n_awd <- periodDataToMatrix(df = tot_counts,
                            date_col = "date",
                            val_col = "n_tot")

map_obs_incid <- map_dbl(surveillance_dates, ~ which(full_dates == .))

N_pos_periods <- nrow(period_data_dates)
pos_periods_tl <- map_dbl(period_data_dates$tl, ~ which(surveillance_dates == .))
pos_periods_tr <- map_dbl(period_data_dates$tr, ~ which(surveillance_dates == .))
pos_periods_incid_tl <- map_dbl(period_data_dates$tl, ~ which(full_dates == .))
pos_periods_incid_tr <- map_dbl(period_data_dates$tr, ~ which(full_dates == .))

y_A <- periodDataToMatrix(df = wide_count_data2,
                          date_col = "date",
                          val_col = "y_A")

y_B <- periodDataToArray(df = wide_count_data2,
                         date_col = "date",
                         val_col = "y_B")

y_C <- periodDataToArray(df = wide_count_data2,
                         date_col = "date",
                         val_col = "y_C")

saveRDS(period_data, "generated_data/period_data.rds")

# Prior probability of cholera for each age class
prior_chol <- c(.2, .2, .2)

# Priors on sensitivity and specificity

# Priors on test performance from Sayeed et al. 2018: https://doi.org/10.1371/journal.pntd.0006286
# The sensitivity of Cholkit, microbiological culture, PCR and Crystal VC was
#  98% (95% CI: 88-100), 71% (95% CI: 59-81), 74% (95% CI: 59-86) and 98% (95% CI: 88-100), respectively. 
#  The specificity for V. cholerae O1 was 97% (95% CI: 89-100), 100%, 97% (95% CI: 93-99) and 98% (95% CI: 92-100), 

sens_bounds <- tribble(
  ~test, ~mean, ~lo, ~hi,
  "RDT", 0.98, 0.88, 0.99,
  "PCR", 0.74, 0.59, 0.86,
  "cul", 0.71, 0.59, 0.81
) %>%
  mutate(across(c("mean", "lo", "hi"), ~ logit(. - 1e-6))) %>%
  mutate(sd = (mean - lo)/2)

## tightening the culture specficity 
spec_bounds <- tribble(
  ~test, ~mean, ~lo, ~hi,
  "RDT", 0.97, 0.89, 1,
  "PCR", 0.97, 0.93, 0.99,
  "cul", 0.999, 0.997, 1
)  %>%
  mutate(across(c("mean", "lo", "hi"), ~ logit(. - 1e-6))) %>%
  mutate(sd = (mean - lo)/2)


ind_vec <- c(1, 1, 2, 3)
map_pos_period_sens <- map_dbl(period_data_dates$tl, ~ ifelse(. < "2021-06-29", 1, 2))
logit_sens_prior_mu <- sens_bounds$mean[ind_vec]
logit_sens_prior_sd <- sens_bounds$sd[ind_vec]

# Set prior on first period of RDT batch
logit_sens_prior_mu[1] <- 0
logit_sens_prior_sd[1] <- .75

logit_spec_prior_mu <- spec_bounds$mean[ind_vec]
logit_spec_prior_sd <- spec_bounds$sd[ind_vec]

age_categories <- unique(period_data$age_cat)

stan_data <- list(
  N_obs = N_obs,      # number of surveillance observations
  N_incid = N_incid,  # number of days for which incidence is estimated
  N_agecat = length(age_categories),
  n_awd = n_awd,      # number of suspected cases
  # n_chol = n_chol,    # number of cholera cases
  map_obs_incid = map_obs_incid,      # map from surveillance observations to incidence days
  sd_sigma_autocorr = 1,         # prior on sd of autocorrelation
  autocorr_order = 1,    # order of autocorrelation (1 or 2)
  mu_beta0_chol = 0,
  mu_beta0_nonchol = 0,
  sd_beta0_chol = 1,
  sd_beta0_nonchol = 1,
  N_pos_periods = N_pos_periods,
  pos_periods_tl = pos_periods_tl,
  pos_periods_tr = pos_periods_tr,
  pos_periods_incid_tl = pos_periods_incid_tl,
  pos_periods_incid_tr = pos_periods_incid_tr,
  y_A = y_A,
  y_B = y_B,
  y_C = y_C,
  logit_sens_prior_mu = logit_sens_prior_mu,
  logit_sens_prior_sd = logit_sens_prior_sd,
  logit_spec_prior_mu = logit_spec_prior_mu,
  logit_spec_prior_sd = logit_spec_prior_sd,
  map_pos_period_sens = map_pos_period_sens,
  prior_chol = prior_chol
)

saveRDS(stan_data, "generated_data/infer_incidence_stan_data.rds")


# Run prior ---------------------------------------------------------------


# Incidence model
incid_prior_model <- cmdstan_model("analysis/stan/infer_incidence_period_multitest_2senspec_mutliage_prior.stan")

# Sample from model
incid_prior <- incid_prior_model$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  init = .1,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1250,
  max_treedepth = 12L,
  adapt_delta = .9,
  save_warmup = F,
  refresh = 100)

incid_prior$save_object("generated_data/incid_prior_survonly_stan_output_updated_multiage_prior.rds")


# Run model ---------------------------------------------------------------

# Incidence model
incid_model <- cmdstan_model("analysis/stan/infer_incidence_period_multitest_2senspec_mutliage.stan")

# Sample from model
incid_samples <- incid_model$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  init = .1,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1250,
  max_treedepth = 12L,
  adapt_delta = .9,
  save_warmup = F,
  refresh = 100)

incid_samples$save_object("generated_data/incid_prior_survonly_stan_output_updated_multiage.rds")

# Run genquant separately
incid_genquant <- incid_model$generate_quantities(fitted_params = incid_samples,
                                                  data = stan_data,
                                                  parallel_chains = 4)

incid_genquant$save_object("generated_data/incid_prior_survonly_stan_output_updated_multiage_genquant.rds")

# Check results -----------------------------------------------------------

incid_samples <- readRDS("generated_data/incid_prior_survonly_stan_output_updated_multiage.rds")

date_limits <- c("2021-01-01", "2022-04-01") %>% as.Date()
incid_samples$summary(variables = c("sens_all", "spec_all"))

incid_samples$draws(variables = c("sens_all", "spec_all")) %>% 
  bayesplot::mcmc_trace()

age_categories <- lab_data$age_cat %>% levels()

incid_traj <- incid_samples$summary(variables = c("I_chol", "I_nonchol"), .cores = 2) %>% 
  mutate(time = str_extract(variable, "[0-9]+(?=,)") %>% as.numeric(),
         date = full_dates[time],
         age_cat_num = str_extract(variable, "(?<=,)[0-9]+") %>% as.numeric(),
         age_cat = age_categories[age_cat_num],
         var = str_extract(variable, "(.)*(?=\\[)"))

if (do_plots) {
  
  p_incid_infer <- incid_traj %>% 
    ggplot(aes(x = date, y = mean)) +
    geom_ribbon(aes(ymin = q5, ymax = q95, fill = age_cat), alpha = .2) +
    geom_line(aes(col = age_cat)) +
    facet_grid(var ~ ., scales = "free")+
    # coord_cartesian(ylim = c(0,2)) +
    theme_bw() 
  
  ggsave(p_incid_infer, filename = "figures/inferred_incidence_by_age.png",
         width = 10, height = 6)
}

probchol_traj <- incid_samples$summary(variables = c("gamma"), .cores = 2) %>% 
  mutate(time = str_extract(variable, "[0-9]+(?=,)") %>% as.numeric(),
         date = full_dates[time],
         age_cat_num = str_extract(variable, "(?<=,)[0-9]+") %>% as.numeric(),
         age_cat = age_categories[age_cat_num],
         var = str_extract(variable, "(.)*(?=\\[)"))

if (do_plots) {
  probchol_traj %>% 
    ggplot(aes(x = date, y = mean)) +
    geom_ribbon(aes(ymin = q5, ymax = q95, fill = age_cat), alpha = .2) +
    geom_line(aes(col = age_cat)) +
    facet_grid(age_cat ~ ., scales = "free")+
    theme_bw() 
}

saveRDS(incid_traj, "generated_data/incid_prior_survonly_traj_multiage.rds")
saveRDS(probchol_traj, "generated_data/est_prob_traj_survonly_multiage.rds")

# Compare generated data to observations ----------------------------------
if (do_plots) {
  wide_period_data <- furrr::future_map_dfr(
    unique(period_data$epiweek),
    function(x) {
      lab_data %>% 
        filter(epiweek == x) %>% 
        group_by(epiweek, age_cat) %>% 
        group_modify(computeCounts)
    })
  
  y_A_obs <- wide_period_data %>%
    select(epiweek, age_cat, count = y_A)
  
  y_B_obs <- wide_period_data %>%
    select(epiweek, age_cat, y_B) %>%
    mutate(counts_1 = map_dbl(y_B, ~ .[1]),
           counts_2 = map_dbl(y_B, ~ .[2])) %>%
    select(-y_B) %>%
    pivot_longer(cols = contains("counts"),
                 names_to = "case",
                 values_to = "count")
  
  
  p_rdt <- y_B_obs %>%
    pivot_wider(names_from = case,
                values_from = count) %>%
    rowwise() %>%
    mutate(tot_rdt_neg = counts_1+counts_2,
           frac_false_neg = counts_1/tot_rdt_neg,
           lo = Hmisc::binconf(counts_1, tot_rdt_neg)[2],
           hi = Hmisc::binconf(counts_1, tot_rdt_neg)[3]) %>%
    inner_join(period_data %>% select(epiweek, age_cat, tl)) %>%
    ggplot(aes(x = tl, y = frac_false_neg)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
    geom_point(aes(size = tot_rdt_neg)) +
    facet_grid(age_cat~.) +
    theme_bw() +
    labs(x = "date", y = "fraction RDT- & PCR+/total RDT-") 
  
  ggsave(p_rdt, filename = "figures/rdt_neg_frac_pos.png", width = 10, height = 6)
  
  y_C_obs <- wide_period_data %>%
    select(epiweek, age_cat, y_C) %>%
    mutate(counts_1 = map_dbl(y_C, ~ .[1]),
           counts_2 = map_dbl(y_C, ~ .[2]),
           counts_3 = map_dbl(y_C, ~ .[3]),
           counts_4 = map_dbl(y_C, ~ .[4])) %>%
    select(-y_C) %>%
    pivot_longer(cols = contains("counts"),
                 names_to = "case",
                 values_to = "count")
  
  y_A_gen <- extractCountStats(genquant = incid_samples,
                               period_data = period_data_dates,
                               variable = "y_A_gen",
                               age_categories = age_categories) %>%
    inner_join(y_A_obs) %>%
    mutate(case = "<RDT-, PCR:NA, CUL:NA>")
  
  ggplot(y_A_gen, aes(x = count, y = mean)) +
    geom_abline(aes(intercept = 0, slope = 1), col = "red") +
    geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
    geom_point() +
    facet_grid(age_cat ~ . ) +
    theme_bw()
  
  ggplot(y_A_gen, aes(x = date, y = mean)) +
    geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
    geom_point() +
    geom_point(aes(y = count), col = "red") +
    facet_grid(age_cat ~ . ) +
    theme_bw()
  
  
  y_B_gen <- extractCountStats(genquant = incid_samples,
                               period_data = period_data_dates,
                               variable = "y_B_gen",
                               age_categories = age_categories) %>%
    inner_join(y_B_obs) %>%
    mutate(case = case_when(case == "counts_1" ~ "<RDT-, PCR+, CUL:NA>",
                            T ~ "<RDT-, PCR-, CUL:NA>"))
  
  ggplot(y_B_gen, aes(x = date, y = mean)) +
    geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
    geom_point() +
    geom_point(aes(y = count), col = "red") +
    theme_bw() +
    facet_grid(age_cat ~ case)
  
  y_C_gen <- extractCountStats(genquant = incid_samples,
                               period_data = period_data_dates,
                               variable = "y_C_gen",
                               age_categories = age_categories) %>%
    inner_join(y_C_obs) %>%
    mutate(case = case_when(case == "counts_1" ~ "<RDT+, PCR+, CUL+>",
                            case == "counts_2" ~ "<RDT+, PCR+, CUL->",
                            case == "counts_3" ~ "<RDT+, PCR-, CUL+>",
                            case == "counts_4" ~ "<RDT+, PCR-, CUL->",
                            T ~ NA_character_))
  
  p_gen <- bind_rows(
    y_A_gen,
    y_B_gen,
    y_C_gen
  ) %>%
    ggplot(aes(x = date, y = mean)) +
    geom_errorbar(aes(ymin = q5, ymax = q95), width = 0, alpha = .5) +
    geom_point(alpha = .7) +
    geom_point(aes(y = count), col = "red", alpha = .7) +
    theme_bw() +
    facet_grid(case ~ age_cat, scales = "free_y")
  
  ggsave(p_gen, filename = "figures/generated_counts.png", width = 10, height = 10)
}

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
