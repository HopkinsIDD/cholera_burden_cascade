# This script prepares data for the full stan model for seroincidence


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(future)
library(furrr)
library(here)
library(optparse)

source(here("analysis/utils.R"))

# User-supplied options
option_list <- list(
  make_option(c("-j", "--num_draws"), 
              default = 100, action ="store", type = "double", help = "Number of draws to use"),
  make_option(c("-a", "--age_grp"), 
              default = "all", 
              action ="store", type = "character", help = "Age group to use")
)

opt <- parse_args(OptionParser(option_list = option_list))

do_plots <- T

# Parallel setup for furrr
future::plan("multisession", workers = 6)

# Parameters in the rest of the code
# spec <- .8
# sens <- .8
J <- opt$num_draws    # This should be changed if we use multiple draws
n_subsample <- NA

set.seed(1234)

# Decay model data --------------------------------------------------------
# Load data from pervious model fits
ogawa_model <- readRDS(here::here("decay_data/model_obj_ogawa.rds"))

# Parameters:
#   omega: baseline
#   D: delay
#   lamda: boost
#   halflife = log2/delta
ogawa_draws <- readRDS(here::here("decay_data/pos_draws_ogawa.rds"))

# Parameters for stan
if (J == 1) {
  # First only use mean. Will use full distribution later
  decay_param_means <- ogawa_draws %>% 
    summarise(across(c("mu_omega", "mu_logD", "mu_lambda", "halflife"), mean))
  
  subpars <- decay_param_means%>% 
    rename(omega = mu_omega,
           lambda = mu_lambda,
           logD = mu_logD)
} else {
  
  subpars <-  posterior::bind_draws(
    rstan::extract(
      ogawa_model$fit, 
      pars = c("omega_ind", "lambda_ind", "logD_ind")) %>% 
      posterior::as_draws() %>% 
      posterior::resample_draws(ndraws = J),
    rstan::extract(
      ogawa_model$fit, 
      pars = c("halflife")) %>% 
      posterior::as_draws() %>% 
      posterior::resample_draws(ndraws = J)) %>% 
    posterior::as_draws_df() %>% 
    as_tibble() %>% 
    rename(omega = omega_ind,
           lambda = lambda_ind,
           logD = logD_ind)
}

# Unpack
gamma <- subpars %>% pull(omega) %>% exp()
beta <- subpars %>% pull(lambda) %>% exp()
alpha <- log(2)/(subpars %>% pull(halflife))
delta <- subpars %>% pull(logD) %>% exp()

# Load vibriocidal titer data and in wide data format to load to stan -----------------------------------------------------------

vc_titier_data <- readRDS(str_glue("data/vc_ogawa_titers_full_agegrpall.rds"))
dates_wide <- readRDS(str_glue("data/dates_full_wide_agegrpall.rds"))
titers_wide <- readRDS(str_glue("data/titers_full_wide_agegrpall.rds"))

# Dilutions levels
dilutions <- vc_titier_data %>% 
  pull(titer) %>% 
  unique() %>% 
  sort() 

# Model with censoring
y_titer <- titers_wide %>% 
  select(-id) %>% 
  rowwise() %>% 
  summarise(across(.fns = function(x) which(dilutions == x))) %>% 
  as.matrix()

y <- titers_wide %>% 
  select(-id) %>% 
  as.matrix() %>% 
  {log(.)}

dilutions <- dilutions %>%
  # Add 9999 as Inf upper bound
  c(., 9999) %>% 
  log()

# Define numbers for stan
N <- nrow(titers_wide)    # Number of participants


# Define times ------------------------------------------------------------
# Set times in model. Note that all times in the model are relative to the 
# sampling dates and are specific to each participant.
# Definitions: 
#   - t_follow_1: dt from baseline sample to round 1 sample
#   - t_follow_2: dt from baseline sample to round 2 sample
#   - mdt_pre0: dt from midpoint of interval for boost pre-baseline to baseline sample
#   - mdt_0_1:  dt from midpoint of interval for boost between baseline and round 1 to round 1 sample
#   - mdt_1_2:  dt from midpoint of interval for boost between round 1 and round 2 to round 2 sample

# --- A. Times to followups ---
t_follow_1 <- difftime(dates_wide$R1, dates_wide$R0, units = "days") %>% as.numeric()
t_follow_2 <- difftime(dates_wide$R2, dates_wide$R0, units = "days") %>% as.numeric()

# --- B. Times to boosts ---
#  Note that each period (pre-baseline, baseline-round 1, round 1-round 2)
# can be subdivided in an arbitrary number of intervals.
t_min <- getTmin()   # minimum date before baseline = 6 months prior to first baseline sample
t_max <- getTmax()

# Consider at first that each period is subdivided into the same number of intervals
M <- c(12, 12, 12)    # Number of intervals into which to devide each period.

# Define intervals
mdt_pre0 <- matrix(0, nrow = N, ncol = M[1])
mdt_0_1 <- matrix(0, nrow = N, ncol = M[2])
mdt_1_2 <- matrix(0, nrow = N, ncol = M[3])

for (i in 1:N) {
  mdt_pre0[i, ] <- computeSeqMidpointsDT(t_left = t_min,
                                         t_right = dates_wide$R0[i],
                                         n_intervals = M[1])
  
  mdt_0_1[i, ] <- computeSeqMidpointsDT(t_left = dates_wide$R0[i],
                                        t_right = dates_wide$R1[i],
                                        n_intervals = M[2])
  
  mdt_1_2[i, ] <- computeSeqMidpointsDT(t_left = dates_wide$R1[i],
                                        t_right = dates_wide$R2[i],
                                        n_intervals = M[3])
}


# Define probabilities of infection ---------------------------------------
# Probabilities of infection for each participant within each boost 
# interval. We assume that the probability of infection is proportional
# to the number of reported cases.


# --- A. compute smooth case date for probabilities ---
## Surveillance data
#case_data <- read_csv("data/clinical_sitakunda.csv") %>% 
#  rename(date = start) %>% 
#  filter(sitakunda == "Yes") %>% 
#  group_by(date) %>% 
#  summarise(n_tot = n(),
#            n_pos = sum(rdt_result == "Positive")) %>% 
#  ungroup()

## Use incubation period + reporting delay to shift reported cases to infection
#inc_period <- 2   
#report_delay <- 2

## Select surveillance data within the time bounds of foi
#case_data_subset <- case_data %>% 
#  complete(date = seq.Date(min(case_data$date), 
#                           max(case_data$date),
#                           by = "1 days")) %>% 
#  replace_na(list(n_tot = 0, n_pos = 0)) %>% 
#  filter(date >= t_min, date <= t_max) %>% 
#  arrange(date)

## Full vector of dates for incidence inference
#full_dates <- makeFullDates()

## Pad data for period prior to start of observations
#case_data_subset_full <- case_data_subset %>%  
#  complete(date = seq.Date(t_min, 
#                           max(case_data$date),
#                           by = "1 days")) %>% 
#  mutate(yday = lubridate::yday(date),
#         year = lubridate::year(date)) %>%
#  add_count(yday) %>%
#  {
#    x <- . 
#    x <- x %>% 
#      group_by(yday) %>% 
#      mutate(n_na = sum(is.na(n_tot)))
    
#    bind_rows(
#      x %>% filter(n == 1 | n_na == 0, date >= min(case_data$date)),
#      x %>% filter(n > 1, n_na > 0) %>% 
#        mutate(n_tot = n_tot[!is.na(n_tot)],
#               n_pos = n_pos[!is.na(n_pos)])
#    )
#  } %>% 
#  ungroup() %>% 
#  complete(date = full_dates) %>% 
#  # Fill missing dates
#  replace_na(list(n_tot = 1, n_pos = 0)) %>% 
#  arrange(date) %>% 
#  filter(date >= t_min, date <= t_max)


## Smooth timeseries
#case_data_smooth <- case_data_subset_full %>%
#  # Compute weekly seropositivity
#  mutate(week = str_c(lubridate::year(date),
#                      lubridate::epiweek(date), sep = "-")) %>% 
#  group_by(week) %>% 
#  mutate(chol_pos = (sum(n_pos)+1)/sum(n_tot),
#         # Invert to true probability based on sens and spec
#         true_pos = pmax(.01, (chol_pos - 1 + .9)/(.9 - 1 + .9))) %>% 
#  ungroup() %>% 
#  mutate(n_chol_est = (n_tot * true_pos)) %>% 
#  mutate(cases_shifted = dplyr::lead(n_chol_est, inc_period + report_delay),
#         cases_smooth = forecast::ma(cases_shifted, 28)) %>% 
#  rowwise() %>% 
#  mutate(cases_shifted = case_when(is.na(cases_shifted) ~ n_chol_est,
#                                   T ~ cases_shifted),
#         cases_smooth = case_when(is.na(cases_smooth) ~ as.numeric(cases_shifted),
#                                  T ~ cases_smooth)) %>% 
#  ungroup()

#if (do_plots) {
#  case_data_smooth %>% 
#    pivot_longer(cols = c("n_tot", "n_pos", "n_chol_est", "cases_shifted", "cases_smooth")) %>% 
#    bind_rows(vc_titier_data %>% 
#                select(date, name = round) %>% 
#                group_by(date, name) %>% 
#                summarise(value = n())) %>% 
#    ggplot(aes(x = date, y = value, fill = name)) + 
#    geom_bar(stat = "identity") +
#    facet_grid(name ~ ., scales = "free_y") +
#    theme_bw()
#  
#  case_data_smooth %>% 
#    mutate(year = lubridate::year(date) %>% factor(),
#           yday = lubridate::yday(date)) %>% 
#    ggplot(aes(x = yday, y = n_chol_est, color = year)) + 
#    geom_line() +
#    theme_bw()
#}

# Read in case_data_smooth (smooth case date for probabilities) previously generated 
case_data_smooth <- readRDS("data/surveillance_data_padded.rds")
case_data_subset_full <- readRDS("data/case_data_subset_full.rds")

# --- B. Compute probabilities  ---
# Using furrr for parallel computations

# Probabilities before baseline
rprob_inf_pre0 <- future_map(
  1:N,
  function(x) {
    computeIntervalProbs(
      intervals = computeIntervals(t_left = t_min,
                                   t_right = dates_wide$R0[x],
                                   n_intervals = M[1]),
      cases = case_data_smooth
    )
  }) %>% 
  Reduce(f = "rbind") %>%
  as.matrix()

if (do_plots) {
  image(rprob_inf_pre0)
}

# Probabilities between baseline and round 1
rprob_inf_0_1 <-  future_map(
  1:N,
  function(x) {
    computeIntervalProbs(
      intervals = computeIntervals(t_left = dates_wide$R0[x],
                                   t_right = dates_wide$R1[x],
                                   n_intervals = M[2]),
      cases = case_data_smooth
    )
  }) %>% 
  Reduce(f = "rbind") %>%
  as.matrix()

if (do_plots) {
  image(rprob_inf_0_1)
}

# Probabilities between round 1 and round 2
rprob_inf_1_2 <- future_map(
  1:N,
  function(x) {
    computeIntervalProbs(
      intervals = computeIntervals(t_left = dates_wide$R1[x],
                                   t_right = dates_wide$R2[x],
                                   n_intervals = M[3]),
      cases = case_data_smooth
    )
  }) %>% 
  Reduce(f = "rbind") %>%
  as.matrix()

if (do_plots) {
  image(rprob_inf_1_2)
}


# Make plot to verify periods
t_m_pre0 <- tibble()
t_m_0_1 <- tibble()
t_m_1_2 <- tibble()

if (do_plots) {
  for (i in 1:M[1]) {
    t_m_pre0 <- bind_rows(t_m_pre0, 
                          tibble(date = dates_wide$R0 - mdt_pre0[,i],
                                 m = i,
                                 id = dates_wide$id,
                                 probs = rprob_inf_pre0[, i]))
  }
  for (i in 1:M[2]) {
    t_m_0_1 <- bind_rows(t_m_0_1, tibble(date = dates_wide$R1 - mdt_0_1[,i],
                                         m = i,
                                         id = dates_wide$id,
                                         probs = rprob_inf_0_1[, i]))
  }
  for (i in 1:M[3]) {
    t_m_1_2 <- bind_rows(t_m_1_2, tibble(date = dates_wide$R2 - mdt_1_2[,i],
                                         m = i,
                                         id = dates_wide$id,
                                         probs = rprob_inf_1_2[, i]))
  }
  
  t_m_1_2 %>%
    filter(id == dates_wide$id[1]) %>% 
    ggplot(aes(x = date, y = 1.5)) + 
    geom_point(data = t_m_pre0 %>%
                 filter(id == dates_wide$id[1]),
               aes(color = probs, size = probs, y = .5)) +
    geom_point(data = t_m_0_1 %>%
                 filter(id == dates_wide$id[1]),
               aes(color = probs, size = probs, y = 1)) +
    geom_point(aes(color = probs, size = probs)) +
    theme_bw() +
    geom_line(data = case_data_smooth, 
              aes(x = date, y = cases_smooth)) +
    scale_color_viridis_b()
}

# Cumulative force of infection -------------------------------------------
# Define cumulative force of infection for prob of infection
# 
# cum_foi <- map(
#   1:nrow(dates_wide), 
#   function(i) {
#     t0 <- dates_wide$R0[i]
#     t1 <- dates_wide$R1[i]
#     t2 <- dates_wide$R2[i]
#     # Sum smoothed cases
#     with(
#       case_data_padded,
#       c(
#         sum(cases_smooth[date > t_min & date <= t0]),
#         sum(cases_smooth[date > t0 & date <= t1]),
#         sum(cases_smooth[date > t1 & date <= t2])
#       ))
#   }) %>% 
#   Reduce(f = "rbind") %>% 
#   as.matrix()

# In new setup define time-left and time-right of force of infection (foi)
foi_tl <- dates_wide %>% 
  filter(id %in% titers_wide$id) %>% 
  mutate(pre0 = R0 - 180) %>% 
  select(pre0, R0, R1) %>%
  # Set period pre-baseline to 6months prior to baseline
  mutate(across(.cols = c("pre0", "R0", "R1"), 
                function(x) {
                  difftime(x, t_min, units = "days") %>% 
                    as.numeric() %>% 
                    {. + 1}
                })) %>% 
  as.matrix()

foi_tr <- dates_wide %>% 
  filter(id %in% titers_wide$id) %>% 
  select(R0, R1, R2) %>%
  # Set period pre-baseline to 6months prior to baseline
  mutate(across(.cols = c("R0", "R1", "R2"), 
                function(x) {
                  difftime(x, t_min, units = "days") %>% 
                    as.numeric()
                })) %>% 
  as.matrix()

if (do_plots) {
  ggplot(case_data_subset_full,
         aes(x = date, y = n_tot)) +
    geom_line() +
    geom_line(data = case_data, col = "red", alpha = .5) +
    theme_bw()
}


# For incidence observations only keep observed data
case_data_subset_filtered <- filter(case_data_subset_full,
                                    date >= min(case_data$date),
                                    date <= max(case_data$date))


# Save section objects
saveRDS(case_data_subset_full, "generated_data/case_data_subset_full.rds")
saveRDS(case_data_subset_filtered, "generated_data/case_data_subset_filtered.rds")


# FOI data from incidence inference ---------------------------------------

# Surveillance-only incidence fit
incid_fit <- readRDS("generated_data/incid_prior_survonly_stan_output_updated_multiage.rds")

if (opt$age_grp != "all") {
  # Index of age category to pull
  age_cat_num <- as.numeric(opt$age_grp)
  
  chol_incid_draws <- incid_fit$draws("I_chol") %>% 
    # Keep only draws for age group of interest
    posterior::as_draws() %>%
    posterior::subset_draws(variable = str_glue("I_chol\\[[0-9]+,{age_cat_num}\\]"), 
                            regex = T,
                            iteration = 1:50) %>% 
    posterior::as_draws_matrix() %>% 
    t()
  
  I_chol <- matrix(chol_incid_draws, 
                   ncol = ncol(chol_incid_draws)) %>% as.array()
  J_foi <- ncol(I_chol)
} else {
  
  
  chol_incid_draws <- map(
    1:3, 
    function(x) {
      incid_fit$draws("I_chol") %>% 
        # Keep only draws for age group of interest
        posterior::as_draws() %>%
        posterior::subset_draws(variable = str_glue("I_chol\\[[0-9]+,{x}\\]"), 
                                regex = T,
                                iteration = 1:50) %>% 
        posterior::as_draws_matrix() %>% 
        t()
    })
  
  I_chol <- map(chol_incid_draws, 
                ~ matrix(., 
                         ncol = ncol(.))) %>% as.array()
  
  J_foi <- ncol(I_chol[[1]]) %>% as.integer()
}

# Vector for age categories
age_cat_vec <- titers_wide %>% 
  inner_join(enroll_dates %>% distinct(id, age_grp)) %>% 
  pull(age_grp)

incidence_stan_data <- readRDS("generated_data/infer_incidence_stan_data.rds")

# Stan data ---------------------------------------------------------------

stan_data <- list(
  N = N,    # number sero tested
  J = J,    # number of draws of decay params
  M = M,    # number of time slices
  K = length(dilutions),    # number of dilutions in series (only used for censored model)
  N_age_cat = ifelse(opt$age_grp == "all", 3, 1) %>% as.integer(),
  y_titer = y_titer,        # titer dilution values
  y = y,                    # titer values
  dilutions = dilutions,    # dilutions series (only used for censored model)
  age_cat_vec = age_cat_vec,
  mdt_pre0 = mdt_pre0,      # time between possible infection dates prior to baseline and baseline survey
  mdt_0_1 = mdt_0_1,        # time between possible infection dates between baseline survey and round 1
  mdt_1_2 = mdt_1_2,        # time between possible infection dates between round 1 and round 2
  rprob_inf_pre0 = rprob_inf_pre0,
  rprob_inf_0_1 = rprob_inf_0_1,
  rprob_inf_1_2 = rprob_inf_1_2,
  gamma = as.array(gamma),    # baseline titer
  beta = as.array(beta),      # boost
  alpha = as.array(alpha),    # decay
  delta = as.array(delta),    # delay
  t_follow = cbind(t_follow_1, t_follow_2), # follow-up time
  sigma = 0.05,       # antibody titer measurement error
  N_incid = incidence_stan_data$N_incid,
  foi_tl = foi_tl,    # left time bound for foi computation
  foi_tr = foi_tr,    # right time bound for foi computation
  J_foi = J_foi,
  I_chol = I_chol,
  sd_lambda = 1
  # N_obs = N_obs,      # number of surveillance observations
  # N_incid = N_incid,  # number of days for which incidence is estimated
  # n_awd = n_awd,      # number of suspected cases
  # n_chol = n_chol,    # number of cholera cases
  # map_obs_incid = map_obs_incid,      # map from surveillance observations to incidence days
  # sd_sigma_autocorr = .075,         # prior on sd of autocorrelation
  # autocorr_order = 1,    # order of autocorrelation (1 or 2)
  # mu_beta0_chol = 0,
  # mu_beta0_nonchol = 0,
  # sd_beta0_chol = 1,
  # sd_beta0_nonchol = 1#,
  # N_pos_periods = N_pos_periods,
  # pos_periods_tl = pos_periods_tl,
  # pos_periods_tr = pos_periods_tr,
  # pos_periods_incid_tl = pos_periods_incid_tl,
  # pos_periods_incid_tr = pos_periods_incid_tr,
  # y_A = y_A,
  # y_B = y_B,
  # y_C = y_C,
  # logit_sens_prior_mu = logit_sens_prior_mu,
  # logit_sens_prior_sd = logit_sens_prior_sd,
  # logit_spec_prior_mu = logit_spec_prior_mu,
  # logit_spec_prior_sd = logit_spec_prior_sd,
  # map_pos_period_sens = map_pos_period_sens
)

# cmdstanr::write_stan_json(stan_data, file = "generated_data/stan_data_vibriocidal_ogawa.json")
saveRDS(stan_data, makeStanDataFile(opt))

# Scraps ------------------------------------------------------------------


# plot trajectories to check
# upped J for this
# trajs_mat <- lapply(1:J,function(i){
#   rbind(
#     data.frame(
#     t=t,
#     titer=gamma[i] + beta[i] * exp(-alpha[i] * t),
#     type="decay",
#     iter=i),
#     data.frame(
#       t=t,
#       titer=rep(gamma[i],length(t)),
#       type="baseline",
#     iter=i))
# }) %>% do.call(rbind,.)
# 
# envs <- trajs_mat %>% group_by(t,type) %>% summarize(ql=quantile(titer,.1),
#                                                      q50=quantile(titer,.5),
#                                                      qu=quantile(titer,.9))
# 
# envs %>% 
#   ggplot(aes(x=t,y=q50,fill = type,color=type,fill =type)) +
#   geom_line() + 
#   geom_ribbon(aes(ymin=ql,
#                   ymax=qu,group=type),alpha=.5) + 
#   scale_y_continuous(trans = "log2") + 
#   scale_fill_brewer(palette ="Dark2") + 
#   scale_color_brewer(palette ="Dark2") + 
#   ylab("vibriocidal titer") + 
#   xlab("time (days)") + theme_minimal() +  theme(legend.position="bottom")

