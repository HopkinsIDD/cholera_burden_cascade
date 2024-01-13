# This script runs inference

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(optparse)
library(here)

source(here("analysis/utils.R"))

# User-supplied options
option_list <- list(
  make_option(c("-s", "--stan_model"), 
              default = here("analysis/stan/seroincidence_estimation_additive_full_censored_rates_foi_envlambda_joint_fullpooling.stan"), 
              action ="store", type = "character", help = "Stan model"),
  make_option(c("-d", "--stan_data"), 
              default = here("generated_data/stan_data_vibriocidal_agegrpall_J100_ogawa.rds"), 
              action ="store", type = "character", help = "Stan input data"),
  make_option(c("-p", "--output_prefix"), 
              default = str_c(getTimestamp(), "_full_model_envlambda"), 
              action ="store", type = "character", help = "prefix for output filenames")
)

opt <- parse_args(OptionParser(option_list = option_list)) %>% 
  updateOptPrefix(.)

cat("---- Running with output prefix: ", opt$output_prefix, "\n")

# Update stan_data with incidence estimates -------------------------------

# Data for initialization
stan_data <- readRDS(opt$stan_data)


# Prior model -------------------------------------------------------------
# seroburden_prior_model <- cmdstan_model(str_replace(opt$stan_model, "\\.stan", "_prior.stan"))
# 
# model_prior <- seroburden_prior_model$sample(
#   data =  stan_data,
#   init = .1,
#   seed = 1234,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 250,
#   iter_sampling = 1000,
#   max_treedepth = 12L,
#   adapt_delta = .9,
#   save_warmup = T,
#   refresh = 10)
# 
# model_prior$save_object(makeStanOutputName(opt = opt) %>% str_replace("\\.rds", "_prior.rds"))

# Stan model --------------------------------------------------------------

seroburden_model <- cmdstan_model(opt$stan_model)

# Run sampler
model_fit <- seroburden_model$sample(
  data =  stan_data,
  init = map(1:4, function(x) {
    list(
      log_lambda = rnorm(stan_data$N_age_cat, -7, 1),
      log_lambda_env = rnorm(stan_data$N_age_cat, -7, 1)
    )
  }),
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 250,
  iter_sampling = 1000,
  max_treedepth = 12L,
  adapt_delta = .9,
  save_warmup = T,
  refresh = 10,
  output_dir = here("generated_data"),
  output_basename = makeStanOutputCSVBasename(opt = opt) %>% 
    str_split("data/") %>% 
    .[[1]] %>% .[2])

# Save output object
model_fit$save_object(makeStanOutputName(opt = opt))

# Save as stan csvs
model_fit$save_output_files(dir  = here("generated_data"), 
                            basename = makeStanOutputCSVBasename(opt = opt) %>% 
                              str_split("data/") %>% 
                              .[[1]] %>% .[2])

# Extract stats of interest -----------------------------------------------
# 
# # Incidence of cholera and non-cholera AWD
# incid_ts <- model_fit$summary(variables = c("I_chol", "I_nonchol"), .cores = 4) %>% 
#   mutate(time = str_extract(variable, "[0-9]+") %>% as.numeric(),
#          # date = full_dates[time],
#          var = str_extract(variable, "(.)*(?=\\[)"))
# 
# # True and adjusted probability of positive cholera test result 
# probchol_ts <- model_fit$summary(variables = c("phi", "p"), .cores = 4) %>% 
#   mutate(time = str_extract(variable, "[0-9]+") %>% as.numeric(),
#          # date = full_dates[time],
#          var = str_extract(variable, "(.)*(?=\\[)"))
# 
# beta_ts <- model_fit$summary(variables = c("log_beta_chol", "log_beta_nonchol"), .cores = 4) %>% 
#   mutate(time = str_extract(variable, "[0-9]+") %>% as.numeric(),
#          # date = full_dates[time],
#          var = str_extract(variable, "(.)*(?=\\[)"))
# 
# model_fit$summary(c("sigma_autocorr_chol", "sigma_autocorr_nonchol"))
# 
# vc_titier_data <- readRDS("generated_data/vc_ogawa_titers_full.rds")
# 
# round_periods <- vc_titier_data %>%
#   group_by(id) %>%
#   mutate(round = case_when(round == "R0" & min(date) < "2021-05-01" ~ "R0a",
#                            round == "R0" ~ "R0b",
#                            T ~ round)) %>%
#   group_by(round) %>%
#   summarise(tl = min(date),
#             tr = max(date))
# 
# dates_wide <- readRDS("generated_data/dates_full_wide.rds")
# t_min <- min(dates_wide$R0) - 180
# 
# bind_rows(incid_ts,
#           # probchol_ts,
#           beta_ts) %>%
#   ggplot(aes(x = time + t_min,  y = mean)) +
#   geom_rect(data = round_periods, inherit.aes = F,
#             aes(xmin = tl, xmax = tr, ymin = 0, ymax = Inf,
#                 fill = round), alpha = .4) +
#   geom_ribbon(aes(ymin = q5, ymax = q95, fill = var), alpha = .2) +
#   geom_line(aes(col = var)) +
#   facet_grid(var ~ ., scale = "free") +
#   theme_bw()

# Probabilities of infection for each participant
exp_periods <- c("Pre-baseline",
                 "Baseline - Round 1",
                 "Round 1 - Round 2")

prob_inf <- model_fit$summary(variables = c("mu_round"), .cores = 4) %>% 
  mutate(id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
         round =  str_extract(variable, "(?<=\\,)[0-9]+") %>% as.numeric(),
         exp_period = exp_periods[round])

# Posterior probabilities of each case
all_cases <- expand.grid(a = c(0,1),
                         b = c(0,1),
                         c = c(0,1)) %>% 
  arrange(a, b, c) %>% 
  mutate(name = str_c("<", a, ",", b, ",", c, ">"))

prob_case <- model_fit$summary(variables = c("ll"), .cores = 4) %>% 
  mutate(id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
         case = str_extract(variable, "(?<=\\,)[0-9]+") %>% as.numeric(),
         case_name = all_cases$name[case]) %>% 
  group_by(id) %>% 
  mutate(case_prob = exp(mean)/sum(exp(mean)))


model_stats <- list(
  prob_inf = prob_inf,
  prob_case = prob_case
)

saveRDS(model_stats, makeStanOutputStatsName(opt = opt))

