library(tidyverse)
library(magrittr)
library(kableExtra)
library(gtsummary)
library(here)
library(posterior)
library(bayesplot)
library(tidybayes)
library(cmdstanr)
source(here("analysis/utils.R"))


# Colors
this_blue <- "#493D99"
this_red <- "#BD3E5C"



# Sketch of modeling framework --------------------------------------------

# Clinical data
lab_data <- readRDS(here("data/lab_data.rds"))
case_data_subset_filtered <- readRDS(here("generated_data/case_data_subset_filtered.rds"))

theme_sketch <- theme_classic() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 14))

# Plot for totals
lab_data <- lab_data %>% 
  filter(date %in% case_data_subset_filtered$date) %>% 
  addEpiWeek() %>% 
  mutate(obs_id = row_number()) %>% 
  group_by(epiweek) %>% 
  arrange(epiweek, desc(rdt_res), desc(pcr_res), desc(cult_res)) %>% 
  mutate(obs_id_week = row_number(),
         tl = min(date),
         tr = max(date)) %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr))) %>% 
  ungroup()

p_overview_surv <- lab_data %>% 
  mutate(rdt_pcr = str_c("RDT", ifelse(rdt_res, "+", "-"), " & PCR", ifelse(pcr_res, "+", "-")),
         rdt_pcr = ifelse(is.na(rdt_pcr), "RDT- & no PCR", rdt_pcr)) %>% 
  ggplot(aes(x = tmid)) +
  geom_bar(aes(fill = rdt_pcr,
               x = tmid), 
           inherit.aes = F,
           width = 6) +
  scale_fill_manual(values = rev(c("#800ABF", "#CF8FF7", "#857E7E", "#DE6B0D", "#F7973E"))) +
  theme_bw() +
  labs(y = "count", x = "date")  +
  theme_sketch +
  theme(legend.position = "top",
        legend.key.size = unit(.1, "in")) + 
  guides(fill = guide_legend(NULL, nrow = 2))


p_overview_surv
ggsave(p_overview_surv, filename = "figures/overview_sketch_surv.pdf", width = 4, height = 3)


# Titer data
vc_titier_data <- readRDS(here("data/vc_ogawa_titers_full_agegrpall.rds"))

# a. Titer trajectories for a few participants
set.seed(1123)
sketch_ids <- sample(unique(vc_titier_data$id), 5)
p_overview_titers <- vc_titier_data %>% 
  filter(id %in% sketch_ids) %>% 
  ggplot(aes(x = date, y = log(titer/5), group = id)) +
  geom_point(size = 2) +
  geom_line(alpha = 1, linewidth = 1) +
  labs(x = "date", y = "titer") +
  theme_sketch

ggsave(p_overview_titers, filename = "figures/overview_sketch_titer.pdf", width = 4, height = 3)

# b. Inferred incidence
p_overview_incid <-  readRDS(here("generated_data/weekly_I_sum.rds")) %>% 
  filter(tmid > getTmin() + 180) %>% 
  ggplot(aes(x = tmid, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = 'gray') +
  geom_line(linewidth = 1) +
  labs(x = "date", y = "clinical cholera\nsymptomatic incidence") +
  coord_cartesian(ylim = c(0, 50)) +
  theme_sketch

ggsave(p_overview_incid, filename = "figures/overview_sketch_incid.pdf", width = 4, height = 3)

# c. Inferred incidence in community
p_overview_incid_true <-  readRDS(here("generated_data/weekly_I_sum_true.rds")) %>% 
  filter(tmid > getTmin() + 180) %>% 
  ggplot(aes(x = tmid, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = 'gray') +
  geom_line(linewidth = 1) +
  labs(x = "date", y = "community cholera \nsymptomatic incidence") +
  coord_cartesian(ylim = c(0, 50)) +
  theme_sketch

ggsave(p_overview_incid_true, filename = "figures/overview_sketch_incid_true.pdf", width = 4, height = 3)

# d. Inferred probability of infection
p_overview_prob <- readRDS(here("generated_data/prob_output_cumulative_overall_stats.rds")) %>% 
  filter(tmid > getTmin() + 180) %>% 
  ggplot(aes(x = tmid, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = 'gray') +
  geom_line(linewidth = 1) +
  labs(x = "date", y = "cumulative\nprobability of infection") +
  theme_sketch

ggsave(p_overview_prob, filename = "figures/overview_sketch_prob.pdf", width = 4, height = 3)

p_overview_col <- cowplot::plot_grid(
  p_overview_incid,
  p_overview_incid_true,
  p_overview_prob,
  ncol = 1,
  align = "v",
  axis = "lr"
)

ggsave(p_overview_col, filename = "figures/overview_sketch_col.pdf", width = 4, height = 10)


# Posterior of decay model Jones et al. 2022 ------------------------------

# Parameters:
#   omega: baseline
#   D: delay
#   lamda: boost
#   halflife = log2/delta
ogawa_draws <- readRDS(here::here("decay_data/pos_draws_ogawa.rds"))

kin_dict <- c(
  "omega" = "Basline titer [U/mL]",
  "delta" = "Delay in response [days]",
  "lambda" = "Antibody boost [U/mL]",
  "halflife" = "Antibody half-life [days]"
)

p_kin <- ogawa_draws %>% 
  mutate(omega = exp(mu_omega),     # baseline titer
         delta = exp(mu_logD),      # delay
         lambda = exp(mu_lambda),   # boost
         halflife                   # halflife of response in days
  ) %>% 
  select(-contains("_")) %>% 
  pivot_longer(cols = -contains(".")) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(fill = this_blue) +
  facet_wrap(~ name, labeller = labeller(name = kin_dict),
             scales = "free") +
  theme_bw() 

ggsave(p_kin, filename = "figures/SF_decay_model_params.png", width = 7, height = 4, dpi = 300)


# PP checks incidence model -----------------------------------------------
# Load draws
genquant_incid <- readRDS(here("generated_data/incid_prior_survonly_stan_output_updated_multiage_genquant.rds"))

lab_data <- readRDS("generated_data/lab_data.rds")
age_categories <- lab_data$age_cat %>% levels()

# Dictionnary for age categories
age_cat_dict <- c("< 5", "5-64", "65+", "overall")
names(age_cat_dict) <- age_categories

period_data <- readRDS("generated_data/period_data.rds")
period_data_dates <- period_data %>% 
  ungroup() %>% 
  select(epiweek, tl, tr) %>% 
  group_by(epiweek) %>% 
  slice_min(tl, n = 1, with_ties = F)

# Extract data
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

y_gen <- extractCountStats(genquant = genquant_incid,
                           period_data = period_data_dates,
                           variable = "y_gen",
                           age_categories = age_categories) 

# Extract generated quantities
y_A_gen <- extractCountStats(genquant = genquant_incid,
                             period_data = period_data_dates,
                             variable = "y_A_gen",
                             age_categories = age_categories) %>%
  inner_join(y_A_obs) %>%
  mutate(case = "<RDT-, PCR:NA, CUL:NA>")


y_B_gen <- extractCountStats(genquant = genquant_incid,
                             period_data = period_data_dates,
                             variable = "y_B_gen",
                             age_categories = age_categories) %>%
  inner_join(y_B_obs) %>%
  mutate(case = case_when(case == "counts_1" ~ "<RDT-, PCR+, CUL:NA>",
                          T ~ "<RDT-, PCR-, CUL:NA>"))

y_C_gen <- extractCountStats(genquant = genquant_incid,
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
  mutate(age_cat = age_cat_dict[age_cat],
         age_cat = factor(age_cat, levels = age_cat_dict)) %>% 
  ggplot(aes(x = date, y = mean)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0, alpha = .5) +
  geom_point(alpha = .7) +
  geom_point(aes(y = count), col = this_red, alpha = .7) +
  theme_bw() +
  facet_grid(case ~ age_cat, scales = "free_y") +
  theme(strip.text = element_text(size = 5))

ggsave(p_gen, filename = "figures/SF_pp_checks_incidence.png", width = 10, height = 8, dpi = 300)


p_gen_scatter <- bind_rows(
  y_A_gen,
  y_B_gen,
  y_C_gen
) %>%
  ggplot(aes(x = count, y = mean)) +
  geom_abline(linetype = 2, linewidth = .8, col = this_red) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0, alpha = .3) +
  geom_point(alpha = .4) +
  theme_bw() +
  facet_grid(case ~ age_cat, scales = "free_y") +
  theme(strip.text = element_text(size = 5))

ggsave(p_gen_scatter, filename = "figures/SF_pp_checks_incidence_scatter.png", width = 10, height = 8, dpi = 300)


# Params incid ------------------------------------------------------------

res_incid <- readRDS(here("generated_data/incid_prior_survonly_stan_output_updated_multiage.rds"))
prior_incid <- readRDS(here("generated_data/incid_prior_survonly_stan_output_updated_multiage_prior.rds"))

# Parameters of interest to plot
incid_params <- c("alpha_chol", "alpha_nonchol", "sens_all", "spec_all")

incid_params_dict <- c(
  "alpha_chol[1]" = "intercept cholera (log)",
  "alpha_chol[2]" = "intercept cholera (log)",
  "alpha_chol[3]" = "intercept cholera (log)",
  "alpha_nonchol[1]" = "intercept non-cholera AWD (log)",
  "alpha_nonchol[2]" = "intercept non-cholera AWD (log)",
  "alpha_nonchol[3]" = "intercept non-cholera AWD (log)",
  "sens_all[1]" = "sensitivity RDT\n[< June 29 2021]",
  "sens_all[2]" = "sensitivity RDT\n[> June 29 2021]",
  "sens_all[3]" = "sensitivity RT-PCR",
  "sens_all[4]" = "sensitivity culture",
  "spec_all[1]" = "specificity RDT\n[< June 29 2021]",
  "spec_all[2]" = "specificity RDT\n[> June 29 2021]",
  "spec_all[3]" = "specificity RT-PCR",
  "spec_all[4]" = "specificity culture"
)

incid_param_draws <- bind_rows(
  res_incid$draws(incid_params) %>% 
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble() %>% 
    mutate(set = "posterior"),
  prior_incid$draws(incid_params) %>% 
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble() %>% 
    mutate(set = "prior")
) %>% 
  pivot_longer(cols = contains(incid_params)) %>% 
  mutate(param = incid_params_dict[name],
         age_cat = case_when(str_detect(name, "alpha") ~ age_cat_dict[as.numeric(str_extract(name, "[1-3]"))],
                             T ~ NA_character_),
         age_cat = factor(age_cat, levels = age_cat_dict))

p_incid_intercept <- incid_param_draws %>% 
  filter(str_detect(param, "intercept")) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = set), position = "dodge") +
  facet_grid(age_cat ~ param) +
  scale_fill_manual(values = c(this_red, this_blue)) +
  theme_bw()

p_incid_tests <- incid_param_draws %>% 
  filter(str_detect(param, "intercept", negate = T)) %>% 
  mutate(perf = str_extract(param, "sensitivity|specificity"),
         param = str_remove(param, perf) %>% str_trim(),
         param = factor(param, levels = c("RDT\n[< June 29 2021]", 
                                          "RDT\n[> June 29 2021]",
                                          "RT-PCR", 
                                          "culture"))) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = set), position = "dodge") +
  facet_grid(param ~ perf, scales = "free") +
  scale_fill_manual(values = c(this_red, this_blue)) +
  theme_bw()

p_params_incid <- cowplot::plot_grid(
  p_incid_intercept + theme(legend.position = "bottom"),
  p_incid_tests + theme(legend.position = "bottom"),
  nrow = 1,
  labels = "auto"
)

ggsave(p_params_incid, filename = "figures/SF_params_incid.png",
       width = 10, height = 6, dpi = 300)

# PP checks serology ------------------------------------------------------
# Titer data
vc_titier_data <- readRDS(here("data/vc_ogawa_titers_full_agegrpall.rds"))
# Exposure classes based on serology
prob_case <-  readRDS(here("generated_data/stan_output_stats_final_envlambda_agegrpall_J100_ogawa.rds"))$prob_case

# From the rethinking package
log_sum_exp <- function(x) {
  xmax <- max(x)
  xsum <- sum(exp(x - xmax))
  xmax + log(xsum)
}

p_traj <- prob_case %>% 
  # recompute probabilities
  group_by(id) %>% 
  mutate(case_prob = exp(mean - log_sum_exp(mean))) %>% 
  summarise(case_max_prob = case_name[which.max(case_prob)],
            max_prob = max(case_prob, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(id = unique(vc_titier_data$id)[id]) %>% 
  inner_join(vc_titier_data) %>% 
  mutate(age_cat = getAgeCatDict()[age_grp],
         age_cat = factor(age_cat, levels = getAgeCatDict())) %>% 
  ggplot(aes(x = date, y = log2(titer/5), group = id, color = max_prob)) +
  geom_point(alpha = .5, size = .6) +
  geom_line(alpha = .2, linewidth = .5) +
  theme_bw() +
  facet_grid(age_cat ~  case_max_prob) +
  scale_color_viridis_b(direction = -1) +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.title = element_text(size = 10)) +
  guides(color = guide_colorbar("Maximum\nposterior\nprobability"))

ggsave(p_traj, filename = "figures/SF_pp_checks_serology.png",
       width = 10, height = 5, dpi = 300)

# Params serology ---------------------------------------------------------

res_sero <- readRDS(here("generated_data/stan_output_20230510_sitakunda_joint_fullpooling_agegrpall_J100_ogawa.rds"))
prior_sero <- readRDS(here("generated_data/stan_output_20230510_sitakunda_joint_fullpooling_agegrpall_J100_ogawa_prior.rds"))

sero_params <- str_c("lambda[", 1:3, "]") %>% 
  c(str_c("lambda_env[", 1:3, "]"))

sero_param_draws <- bind_rows(
  res_sero$draws(sero_params) %>% 
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble() %>% 
    mutate(set = "posterior"),
  prior_sero$draws(sero_params) %>% 
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble() %>% 
    mutate(set = "prior")
) %>% 
  pivot_longer(cols = contains(sero_params)) %>% 
  mutate(param = case_when(str_detect(name, "env") ~ "constant FOI", 
                           T ~ "time-varying FOI"),
         age_cat = age_cat_dict[as.numeric(str_extract(name, "[1-3]"))],
         age_cat = factor(age_cat, levels = age_cat_dict))

p_sero_draws <- sero_param_draws %>% 
  filter(value < 1) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = set), position = "dodge") +
  facet_grid(age_cat ~ param) +
  scale_fill_manual(values = c(this_red, this_blue)) +
  theme_bw() +
  scale_x_log10()

ggsave(p_sero_draws, filename = "figures/SF_params_sero.png",
       width = 8, height = 5, dpi = 300)


# Demographic pyramids ---------------------------------------------------------

survey <- read_csv("data/demographic_pyramid_survey_data.csv")
clinical <- read_csv("data/demographic_pyramid_clinical_data.csv")

cdat <- read_csv("data/census_data_sex_age_sitakunda.csv") %>% 
  filter(Age != "100+") %>% 
  mutate(prop_male = 100*`Male Population`/sum(`Both Sexes Population`),
         prop_female = 100*`Female Population`/sum(`Both Sexes Population`))

pop_levels <- cut(c(1:99), c(seq(0,15,by=5), seq(20, 80, by = 10), 100), right=FALSE, include.lowest = TRUE) %>% ordered() %>% levels()

pop_labels <- c("0-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

cdat <- cdat %>% mutate(Age = str_replace(Age, "\\-0", "\\-"),
                        age = case_when(Age %in% c("0-4", "5-9", "10-14", "15-19", "80+") ~ Age, 
                                        TRUE ~ paste0(str_extract(Age, "^\\d"), "0-",
                                                      str_extract(Age, "^\\d"), "9"))) %>% 
  group_by(age) %>% 
  summarize(total = sum(`Both Sexes Population`), 
            male = sum(`Male Population`), 
            female = sum(`Female Population`)) %>% 
  mutate(age = factor(age, levels = pop_labels)) %>%
  arrange(age) %>%
  mutate(prop_male = -male/sum(male)*100, 
         prop_female = female/sum(female)*100)

pyramid <- bind_rows(
  survey %>% 
    filter(num_rounds == "Serological Data") %>%     
    dplyr::select(cat_sex, age = n_age) %>%
    mutate(type = "A"), 
  clinical %>% 
    dplyr::select(cat_sex, age = n_criteria_age) %>% 
    mutate(type = "B"), 
  clinical %>% 
    filter(rdt_result=="Positive") %>% 
    dplyr::select(cat_sex, age = n_criteria_age) %>% 
    mutate(type = "C")
) %>% 
  mutate(cat_age = cut(age, breaks = c(seq(0,15,by=5), seq(20, 80, by = 10), 100), 
                       right=FALSE, include.lowest = TRUE, labels = pop_labels)) %>%
  count(type, cat_sex, cat_age) %>% 
  dplyr::group_by(type, cat_sex) %>% 
  mutate(n = n/sum(n)*100) %>% 
  mutate(n = if_else(cat_sex=="Male"|cat_sex=="male", -n, n)) 

pyramid_a <- pyramid %>% filter(type == "A") %>%
  ggplot()+
  geom_col(aes(x=cat_age, y=n, fill=cat_sex))+
  geom_point(data = cdat, aes(x=age, y=prop_male),color = "black")+
  geom_point(data = cdat, aes(x=age, y=prop_female),color = "black")+
  coord_flip()+
  scale_y_continuous(labels=abs, breaks = c(-20, -15, -10, -5, 0, 5, 10, 15, 20))+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "bottom", panel.grid.minor = element_blank())+
  xlab("Age Category")+
  ylab("Proportion of survey participants")+
  scale_fill_manual(values = c("steelblue", "pink"))

pyramid_b <- pyramid %>% filter(type == "B") %>%
  ggplot()+
  geom_col(aes(x=cat_age, y=n, fill=cat_sex))+
  geom_point(data = cdat, aes(x=age, y=prop_male),color = "black")+
  geom_point(data = cdat, aes(x=age, y=prop_female),color = "black")+
  coord_flip()+
  scale_y_continuous(labels=abs, breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40))+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "bottom", panel.grid.minor = element_blank())+
  xlab("Age Category")+
  ylab("Proportion of suspected cases (clinical surveillance)")+
  scale_fill_manual(values = c("steelblue", "pink"))

pyramid_c <- pyramid %>% filter(type == "C") %>%
  ggplot()+
  geom_col(aes(x=cat_age, y=n, fill=cat_sex))+
  geom_point(data = cdat, aes(x=age, y=prop_male),color = "black")+
  geom_point(data = cdat, aes(x=age, y=prop_female),color = "black")+
  coord_flip()+
  scale_y_continuous(labels=abs, breaks = c(-30, -20, -10, 0, 10, 20, 30))+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "bottom", panel.grid.minor = element_blank())+
  xlab("Age Category")+
  ylab("Proportion of RDT-confirmed cases (clinical surveillance)")+
  scale_fill_manual(values = c("steelblue", "pink"))

pyramid_clinical <- cowplot::plot_grid(pyramid_b+theme(legend.position="none"), pyramid_c+theme(axis.title.y = element_blank(), legend.position="none"), 
                                       rel_widths = c(0.52, 0.48), 
                                       labels = c("B", "C"))    

pyramid_survey <- cowplot::plot_grid(pyramid_a+ theme(legend.position="none"),
                                     rel_widths = c(0.52, 0.48), 
                                     ncol = 2, 
                                     labels = c("A", ""))

pyramid_final2 <- cowplot::plot_grid(pyramid_survey, 
                                     pyramid_clinical,
                                     ncol=1)

png("figures/pyramid2.png", width = 2302, height = 1740, res = 200)
cowplot::plot_grid(pyramid_final2, 
                   cowplot::get_legend(pyramid_b), 
                   ncol=1,
                   rel_heights = c(0.95, 0.05))
dev.off()


