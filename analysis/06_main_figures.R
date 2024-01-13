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


files <- dir("generated_data/", pattern = "for_figures_diarrhea2", full.name = T)
# files <- dir("generated_data/", pattern = "for_figure", full.name = T)
for (f in files) {
  load(f)
}


# Figure specs
caselim <- c(0, 45)
t_min_plot <- as.Date("2021-01-24") # tmin for plotting
t_max_plot <- getTmax()   # tmin for plotting
datelim <- c(t_min_plot, t_max_plot)
date_subset <- seq.Date(t_min_plot, t_max_plot, by = "1 days")

# Colors
this_blue <- "#493D99"
this_red <- "#BD3E5C"

# Update data -------------------------------------------------------------


# Modify for plotting
lab_data <- lab_data %>% 
  mutate(rdt_pcr = str_c("RDT", ifelse(rdt_res, "+", "-"), " & PCR", ifelse(pcr_res, "+", "-")),
         rdt_pcr = ifelse(is.na(rdt_pcr), "RDT- & no PCR", rdt_pcr),
         rdt_pcr = factor(rdt_pcr, levels = c("RDT+ & PCR+",
                                              "RDT+ & PCR-",
                                              "RDT- & PCR+",
                                              "RDT- & PCR-",
                                              "RDT- & no PCR") %>% 
                            rev())) 

# Figure 2 - no agecats ---------------------------------------------------

p_data <- lab_data %>%
  ggplot(aes(x = tmid)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_bar(aes(fill = rdt_pcr,
               x = tmid), 
           inherit.aes = F,
           width = 6) +
  geom_label(data = round_dates,
             inherit.aes = F,
             aes(label = round, x = tmid, y = 40),
             size = 2) +
  scale_fill_manual(values = rev(c("#800ABF", "#CF8FF7", "#857E7E", "#DE6B0D", "#F7973E"))) +
  guides(fill = guide_legend("RDT and PCR result")) +
  theme_bw() +
  labs(y = "Weekly count", x = "Date")  +
  theme_bw()  +
  coord_cartesian(xlim = datelim) 

p_infer_sum <- weekly_I_sum_true %>% 
  mutate(set = "community") %>% 
  bind_rows(weekly_I_sum %>% mutate(set = "clinical")) %>% 
  filter(tmid >= min(lab_data$tmid)) %>% 
  ggplot(aes(x = tmid, y  = mean, ymin = q5, ymax = q95,
             group = set)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_errorbar(size = .35, width = 0, color = "black") +
  geom_point(aes(pch = set), fill = "gray", stroke = .35) +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend("Symptomatic infections")) +
  theme_bw() +
  labs(y = "Weekly symptomatic infections", x = "Date")  +
  coord_cartesian(xlim = datelim, ylim = caselim)

p_exp_sum <- tot_exposed %>% 
  filter(tmid >= min(lab_data$tmid)) %>% 
  ggplot(aes(x = tmid, y  = mean, ymin = q5, ymax = q95)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_errorbar(size = .35, width = 0, color = "black") +
  geom_point(pch = 22, fill = "gray", stroke = .35) +
  theme_bw() +
  labs(y = "Weekly exposures", x = "Date")  +
  coord_cartesian(xlim = datelim, ylim = c(0, 20000))


remove_x_axis <- theme(axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       axis.ticks.x = element_blank())

scale_dates_axis <- scale_x_date(date_breaks = "2 months",
                                 date_labels = "%Y-%m")

p_figure2_noage <- cowplot::plot_grid(
  p_data +
    remove_x_axis +
    scale_dates_axis,
  p_infer_sum +
    remove_x_axis +
    scale_dates_axis,
  p_exp_sum +
    scale_dates_axis,
  ncol = 1, align = "v", axis = "lr")

ggsave(p_figure2_noage, 
       filename = "figures/figure2_noage.png", 
       width = 8, height = 7)

ggsave(p_figure2_noage, 
       filename = "figures/figure2_noage.pdf", 
       width = 8, height = 7)

# weekly infections and annualized incidence when no medically attended cholera cases occur
epiweek_no_clinical_cases <- weekly_I_sum %>% 
  filter(tmid >= min(lab_data$tmid)) %>% 
  filter(mean<1) 

epiweek_infections_no_clinical_cases <- tot_exposed %>% 
  filter(tmid >= min(lab_data$tmid)) %>%
  addEpiWeek(., date_col = "tmid") %>% 
  filter(epiweek %in% epiweek_no_clinical_cases$epiweek)

infections_no_clinical_cases <- range(epiweek_infections_no_clinical_cases$mean)

poi_dt <- getOutputPeriodOverall() %>% # length of modeling period
  unlist() %>% 
  diff()

pop_overall <- 450094.9

incidence_annualized_no_clinical_cases <- (infections_no_clinical_cases/pop_overall)*(365.25/poi_dt)*1e3

# Figure 2 - with agecats ------------------------------------------------

age_cat_dict <- getAgeCatDictFigures() 

p_data_age <- lab_data %>% 
  ggplot(aes(x = tmid)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_bar(aes(fill = rdt_pcr,
               x = tmid), 
           inherit.aes = F,
           width = 6) +
  geom_label(data = round_dates,
             inherit.aes = F,
             aes(label = round, x = tmid, y = 35),
             size = 2) +
  facet_grid(.~age_cat, labeller = labeller(age_cat = age_cat_dict)) +
  scale_fill_manual(values = rev(c("#800ABF", "#CF8FF7", "#857E7E", "#DE6B0D", "#F7973E"))) +
  guides(fill = guide_legend("RDT and PCR result")) +
  theme_bw() +
  labs(y = "Weekly count", x = "Date")  +
  theme_bw()  +
  coord_cartesian(xlim = datelim) +
  theme(legend.position = c(.89, .6),
        legend.key.size = unit(.15, "in"),
        legend.title = element_text(size = 10),
        legend.box.background = element_rect(fill="white", size=0, color="white"),
        legend.box.margin = margin(1, 8, 1, 1))

p_infer_sum_age <- weekly_I_sum_age_true %>% 
  mutate(set = "Overall community\nsymptomatic infections") %>% 
  bind_rows(weekly_I_sum_age %>% mutate(set = "Medically-attented\ncholera cases")) %>% 
  filter(tmid >= min(lab_data$tmid)) %>% 
  ggplot(aes(x = tmid, y  = mean, ymin = q5, ymax = q95,
             group = set)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_errorbar(size = .35, width = 0, color = "black", alpha = .5) +
  geom_point(aes(pch = set), fill = "gray", stroke = .35) +
  facet_grid(.~age_cat, labeller = labeller(age_cat = age_cat_dict)) +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend("Cholera cases", byrow = TRUE)) +
  theme_bw() +
  labs(y = "Weekly incidence", x = "Date")  +
  coord_cartesian(xlim = datelim, ylim = caselim) +
  theme(legend.position = c(.87, .6),
        legend.key.size = unit(.15, "in"),
        legend.title = element_text(size = 10),
        legend.box.background = element_rect(fill="white", size=0, color="white"),
        legend.box.margin = margin(1, 15, 1, 1),
        legend.spacing.y = unit(.07, "in"))

p_exp_sum_age <- tot_exposed_age %>%
  inner_join(pop_dat) %>% 
  mutate(across(c("mean", "q5", "q95"), .fns = ~ .x/pop*1000)) %>% 
  filter(age_cat != "overall") %>% 
  filter(tmid >= min(lab_data$tmid)) %>% 
  ggplot(aes(x = tmid, y  = mean, ymin = q5, ymax = q95)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_errorbar(size = .35, width = 0, color = "black") +
  geom_point(pch = 22, fill = "gray", stroke = .35) +
  facet_grid(.~age_cat, labeller = labeller(age_cat = age_cat_dict)) +
  theme_bw() +
  labs(y = "Weekly incidence rate\n per 1,000 person", x = "Date")#  +
  # coord_cartesian(xlim = datelim, ylim = c(0, 20000))

p_figure2_age <- cowplot::plot_grid(
  p_data_age +
    remove_x_axis +
    scale_dates_axis,
  p_infer_sum_age +
    remove_x_axis +
    scale_dates_axis,
  p_exp_sum_age +
    scale_dates_axis +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1, 1.2))

ggsave(p_figure2_age, 
       filename = "figures/figure2_age.png", 
       width = 10, height = 8)

ggsave(p_figure2_age, 
       filename = "figures/figure2_age.pdf", 
       width = 10, height = 8)


# Figure 2 symptomatic infections free scale inset - with agecats ------------------------------------------------

p_infer_sum_age_freescale <- weekly_I_sum_age_true %>% 
  mutate(set = "Overall community\nsymptomatic infections") %>% 
  bind_rows(weekly_I_sum_age %>% mutate(set = "Medically-attented\ncholera cases")) %>% 
  filter(tmid >= min(lab_data$tmid)) %>% 
  ggplot(aes(x = tmid, y  = mean, ymin = q5, ymax = q95,
             group = set)) +
  geom_rect(data = round_dates,
            inherit.aes = F,
            aes(xmin = tl, xmax = tr, group = round, ymin = -Inf, max = Inf),
            fill = "gray", alpha = .3) +
  geom_errorbar(size = .35, width = 0, color = "black", alpha = .5) +
  geom_point(aes(pch = set), fill = "gray", stroke = .35) +
  coord_cartesian(xlim = datelim) +
  facet_grid(age_cat~., scales = "free", labeller = labeller(age_cat = age_cat_dict)) +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend("Cholera cases", byrow = TRUE)) +
  theme_bw() +
  labs(y = "Weekly incidence", x = "Date") +
  theme(legend.position = c(.85, .54),
        legend.key.size = unit(.15, "in"),
        legend.title = element_text(size = 10),
        legend.box.background = element_rect(fill="white", size=0, color="white"),
        legend.box.margin = margin(1, 15, 1, 1),
        legend.spacing.y = unit(.07, "in"))


p_figure2_age_freescale <- cowplot::plot_grid(
  p_infer_sum_age_freescale +
    scale_dates_axis +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  ncol = 1, align = "v", axis = "lr")

ggsave(p_figure2_age_freescale, 
       filename = "figures/figure2_age_symptomatic_freescale.png", 
       width = 8, height = 6)

ggsave(p_figure2_age_freescale, 
       filename = "figures/figure2_age_symptomatic_freescale.pdf", 
       width = 8, height = 6)

# Figure 3 - conversion factors -------------------------------------------

# Combine numbers
num_stats <- bind_rows(
  num_case_age_stats %>% 
    mutate(var = "exp_to_case"),
  # num_culture_age_stats %>% 
  #   mutate(var = "exp_to_culture"),
  num_symp_age_stats %>% 
    mutate(var = "exp_to_symp"),
  num_seek_age_stats %>% 
    mutate(var = "symp_to_case")
)

write_csv(num_stats, "generated_data/num_stats.csv")

num_stats <-  num_stats %>% 
  mutate(label = str_c(
    formatC(ifelse(mean > 100, round(mean/10)*10, mean), digits = 0, format = "f"),
    "\n(",
    formatC(ifelse(q5 > 100, round(q5/10)*10, q5), digits = 0, format = "f"),
    "-",
    formatC(ifelse(q95 > 100, round(q95/10)*10, q95), digits = 0, format = "f"),
    ")"
  ))

num_stats <- num_stats %>% 
  mutate(age_cat = factor(age_cat_dict[age_cat], levels = age_cat_dict)) 

p_state <- num_stats %>%
  ggplot(aes(x = var, y = mean, 
             ymin = q5, ymax = q95, color = age_cat)) +
  geom_point(position = position_dodge(width = .3), key_glyph = draw_key_point) +
  geom_errorbar(position = position_dodge(width = .3), 
                width = 0, 
                key_glyph = NULL) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  ggthemes::scale_color_few("Dark") + scale_y_continuous(breaks = c(1,1000,2000,3000,4000))

for (i in 1:length(unique(num_stats$var))) {
  # https://stackoverflow.com/questions/61505768/use-position-dodge-for-arrowheads-but-nudge-x-for-label-boxes-in-ggreppelgeom
  p_state <- p_state + 
    ggrepel::geom_label_repel(data = filter(num_stats,
                                            var == unique(num_stats$var)[i],
                                            #age_cat == "overall"
                                            ),
                              aes(label = label), 
                              position = position_dodge(width = .3), 
                              size = 2,
                              min.segment.length = 0, 
                              seed = 42,
                              force = 1, 
                              direction = 'y',
                              xlim = c(i+.2, i+2),
                              segment.size = 0,
                              segment.linetype = 2,
                              segment.alpha = .5,
                              label.size = .2,
                              alpha = .8, key_glyph = NULL)
}


p_state <- p_state +
  guides(color = guide_legend("Age category")) +
  scale_x_discrete(labels = c("Infections per\nmedically-attended cholera case",
                              # "infections per\npositive\nculture",
                              "Infections per\nsymptomatic infection",
                              "Symptomatic infections per\nmedically-attended case")) +
  theme(legend.position = c(.85, .8))

ggsave(p_state, filename = "figures/figure_3.png",
       width = 7, height = 5)

ggsave(p_state, filename = "figures/figure_3.pdf",
       width = 7, height = 5)
