# This script compares model fits between different seroincidence formulations

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)
library(here)

getLL <- function(chain_file) {
  print(chain_file)
  chain <- try(read_cmdstan_csv(chain_file, variables = "case_lik"))
  
  if (!inherits(chain, "try-error")) {
    chain$post_warmup_draws %>% 
      as_draws() %>% 
      as_draws_df() %>% 
      as_tibble()
  }
}

# Load fits ---------------------------------------------------------------
all_draw_files <- dir("generated_data/", pattern = "draws", full.names = T)

# Cnst FOI and full pulling
draws_cnst_fullpool <- all_draw_files %>% 
  str_subset("20230510") %>% 
  map_df(~ getLL(.)) 

# Cnst FOI and no pulling
draws_cnst_nopool <- all_draw_files %>% 
  str_subset("20230512") %>% 
  str_subset("joint_nopooling") %>% 
  map_df(~ getLL(.)) 

# No cnst FOI and pulling
draws_nocnst_pool <- all_draw_files %>% 
  str_subset("20230512") %>% 
  str_subset("novenv_foipooling")%>% 
  map_df(~ getLL(.)) 

# No cnst FOI and no pulling
draws_nocnst_nopool <- all_draw_files %>% 
  str_subset("20230512") %>% 
  str_subset("novenv_nopooling")%>% 
  map_df(~ getLL(.)) 



# Make ll matrices --------------------------------------------------------
comb_draws <- list(
  draws_cnst_fullpool,
  draws_cnst_nopool,
  draws_nocnst_pool,
  draws_nocnst_nopool
)

max_draws <- comb_draws %>% 
    map_dbl(~ filter(., !is.na(`case_lik[1]`)) %>% nrow()) %>% 
  min()

loo_list <- map(comb_draws, ~ filter(., !is.na(`case_lik[1]`)) %>% 
                  .[1:max_draws,] %>% 
                  as.matrix() %>% 
                  loo()) 

loo_comp <- loo_list %>% 
  set_names(c("cnstfoi_pool",
              "cnstfoi_nopool",
              "nocnstfoi_pool",
              "nocnstfoi_nopool")) %>% 
  loo::loo_compare()

saveRDS(loo_comp, file = here("generated_data/model_comparison.rds"))
