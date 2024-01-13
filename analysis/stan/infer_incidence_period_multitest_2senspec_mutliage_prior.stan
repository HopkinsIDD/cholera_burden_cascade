// This Stan script aims at inferring time-varying cholera incidence based
// on suspected and confirmed cases.
// 
// The model assumes that log-incidence rates of cholera and non-cholera AWD
// follow a Brownian motion.
// 
// This version of the model uses aggregated data over user-specified periods
// for AWD vs. cholera binomial draws.
// 
// This version of the model uses multiple test results including the RDT, the PCR
// and culture results.
// 

functions {
  // @title: Autocorrelation prior
  // @description Prior to model autocorrelated random vector following
  //              Sorbye and Rue (2014) https://doi.org/10.1016/j.spasta.2013.06.004
  // @param t: time to observation
  // @param gamma: baseline value
  // @param delta: boost delay
  // @param order
  real autocorrPrior(int n, real tau, vector x, int order) {
    real ll;
    
    if (order == 1) {
      // First order autoregressive model
      ll = (n-1.0)/2.0 * log(tau) - tau/2 * (dot_self(x[2:n] - x[1:(n-1)]));
      // Soft sum to 0 constraint
      ll += normal_lpdf(sum(x) | 0, 0.001 * n);
    } else {
      vector[n] ix;
      for (i in 1:n) {
        ix[i] = i * x[i];
      }
      // Second order autoregressive model
      ll = (n-2.0)/2.0 * log(tau) - tau/2 * (dot_self(x[1:(n-2)] - 2*x[2:(n-1)] + x[3:n]));
      // Soft sum to 0 constraint
      ll += normal_lpdf(sum(x) | 0, 0.001 * n);
      // Soft sum to 0 constraint on sum_i ix
      ll += normal_lpdf(sum(ix) | 0, 0.001 * n);
    }
    
    return(ll);
  }
}

data {
  int<lower=0> N_obs;       // total number of observation periods
  int<lower=0> N_incid;     // number of observed incidence time slices
  int<lower=0> N_agecat;    // number of age categories
  
  array[N_obs, N_agecat] int n_awd;    // number of suspected AWD cases
  // int n_chol[N_obs, N_agecat];    // number of positive tests
  
  // There are three testing situations:
  // A. RDT only (4 in 5 RDT-)
  // B. RDT + PCR (1 in 5 RDT-)
  // C. RDT + PCR + CUL (all RDT+)
  
  array[N_obs, N_agecat] int y_A; // number of test results when only RDT available (only the case for 4 in 5 neg RDT)
  array[N_obs, 2, N_agecat] int y_B; // number of test results when both RDT and PCR available (this is only the case for 1 in 5 neg RDT):
  //   y_B[i, 1]: number of samples with RDT- and PCR+
  //   y_B[i, 2]: number of samples with RDT- and PCR-
  array[N_obs, 4, N_agecat] int y_C; // number of test results when RDT, PCR and Culture available
  
  // int n_tests[N_obs, 3]; // number of tests peformed in testing setting
  
  array[N_obs] int<lower=1, upper=N_incid> map_obs_incid;
  
  array[4] real logit_sens_prior_mu;    // chol test sensitivity prior{RDT first half, RDT second half, PCR, culture}
  array[4] real<lower=0> logit_sens_prior_sd;   
  array[4] real logit_spec_prior_mu;    // chol test specificity prior{RDT first half, RDT second half, PCR, culture}
  array[4] real<lower=0> logit_spec_prior_sd;   
  
  real<lower=0> sd_sigma_autocorr;    // sd of prior on sd of Brownian motion
  int<lower=1, upper=2> autocorr_order;     // order of the autogressive model 
  real mu_beta0_chol;
  real mu_beta0_nonchol;
  real<lower=0> sd_beta0_chol;
  real<lower=0> sd_beta0_nonchol;
  
  // Periods to check for positivity
  int<lower=1> N_pos_periods;
  array[N_pos_periods] int<lower=1, upper=N_obs> pos_periods_tl;
  array[N_pos_periods] int<lower=1, upper=N_obs> pos_periods_tr;
  array[N_pos_periods] int<lower=1, upper=N_incid> pos_periods_incid_tl;    // bounds in full incidence vector
  array[N_pos_periods] int<lower=1, upper=N_incid> pos_periods_incid_tr;
  array[N_pos_periods] int<lower=1, upper=2> map_pos_period_sens;    // map from surveillance period to RDT batch
  
  // Prior of cholera for each age class
  array[N_agecat] real<lower=0, upper=1> prior_chol;
  
}
transformed data {
  array[N_pos_periods, N_agecat] int n_tot_period;
  array[N_pos_periods, N_agecat] int y_period_A;       // number of test results when only RDT available
  array[N_pos_periods, 2, N_agecat] int y_period_B;    // number of test results when both RDT and PCR available
  array[N_pos_periods, 8, N_agecat] int y_period_C;    // number of test results when RDT, PCR and Culture available
  array[N_pos_periods, N_agecat] real frac_B;           // fraction of RDT- samples that were tested for PCR
  array[4, 8] int fake_data_rdt;       // this is for unknwon pcr and culture tests
  array[4, N_agecat] real priors_A;    // priors for culture results when unavailable 
  array[4, N_agecat] real priors_B;    // priors for PCR and culture results when unavailable
  array[4, N_agecat] real log_priors_A;
  array[4, N_agecat] real log_priors_B;
  array[4] real sens_prior_mu;        // indices: 1:RDT first period, 2:RDT second period, 3:PCR, 4:culture
  array[4] real spec_prior_mu;
  
  
  for (i in 1:4) {
    sens_prior_mu[i] = inv_logit(logit_sens_prior_mu[i]);
    spec_prior_mu[i] = inv_logit(logit_spec_prior_mu[i]);
  }
  
  for (i in 1:N_agecat) {
    // Prior probability of a neg RDT
    real p_rdtneg = (1-sens_prior_mu[2]) * prior_chol[i] + spec_prior_mu[1] * (1-prior_chol[i]);
    // Prior probability of a neg RDT and pos PCR
    real p_rdtneg_pcrpos = (1-sens_prior_mu[2]) * sens_prior_mu[3] * prior_chol[i] + spec_prior_mu[1] * (1-spec_prior_mu[3]) * (1-prior_chol[i]);
    // Prior probability of a neg RDT and neg PCR
    real p_rdtneg_pcrneg = (1-sens_prior_mu[2]) * (1-sens_prior_mu[3]) * prior_chol[i] + spec_prior_mu[1] * spec_prior_mu[3] * (1-prior_chol[i]);
    
    // Case A
    // 1. RDT- PCR+ CUL+
    priors_A[1, i] = sens_prior_mu[3] * sens_prior_mu[4] * (1-sens_prior_mu[2]) * prior_chol[i] +
    (1-spec_prior_mu[3]) * (1-spec_prior_mu[4]) * spec_prior_mu[1] * (1-prior_chol[i]);
    // 2. RDT- PCR+ CUL-
    priors_A[2, i] = sens_prior_mu[3] * (1-sens_prior_mu[4]) * (1-sens_prior_mu[2]) * prior_chol[i] +
    (1-spec_prior_mu[3]) * spec_prior_mu[4] * spec_prior_mu[1] * (1-prior_chol[i]);
    // 3. RDT- PCR- CUL+
    priors_A[3, i] = (1-sens_prior_mu[3]) * sens_prior_mu[4] * (1-sens_prior_mu[2]) * prior_chol[i] +
    spec_prior_mu[3] * (1-spec_prior_mu[4]) * spec_prior_mu[1] * (1-prior_chol[i]); 
    // 4. RDT- PCR- CUL-
    priors_A[4, i] = (1-sens_prior_mu[3]) * (1-sens_prior_mu[4]) * (1-sens_prior_mu[2]) * prior_chol[i] +
    spec_prior_mu[3] * spec_prior_mu[4] * spec_prior_mu[1] * (1-prior_chol[i]);
    
    
    // Case B
    // Initialize at unnormalized values
    for (k in 1:4) {
      priors_B[k, i] = priors_A[k, i];
    }
    
    // Normalize by respective values
    for (k in 1:4) {
      priors_A[k, i] = priors_A[k, i]/p_rdtneg;
      if (k <= 2) {
        priors_B[k, i] = priors_B[k, i]/p_rdtneg_pcrpos;
      } else {
        priors_B[k, i] = priors_B[k, i]/p_rdtneg_pcrneg;
      }
      
      // Pre-compute logs
      log_priors_A[k, i] = log(priors_A[k, i]);
      log_priors_B[k, i] = log(priors_B[k, i]);
    }
  }
  
  // Initialize fake data for rdt only measurements
  // We know these are only for negative samples
  for (i in 1:4) {
    for (j in 1:8) {
      fake_data_rdt[i, j] = 0;
    }
    fake_data_rdt[i, i+4] = 1;
  }
  
  for (i in 1:N_pos_periods) {
    for (j in 1:N_agecat) {
      n_tot_period[i, j] = sum(n_awd[pos_periods_tl[i]:pos_periods_tr[i], j]);
      y_period_A[i, j] = sum(y_A[pos_periods_tl[i]:pos_periods_tr[i], j]);
      for (k in 1:2) {
        y_period_B[i, k, j] = sum(y_B[pos_periods_tl[i]:pos_periods_tr[i], k, j]);
      }
      for (k in 1:4) {
        y_period_C[i, k, j] = sum(y_C[pos_periods_tl[i]:pos_periods_tr[i], k, j]);
        y_period_C[i, k+4, j] = 0;
      }
      {
        real num = sum(y_period_B[i, , j]);
        real denom = num + y_period_A[i, j];
        if (denom > 0) {
          frac_B[i, j] = num/denom;
        } else {
          frac_B[i, j] = 0;
        }
      }
    }
  }
}
parameters {
  matrix[N_incid, N_agecat] log_beta_chol; // time/age specific incidence rates
  matrix[N_incid, N_agecat] log_beta_nonchol;
  array[N_agecat] real alpha_chol;        // age-specific intercept of cholera log_rate
  array[N_agecat] real alpha_nonchol;     // age-specific intercept of non-cholera log_rate
  array[N_agecat] real<lower=0> sigma_chol;       // brownian motion sd
  array[N_agecat] real<lower=0> sigma_nonchol;    // brownian motion sd
  array[4] real logit_sens;
  array[4] real logit_spec;
}

transformed parameters {
  array[N_agecat] real tau_chol; // precision of the brownian motion
  array[N_agecat] real tau_nonchol; // precision of the brownian motion
  // 
  // matrix[N_incid, N_agecat] I_chol;
  // matrix[N_incid, N_agecat] I_nonchol;
  // matrix<lower=0, upper=1>[N_incid, N_agecat] gamma;       // True probability of observing cholera
  // matrix<lower=0>[N_incid, N_agecat] I_tot;            // Total AWD incidence
  // matrix<lower=0>[N_pos_periods, N_agecat] gamma_avg;      // average true probability of positivity
  // real<lower=0, upper=1> p[N_pos_periods, 8, N_agecat];     // Probability of observing cholera accounting for test performance
  array[4] real sens_all = inv_logit(logit_sens);
  array[4] real spec_all = inv_logit(logit_spec);
  
  for (i in 1:N_agecat) {
    tau_chol[i] = 1/sigma_chol[i];
    tau_nonchol[i] = 1/sigma_nonchol[i];
    // I_chol[,i] = exp(alpha_chol[i] + log_beta_chol[,i]);
    // I_nonchol[,i] = exp(alpha_nonchol[i] + log_beta_nonchol[,i]);
    // 
    // // Incidence
    // I_tot[,i] = I_chol[,i] + I_nonchol[,i];
    // 
    // // Observation probability (proportion of awd that is cholera)
    // gamma[,i] = I_chol[,i] ./(I_chol[,i] + I_nonchol[,i]);
    // 
    // // Compute average
    // for (j in 1:N_pos_periods) {
    //   gamma_avg[j, i] = mean(gamma[pos_periods_incid_tl[j]:pos_periods_incid_tr[j], i]);
    // }
    
    // Compute probabilities for each test result combination
    // We here use all three test type results:
    // 1. RDT+ PCR+ CUL+
    // 2. RTD+ PCR+ CUL-
    // 3. RDT+ PCR- CUL+
    // 4. RDT+ PCR- CUL-
    // 5. RDT- PCR+ CUL+
    // 6. RDT- PCR+ CUL-
    // 7. RDT- PCR- CUL+
    // 8. RDT- PCR- CUL-
    // for (j in 1:N_pos_periods) {
    //   row_vector[3] sens = [sens_all[map_pos_period_sens[j]], sens_all[3], sens_all[4]];
    //   row_vector[3] spec = [spec_all[map_pos_period_sens[j]], spec_all[3], spec_all[4]];
    //   real x = gamma_avg[j, i];    // cache to true infection probability
    //   p[j, 1, i] = prod(sens) * x + prod(1-spec) * (1-x);
    //   p[j, 2, i] = sens[1]*sens[2]*(1-sens[3]) * x + (1-spec[1])*(1-spec[2])*spec[3] * (1-x);
    //   p[j, 3, i] = sens[1]*(1-sens[2])*sens[3] * x + (1-spec[1])*spec[2]*(1-spec[3]) * (1-x);
    //   p[j, 4, i] = sens[1]*(1-sens[2])*(1-sens[3]) * x + (1-spec[1])*spec[2]*spec[3] * (1-x);
    //   p[j, 5, i] = (1-sens[1])*sens[2]*sens[3] * x + spec[1]*(1-spec[2])*(1-spec[3]) * (1-x);
    //   p[j, 6, i] = (1-sens[1])*sens[2]*(1-sens[3]) * x + spec[1]*(1-spec[2])*spec[3] * (1-x);
    //   p[j, 7, i] = (1-sens[1])*(1-sens[2])*sens[3] * x + spec[1]*spec[2]*(1-spec[3]) * (1-x);
    //   p[j, 8, i] = prod(1-sens) * x + prod(spec) * (1-x);
    // }
  }
}

model {
  
  // // Multinomial log-likelihoods
  // for (i in 1:N_pos_periods) {
  //   for (j in 1:N_agecat) {
  //     array[4] real p_rdt;
  //     matrix[2, 2] p_rdt_pcr;
  //     vector[8] p_vec = to_vector(p[i, ,j]);
  //     // note: remember that we are treating missing data are non-inofrmative which isn't quite true
  //     // Likelihood of RDT only observations (missing PCR and culture)
  //     for (k in 1:4) {
  //       p_rdt[k] = multinomial_lpmf(fake_data_rdt[k, ]| p_vec) + log_priors_A[k, j];
  //     }
  //     target += y_period_A[i, j] * log_sum_exp(p_rdt);
  //     
  //     // Likelihood of RDT and PCR observations (missing culture)
  //     for (k in 1:2) {
  //       p_rdt_pcr[k, 1] = multinomial_lpmf(fake_data_rdt[(k-1)*2+1, ]| p_vec) + log_priors_B[k, j];
  //       p_rdt_pcr[k, 2] = multinomial_lpmf(fake_data_rdt[(k-1)*2+2, ]| p_vec) + log_priors_B[k+2, j];
  //       target += y_period_B[i, k, j] * log_sum_exp(p_rdt_pcr[k, ]);
  //     }
  //     
  //     // Complete observations (RDT and PCR and culture)
  //     target += multinomial_lpmf(y_period_C[i, ,j]| p_vec);
  //   }
  // }
  
  for (i in 1:N_agecat) {
    // // Incidence
    // n_awd[,i] ~ poisson(I_tot[map_obs_incid, i]);
    
    // Autocorrelation following Brownian motion
    target += autocorrPrior(N_incid, tau_chol[i], log_beta_chol[, i], autocorr_order);
    target += autocorrPrior(N_incid, tau_nonchol[i], log_beta_nonchol[, i], autocorr_order);
    
  }
  
  // Priors on intercept
  alpha_chol ~ normal(mu_beta0_chol, sd_beta0_chol);
  alpha_nonchol ~ normal(mu_beta0_nonchol, sd_beta0_nonchol);
  
  log_beta_chol[1,] ~ normal(0, .5);
  log_beta_nonchol[1,] ~ normal(0, .5);
  
  // Prior on sigma
  sigma_chol ~ normal(0, sd_sigma_autocorr);
  sigma_nonchol ~ normal(0, sd_sigma_autocorr);
  
  // Prior on sens/spec
  logit_sens ~ normal(logit_sens_prior_mu, logit_sens_prior_sd);
  logit_spec ~ normal(logit_spec_prior_mu, logit_spec_prior_sd);
}

generated quantities {
  // matrix[N_incid, N_agecat] n_awd_gen;
  // array[N_pos_periods, 8, N_agecat] real y_gen;
  // matrix[N_pos_periods, N_agecat] gamma_avg_output = gamma_avg;
  // matrix[N_pos_periods, N_agecat] y_A_gen;
  // array[N_pos_periods, 2, N_agecat] real y_B_gen;
  // array[N_pos_periods, 4, N_agecat] real y_C_gen;
  // 
  // for (j in 1:N_agecat) {
  //   for (i in 1:N_incid) {
  //     n_awd_gen[i, j] = poisson_rng(I_tot[i, j]);
  //   }
  //   
  //   for (i in 1:N_pos_periods) {
  //     if (n_tot_period[i, j] > 0) {
  //       y_gen[i, ,j] = multinomial_rng(to_vector(p[i, , j]), n_tot_period[i, j]);
  //     } else {
  //       for (k in 1:8) {
  //         y_gen[i, k, j] = 0;
  //       }
  //     }
  //     y_A_gen[i, j] = sum(y_gen[i, 5:8, j]) * (1-frac_B[i, j]);     // This is 4 in 5 of all RDT-
  //     y_B_gen[i, 1, j] = sum(y_gen[i, 5:6, j]) * frac_B[i, j];    // This is 1 in 5 of all RDT- that also has an PCR+ test
  //     y_B_gen[i, 2, j] = sum(y_gen[i, 7:8, j]) * frac_B[i, j];    // This is 1 in 5 of all RDT- that also has an PCR- test
  //     y_C_gen[i, , j] = y_gen[i, 1:4, j];
  //   }
  // }
  // 
}
