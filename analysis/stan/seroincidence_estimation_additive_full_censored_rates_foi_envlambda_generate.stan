// This Stan script aims at inferring seroconversion across three rounds
// assuming and exponential decay model for antibody kinetics
//
// The model assumes that prior infection status is unknown. 
// There are eight possible states for each participant:
// 1. Not-infected prior & not-infected round 1 & not-infected round 2
// 2. Not-infected prior & not-infected round 1 & infected round 2
// 3. Not-infected prior & infected round 1 & not-infected round 2
// 4. Not-infected prior & infected round 1 & infected round 2
// 5. Infected prior & not-infected round 1 & not-infectedion round 2
// 6. Infected prior & not-infected round 1 & infected round 2
// 7. Infected prior & infected round 1 & not-infected round 2
// 8. Infected prior & infected round 1 & infected round 2
// 
// This version of the model assumes interval-censored titer observations
// 
// This version of the model uses infection rates instead of probabilities by
// serology round.
// 

data {
  int<lower=1> N;    // Number of participants
  int<lower=1> J;    // Number of samples of decay parameters
  int<lower=1> M[3];    // Number of time slices at which infection may have occurred (can be days but probably better to be more coarse for computation)
  int<lower=1> K;       // number of dilutions in series
  int<lower=1> N_age_cat;   // number of age categories
  
  real<lower=0> dilutions[K];          // sequence of dilutions for antibody titers
  real y[N, 3];      // antibody titer observation at baseline and followups. This matrice contains the dilution factor index.
  int<lower=1, upper=K> y_titer[N, 3];      // antibody titer observation at baseline and followups. This matrice contains the dilution factor index.
  real<lower=0> t_follow[N,2];     // time from baseline to visit (baseline=0)
  
  real mdt_pre0[N, M[1]];     // time difference between mid-point of M and meaasurement at t_0 (this is positive)
  real mdt_0_1[N, M[2]];     // time difference between mid-point of M and meaasurement at t_1
  real mdt_1_2[N, M[3]];     // time difference between mid-point of M and meaasurement at t_2
  
  real rprob_inf_pre0[N, M[1]];    // relative probability of being infected (relative to each person's time slice!!!!)
  real rprob_inf_0_1[N, M[2]];     // relative probability of being infected in each time slice [assuming this is proportional to clinical data]
  real rprob_inf_1_2[N, M[3]];     // relative probability of being infected in each time slice [assuming this is proportional to clinical data]
  
  // Boosting decay model parameters
  real<lower=0> gamma[J];    // baseline level prior to infection
  real<lower=0> beta[J];     // boosting level parameter 
  real<lower=0> alpha[J];    // boosting decay parameter
  real<lower=0> delta[J];    // boosting time delay parameter
  real<lower=0> sigma;       // sd of measurement error (log-scale)
  
  // Data for inference of chol/nonchol incidence
  int<lower=1> J_foi;
  int<lower=1> N_incid;
  matrix[N_incid, J_foi] I_chol[N_age_cat];          // Inicidence of cholera AWD
  
  int<lower=1, upper=N_incid> foi_tl[N, 3];    // left time bound to compute foi for each round
  int<lower=1, upper=N_incid> foi_tr[N, 3];    // right time bound to compute foi for each round
  
  // Output 
  int<lower=1> N_output;  
  int<lower=1, upper=N_incid> foi_tl_output[N_output];    // left time bound to compute foi for each round
  int<lower=1, upper=N_incid> foi_tr_output[N_output];    // right time bound to compute foi for each round
}

transformed data {
  matrix[N_output, J_foi] cum_foi_output[N_age_cat];
  real<lower=0> dt_output[N_output];               // Time difference between output bounds
  
  // Cumulative foi
  for (i in 1:N_output) {
    dt_output[i] = foi_tr_output[i]-foi_tl_output[i]+1;
    for (j in 1:N_age_cat) {
      
      for (k in 1:J_foi) {
        cum_foi_output[j][i, k] = sum(I_chol[j][foi_tl_output[i]:foi_tr_output[i], k]);
      }
    }
  }
}

parameters {
  vector[N_age_cat] log_lambda;   // rate of infection multiplier
  vector[N_age_cat] log_lambda_env;
  real mu_log_lambda;
  real<lower=0> sd_log_lambda;
  real mu_log_lambda_env;
  real<lower=0> sd_log_lambda_env;
}

transformed parameters {
  vector[N_age_cat] lambda = exp(log_lambda);
  vector[N_age_cat] lambda_env = exp(log_lambda_env);
}

generated quantities {
  matrix[N_output, N_age_cat] mu_output_foi;      // Inferred probabilities of interest
  matrix[N_output, N_age_cat] mu_output_total;    // Inferred probabilities of interest
  matrix[N_output, N_age_cat] marginal_ratio;     // Ratio of marginal effects of time-varying vs. constant
  
  for (i in 1:N_output) {
    for (j in 1:N_age_cat) {
      mu_output_total[i, j] = mean(1 - exp(-lambda[j] * cum_foi_output[j][i, ] - lambda_env[j] * dt_output[i]));
      mu_output_foi[i, j] = mean(1 - exp(-lambda[j] * cum_foi_output[j][i, ]));
      marginal_ratio[i, j] = mean(lambda[j] * cum_foi_output[j][i, ]/(lambda_env[j] * dt_output[i]));
    }
  }
}
