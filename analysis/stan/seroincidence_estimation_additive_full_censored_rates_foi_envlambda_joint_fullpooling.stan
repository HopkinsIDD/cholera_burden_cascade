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


functions {
  // @title Boost decay model
  // @description This function models the boost-exponential decay antibody dynamics
  //              Setting beta to 0 is equivalent to an exponential decay
  //              The function returns the log of the predicted antibody level
  // @param t: time to observation
  // @param gamma: baseline value
  // @param beta: boost value
  // @param alpha: decay rate
  // @param delta: boost delay
  real boostModel(real t, real gamma, real beta, real alpha, real delta) {
    real y_gen;
    real epsilon = t - delta;
    real loggamma = log(gamma);
    real logbeta = log(beta);
    if (epsilon < 0) {
      y_gen = loggamma;
    } else {
      real dy;
      dy = logbeta - alpha * epsilon;
      y_gen = log_sum_exp([loggamma, dy]);
    }
    return y_gen;
  }
  
  // @title Interval normal CDF
  // @description Function to compute the interval cdf of normal distribution.
  //              Following https://github.com/stan-dev/stan/issues/1154
  // @param mu
  // @param sigma
  // @param L
  // @param U
  // @details Use 9999 as flag for Inf upper bound
  real intervalNormalCDF(real mu, real sigma, real L, real U) {
    real intervalcdf;
    
    if (U == log(9999.0)) {
      intervalcdf = normal_lccdf(L| mu, sigma);
    } else {
      intervalcdf = log_diff_exp(normal_lcdf(U| mu, sigma), normal_lcdf(L| mu, sigma));
    }
    
    // intervalcdf = normal_lpdf(L| mu, sigma);
    
    if (is_inf(intervalcdf)) {
      intervalcdf = -1e3;
    }
    return intervalcdf;
  }
  
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
  int<lower=1> N;    // Number of participants
  int<lower=1> J;    // Number of samples of decay parameters
  int<lower=1> M[3];    // Number of time slices at which infection may have occurred (can be days but probably better to be more coarse for computation)
  int<lower=1> K;       // number of dilutions in series
  int<lower=1> N_age_cat;   // number of age categories
  
  real<lower=0> dilutions[K];          // sequence of dilutions for antibody titers
  real y[N, 3];      // antibody titer observation at baseline and followups. This matrice contains the dilution factor index.
  int<lower=1, upper=K> y_titer[N, 3];      // antibody titer observation at baseline and followups. This matrice contains the dilution factor index.
  real<lower=0> t_follow[N,2];     // time from baseline to visit (baseline=0)
  int<lower=1, upper=N_age_cat> age_cat_vec[N];
  
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
  
  real<lower=0> sd_lambda;
}

transformed data {
  real J_real = J;    // cast to real
  real logj = log(1/J_real);
  real J_real_foi = J_foi;
  real logjfoi = log(1/J_real_foi);
  
  real log_gamma[J];
  real y_L[N, 3];    // Lower bound of dilution series
  real<lower=min(dilutions), upper=max(dilutions)> y_U[N, 3];    // Upper bound of dilution series
  // Log-likelihoods of observations all combinations
  real ll_obs[N, 8];
  // Pre-computed antibody trajectories
  real y_star_pre0[N, M[1], J, 3];    // Modeled trajectories of infetion before baseline at t0, t1, t2
  real y_star_0_1[N, M[2], J, 2];     // Modeled trajectories of infetion between baseline and round 1 at t1 and t2
  real y_star_1_2[N, M[3], J];        // Modeled trajectories of infetion between baseline and round 1 at t1 and t2
  real y_star_0_1_boost[N, M[2], J, 2];     // Modeled trajectories of infetion between baseline and round 1 at t1 and t2 without baseline
  real y_star_1_2_boost[N, M[3], J];        // Modeled trajectories of infetion between baseline and round 1 at t1 and t2 without baseline
  real log_rprob_inf_pre0[N, M[1]];   
  real log_rprob_inf_0_1[N, M[2]];
  real log_rprob_inf_1_2[N, M[3]];
  real<lower=0> cum_foi[N, 3, J_foi];   // Cumulative force of infection
  real<lower=0> dt[N, 3];               // Time difference between rounds in days
  
  // Cumulative foi
  for (i in 1:N) {
    for (j in 1:3) {
      
      dt[i, j] = foi_tr[i, j]-foi_tl[i, j]+1;
      
      for (k in 1:J_foi) {
        cum_foi[i, j, k] = sum(I_chol[age_cat_vec[i]][foi_tl[i, j]:foi_tr[i, j], k]);
      }
    }
  }
  
  for (j in 1:J) {
    log_gamma[j] = log(gamma[j]);
  }
  
  // Set bounds of dilution series observations
  for (i in 1:N) {
    for (j in 1:3) {
      y_L[i, j] = dilutions[y_titer[i, j]]; 
      y_U[i, j] = dilutions[y_titer[i, j] + 1]; 
    }
  }
  
  for (i in 1:N) {
    for (m in 1:M[1]) {
      log_rprob_inf_pre0[i,m] = log(rprob_inf_pre0[i,m]);
    }
    for (m in 1:M[2]) {
      log_rprob_inf_0_1[i,m] = log(rprob_inf_0_1[i,m]);
    }
    for (m in 1:M[3]) {
      log_rprob_inf_1_2[i,m] = log(rprob_inf_1_2[i,m]);
    }
  }
  
  // Pre-compute modeled antibody levels at each serosurvey round across possible
  // Unknown infection dates
  for (j in 1:J) {
    for (i in 1:N) {
      for (m in 1:M[1]) {
        y_star_pre0[i, m, j, 1] = boostModel(mdt_pre0[i, m], gamma[j], beta[j], alpha[j], delta[j]);
        y_star_pre0[i, m, j, 2] = boostModel(t_follow[i, 1] + mdt_pre0[i, m], gamma[j], beta[j], alpha[j], delta[j]);
        y_star_pre0[i, m, j, 3] = boostModel(t_follow[i, 2] + mdt_pre0[i, m], gamma[j], beta[j], alpha[j], delta[j]);
        
      }
      for (m in 1:M[2]) {
        y_star_0_1[i, m, j, 1] = boostModel(mdt_0_1[i, m], gamma[j], beta[j], alpha[j], delta[j]);
        y_star_0_1[i, m, j, 2] = boostModel(t_follow[i, 2] - t_follow[i, 1] + mdt_0_1[i, m], gamma[j], beta[j], alpha[j], delta[j]);
        y_star_0_1_boost[i, m, j, 1] = boostModel(mdt_0_1[i, m], 0, beta[j], alpha[j], delta[j]);
        y_star_0_1_boost[i, m, j, 2] = boostModel(t_follow[i, 2] - t_follow[i, 1] + mdt_0_1[i, m], 0, beta[j], alpha[j], delta[j]);
      }
      for (m in 1:M[3]) {
        y_star_1_2[i, m, j] = boostModel(mdt_1_2[i, m], gamma[j], beta[j], alpha[j], delta[j]);
        y_star_1_2_boost[i, m, j] = boostModel(mdt_1_2[i, m], 0, beta[j], alpha[j], delta[j]);
      }
    }
  }
  
  for (i in 1:N) {
    // Log-likelihoods of observations all combinations
    real ll_obs_tmp[J, 8];
    for (j in 1:J){
      vector[M[1]] ll_tmp_0; 
      vector[M[2]] ll_tmp_1; 
      vector[M[3]] ll_tmp_2; 
      matrix[M[1], M[2]] ll_tmp_0_1; 
      matrix[M[1], M[3]] ll_tmp_0_2;
      matrix[M[2], M[3]] ll_tmp_1_2;
      real ll_tmp_0_1_2[M[1], M[2], M[3]];
      
      // Initialize observation likelihoods
      for (k in 1:8) {
        ll_obs_tmp[j, k] = logj;
      }
      
      
      // ---- 1. obs loglik of <0,0,0> ----
      // When no boosting, only account for baseline antibody levels 
      for (t in 1:3){
        ll_obs_tmp[j, 1] += intervalNormalCDF(log_gamma[j], sigma, y_L[i, t], y_U[i, t]);
      }
      
      
      // ---- 2. obs loglik of <0,0,1> ----
      // For first two observations only baseline antibody levels
      for (t in 1:2){
        ll_obs_tmp[j, 2] += intervalNormalCDF(log_gamma[j], sigma, y_L[i, t], y_U[i, t]);
      }
      
      // For infection between round 1 and round 2 account for possible
      // infection dates
      for(m in 1:M[3]){
        ll_tmp_2[m] = log_rprob_inf_1_2[i, m] + intervalNormalCDF(y_star_1_2[i, m, j], sigma, y_L[i, 3], y_U[i, 3]);
      }
      ll_obs_tmp[j, 2] += log_sum_exp(ll_tmp_2);
      
      
      // ---- 3. obs loglik of <0,1,0> ----
      // The baseline serology only based on baseline antibody level
      ll_obs_tmp[j, 3] += intervalNormalCDF(log_gamma[j], sigma, y_L[i, 1], y_U[i, 1]);
      
      // For rounds 1 and two account for boosting that occurred between baseline and round 1
      for(m in 1:M[2]){
        ll_tmp_1[m] = log_rprob_inf_0_1[i, m] + 
        intervalNormalCDF(y_star_0_1[i, m, j, 1], sigma, y_L[i, 2], y_U[i, 2]) +
        intervalNormalCDF(y_star_0_1[i, m, j, 2], sigma, y_L[i, 3], y_U[i, 3]);
      }
      ll_obs_tmp[j, 3] += log_sum_exp(ll_tmp_1);
      
      
      // ---- 4. obs loglik of <0,1,1> ----
      // The baseline serology only based on baseline antibody level
      ll_obs_tmp[j, 4] += intervalNormalCDF(log_gamma[j], sigma, y_L[i, 1], y_U[i, 1] );
      
      for(m in 1:M[2]){
        // For round 1 account for boosting that occurred between baseline and round 1
        ll_tmp_1[m] = log_rprob_inf_0_1[i, m] + 
        intervalNormalCDF(y_star_0_1[i, m, j, 1], sigma, y_L[i, 2], y_U[i, 2]);
        
        for (m2 in 1:M[3]) {
          // For round 2 account for boosting that occurred between baseline and round 1
          // as well as round 1
          real ys = log_sum_exp([y_star_0_1[i, m, j, 2], y_star_1_2_boost[i, m2, j]]);
          ll_tmp_1_2[m, m2] = log_rprob_inf_0_1[i, m] + log_rprob_inf_1_2[i, m2] +
          intervalNormalCDF(ys, sigma, y_L[i, 3], y_U[i, 3]);
        }
      }
      
      {
        // intermediate computation
        vector[M[2]] x;
        for (m in 1:M[2]) {
          x[m] =  ll_tmp_1[m] + log_sum_exp(to_vector(ll_tmp_1_2[m, ]));
        }
        // Add likelihoods
        ll_obs_tmp[j, 4] += log_sum_exp(x);
      }
      
      
      // ---- 5. obs loglik of <1,0,0> ----
      // Account for infections dates prior to baseline
      for(m in 1:M[1]){
        ll_tmp_0[m] = log_rprob_inf_pre0[i, m] + 
        intervalNormalCDF(y_star_pre0[i, m, j, 1], sigma, y_L[i, 1], y_U[i, 1]) + 
        intervalNormalCDF(y_star_pre0[i, m, j, 2], sigma, y_L[i, 2], y_U[i, 2]) + 
        intervalNormalCDF(y_star_pre0[i, m, j, 3], sigma, y_L[i, 3], y_U[i, 3]);
        
      }
      ll_obs_tmp[j, 5] += log_sum_exp(ll_tmp_0);
      
      
      // ---- 6. obs loglik of <1,0,1> ----
      for(m in 1:M[1]){
        // Account for infections dates prior to baseline
        ll_tmp_0[m] = log_rprob_inf_pre0[i, m] + 
        intervalNormalCDF(y_star_pre0[i, m, j, 1], sigma, y_L[i, 1], y_U[i, 1]) +
        intervalNormalCDF(y_star_pre0[i, m, j, 2], sigma, y_L[i, 2], y_U[i, 2]);
        
        // Account for infection between round 1 and round 2
        
        for (m2 in 1:M[3]) {
          // For round 2 account for boosting that occurred between baseline and round 1
          // as well as round 1
          real ys = log_sum_exp([y_star_pre0[i, m, j, 3], y_star_1_2_boost[i, m2, j]]);
          ll_tmp_0_2[m, m2] = log_rprob_inf_pre0[i, m] + log_rprob_inf_1_2[i, m2] +
          intervalNormalCDF(ys, sigma, y_L[i, 3], y_U[i, 3]);
        }
      }
      
      // Add likelihoods
      {
        // intermediate computation
        vector[M[1]] x;
        for (m in 1:M[1]) {
          x[m] =  ll_tmp_0[m] + log_sum_exp(to_vector(ll_tmp_0_2[m, ]));
        }
        // Add likelihoods
        ll_obs_tmp[j, 6] += log_sum_exp(x);
      }
      
      
      // ---- 7. obs loglik of <1,1,0> ----
      for(m in 1:M[1]){
        // Account for infections dates prior to baseline
        ll_tmp_0[m] = log_rprob_inf_pre0[i, m] + 
        intervalNormalCDF(y_star_pre0[i, m, j, 1], sigma, y_L[i, 1], y_U[i, 1]);
        
        // Account for infection between round 1 and round 2
        
        for (m2 in 1:M[2]) {
          // For round 2 account for boosting that occurred between baseline and round 1
          // as well as round 1
          real ys = log_sum_exp([y_star_pre0[i, m, j, 2], y_star_0_1_boost[i, m2, j, 1]]);
          real ys2 = log_sum_exp([y_star_pre0[i, m, j, 3], y_star_0_1_boost[i, m2, j, 2]]);
          
          ll_tmp_0_1[m, m2] = log_rprob_inf_pre0[i, m] + log_rprob_inf_0_1[i, m2] +
          intervalNormalCDF(ys, sigma, y_L[i, 2], y_U[i, 2]) +
          intervalNormalCDF(ys2, sigma, y_L[i, 3], y_U[i, 3]);
        }
      }
      
      {
        // intermediate computation
        vector[M[1]] x;
        for (m in 1:M[1]) {
          x[m] =  ll_tmp_0[m] + log_sum_exp(to_vector(ll_tmp_0_1[m, ]));
        }
        // Add likelihoods
        ll_obs_tmp[j, 7] += log_sum_exp(x);
      }
      
      
      // ---- 8. obs loglik of <1,1,1> ----
      for(m in 1:M[1]){
        // Account for infections dates prior to baseline
        ll_tmp_0[m] = log_rprob_inf_pre0[i, m] + 
        intervalNormalCDF(y_star_pre0[i, m, j, 1], sigma, y_L[i, 1], y_U[i, 1]);
        
        // Account for infection between round 1 and round 2
        for (m2 in 1:M[2]) {
          // For round 1 account for boosting that occurred between baseline and round 1
          real ys = log_sum_exp([y_star_pre0[i, m, j, 2], y_star_0_1_boost[i, m2, j, 1]]);
          
          ll_tmp_0_1[m, m2] = log_rprob_inf_pre0[i, m] + log_rprob_inf_0_1[i, m2] +
          intervalNormalCDF(ys, sigma, y_L[i, 2], y_U[i, 2]);
          
          for (m3 in 1:M[3]) {
            // For round 2 account for boosting that occurred between baseline and round 1
            // as well as between round 1 and round 2
            real ys2 = log_sum_exp([y_star_pre0[i, m, j, 3], y_star_0_1_boost[i, m2, j, 2],
            y_star_1_2_boost[i, m3, j]]);
            
            ll_tmp_0_1_2[m, m2, m3] = log_rprob_inf_pre0[i, m] + log_rprob_inf_0_1[i, m2] +
            log_rprob_inf_1_2[i, m3] +
            intervalNormalCDF(ys2, sigma, y_L[i, 3], y_U[i, 3]);
          }
        }
      }
      
      {
        // intermediate computation
        vector[M[1]] x;
        for (m in 1:M[1]) {
          vector[M[2]] x2;
          for (m2 in 1:M[2]) {
            x2[m2] =  ll_tmp_0_1[m, m2] + log_sum_exp(to_vector(ll_tmp_0_1_2[m, m2, ]));
          }
          x[m] = ll_tmp_0[m] + log_sum_exp(x2);
        }
        
        // Add likelihoods
        ll_obs_tmp[j, 8] += log_sum_exp(x);
      }
    }
    
    // Combine likelihoods across parameter samples
    for (k in 1:8) {
      ll_obs[i, k] = log_sum_exp(ll_obs_tmp[, k]);
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

model {
  real mu[N, 3, J_foi];
  matrix[N, 8] ll;
  
  {
    real log_mu[N, 3];
    real log1m_mu[N, 3];
    matrix[8, J_foi] ll_tmp;
    
    // Compute the log-likelihoods for alternative cases
    for (i in 1:N) {
      
      for (k in 1:J_foi) {
        for (j in 1:3){
          mu[i, j, k] = 1-exp(-(lambda[age_cat_vec[i]] * cum_foi[i, j, k] + lambda_env[age_cat_vec[i]] * dt[i, j]));
          log_mu[i, j] = log(mu[i, j, k]);
          log1m_mu[i, j] = log1m(mu[i, j, k]);
        }
        
        // Contributions of each case
        //  case <0,0,0>
        ll_tmp[1, k] = sum(log1m_mu[i, ]) + ll_obs[i, 1] + logjfoi;
        //  case <0,0,1>
        ll_tmp[2, k] = log1m_mu[i, 1] + log1m_mu[i, 2] + log_mu[i, 3] + ll_obs[i, 2] + logjfoi;
        //  case <0,1,0>
        ll_tmp[3, k] = log1m_mu[i, 1] + log_mu[i, 2] + log1m_mu[i, 3] + ll_obs[i, 3] + logjfoi;
        //  case <0,1,1>
        ll_tmp[4, k] = log1m_mu[i, 1] + log_mu[i, 2] + log_mu[i, 3] + ll_obs[i, 4] + logjfoi;
        //  case <1,0,0>
        ll_tmp[5, k] = log_mu[i, 1] + log1m_mu[i, 2] + log1m_mu[i, 3] + ll_obs[i, 5] + logjfoi;
        //  case <1,0,1>
        ll_tmp[6, k] = log_mu[i, 1] + log1m_mu[i, 2] + log_mu[i, 3] + ll_obs[i, 6] + logjfoi;
        //  case <1,1,0>
        ll_tmp[7, k] = log_mu[i, 1] + log_mu[i, 2] + log1m_mu[i, 3] + ll_obs[i, 7] + logjfoi;
        //  case <1,1,1>
        ll_tmp[8, k] = sum(log_mu[i, ]) + ll_obs[i, 8] + logjfoi;
      }
      // Marginalize out uncertainty on foi
      for (j in 1:8) {
        ll[i, j] = log_sum_exp(to_vector(ll_tmp[j, ]));
      }  
    }
    
    // Likelihood serology data
    for (i in 1:N) {
      target += log_sum_exp(ll[i,]);
    }
  }
  
  log_lambda ~ normal(mu_log_lambda, sd_log_lambda);
  log_lambda_env ~ normal(mu_log_lambda_env, sd_log_lambda_env);
  
  // Hierarchical prior on log_lambda_env
  sd_log_lambda ~ normal(0, sd_lambda);
  mu_log_lambda ~ normal(-4, 1);
  sd_log_lambda_env ~ normal(0, sd_lambda);
  mu_log_lambda_env ~ normal(-5, 1);
  
}

generated quantities {
  vector[N] case_lik;
  matrix[N, 8] ll;
  real mu_round[3];
  
  {
    real mu_ind[N, 3];
    real mu_tmp[N, 3, J_foi];
    real log_mu[N, 3];
    real log1m_mu[N, 3];
    matrix[8, J_foi] ll_tmp;
    
    // Compute the log-likelihoods for alternative cases
    for (i in 1:N) {
      
      for (k in 1:J_foi) {
        for (j in 1:3){
          mu_tmp[i, j, k] = 1-exp(-(lambda[age_cat_vec[i]] * cum_foi[i, j, k] + lambda_env[age_cat_vec[i]] * dt[i, j]));
          log_mu[i, j] = log(mu_tmp[i, j, k]);
          log1m_mu[i, j] = log1m(mu_tmp[i, j, k]);
        }
        
        // Contributions of each case
        //  case <0,0,0>
        ll_tmp[1, k] = sum(log1m_mu[i, ]) + ll_obs[i, 1] + logjfoi;
        //  case <0,0,1>
        ll_tmp[2, k] = log1m_mu[i, 1] + log1m_mu[i, 2] + log_mu[i, 3] + ll_obs[i, 2] + logjfoi;
        //  case <0,1,0>
        ll_tmp[3, k] = log1m_mu[i, 1] + log_mu[i, 2] + log1m_mu[i, 3] + ll_obs[i, 3] + logjfoi;
        //  case <0,1,1>
        ll_tmp[4, k] = log1m_mu[i, 1] + log_mu[i, 2] + log_mu[i, 3] + ll_obs[i, 4] + logjfoi;
        //  case <1,0,0>
        ll_tmp[5, k] = log_mu[i, 1] + log1m_mu[i, 2] + log1m_mu[i, 3] + ll_obs[i, 5] + logjfoi;
        //  case <1,0,1>
        ll_tmp[6, k] = log_mu[i, 1] + log1m_mu[i, 2] + log_mu[i, 3] + ll_obs[i, 6] + logjfoi;
        //  case <1,1,0>
        ll_tmp[7, k] = log_mu[i, 1] + log_mu[i, 2] + log1m_mu[i, 3] + ll_obs[i, 7] + logjfoi;
        //  case <1,1,1>
        ll_tmp[8, k] = sum(log_mu[i, ]) + ll_obs[i, 8] + logjfoi;
      }
      // Marginalize out uncertainty on foi
      for (j in 1:8) {
        ll[i, j] = log_sum_exp(to_vector(ll_tmp[j, ]));
      }  
      
      // Compute mean mu
      for (j in 1:3){
        mu_ind[i, j] = mean(mu_tmp[i, j, ]);
      }
    }
    
    for (j in 1:3) {
      mu_round[j] = mean(mu_ind[, j]);
    }
    
    for (i in 1:N) {
      case_lik[i] = log_sum_exp(ll[i, ]);
    }
  }
}
