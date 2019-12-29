// factorization machine with hierarchical priors on latent factors
data{
  int<lower = 0> n_obs ; // number of group 1 observations
  int<lower = 0> N ; // number of group 1 observations
  int<lower = 0> J ; // number of group 2 observations
  int<lower = 0> K ; // number of latent dimensions
  int X[(n_obs), 2]  ; // covariate matrix
  matrix [(n_obs), 6] X2 ; //
  int<lower = 0, upper = 1> y[(n_obs)] ; // outcome
  real<lower = 0> beta_sigma ; // sd on regression coefficients
  real<lower = 0> gamma_sigma_prior ; //sd on gamma factors
  real<lower = 0> delta_sigma_prior ; //sd on gamma factors
  real<lower = 0> gamma_omega_prior ; // prior on omega value for gamma
  real<lower = 0> delta_omega_prior ; // prior on omega value for gamma

}

parameters{
  vector[N] group_1_betas; // non-interacted coefficients
  vector[J] group_2_betas; // non-interacted coefficients
  
  vector[6] coef_betas ; // betas on additional coefficients
  
  matrix[N, K] gamma_mu ; //gamma prior mean
  vector<lower=0>[K] gamma_sigma ; // embedding SD
  cholesky_factor_corr[K] gamma_omega; // correlation matrix
  matrix[K, N] gamma_a ; // for non-centered parameterization
  
  matrix[J , K] delta_mu ; //delta prior mean
  vector<lower=0>[K] delta_sigma ; // embedding SD
  cholesky_factor_corr[K] delta_omega; // correlation matrix
  matrix[K, J] delta_a ; // for non-centered parameterization
}

transformed parameters{
  real linear_predictor[(n_obs)] ;
  matrix[N, K] gammas ;
  matrix[J, K] deltas ;
  vector[n_obs] coef_lp ;
  
  coef_lp = X2 * coef_betas ;
    
  gammas = gamma_mu + (diag_pre_multiply(gamma_sigma, gamma_omega) * gamma_a)' ;
  deltas = delta_mu + (diag_pre_multiply(delta_sigma, delta_omega) * delta_a)' ;

  for(i in 1:(n_obs)){
    linear_predictor[i] = 
      group_1_betas[X[i, 1]] + group_2_betas[X[i, 2]] + coef_lp[i] +
      (gammas[X[i, 1], ]  * deltas[ X[i, 2], ]');
  }

  
}

model{
  
  // regression coefficients
  group_1_betas ~ normal(0, beta_sigma) ;
  group_2_betas ~ normal(0, beta_sigma) ;
  coef_betas ~ normal(0, beta_sigma) ;


  gamma_sigma ~ normal(0, gamma_sigma_prior) ; 
  delta_sigma ~ normal(0, delta_sigma_prior) ; 

  // correlation matrices
  gamma_omega ~ lkj_corr_cholesky(gamma_omega_prior) ;
  delta_omega ~ lkj_corr_cholesky(delta_omega_prior) ;
  
  // for non-centered parameterization
  to_vector(gamma_a) ~ std_normal() ;
  to_vector(delta_a) ~ std_normal() ;
  
  
  // hierarchical means
  for(n in 1:N){
    gamma_mu[n, ] ~ normal(0, 1) ;
  }  
  for(j in 1:J){
    delta_mu[j, ] ~ normal(0, 1) ;
  }
  
  // outcome
  y ~ bernoulli_logit(linear_predictor) ;
}
generated quantities{
  real y_pred[(n_obs)]  ;
  y_pred = inv_logit(linear_predictor) ;
}
