// factorization machine with hierarchical priors on latent factors
data{
  int<lower = 0> N ; // number of group 1 observations
  int<lower = 0> J ; // number of group 2 observations
  int<lower = 0> K ; // number of latent dimensions
  int X[(N*J), 2]  ; // covariate matrix
  vector[(N*J)] y ; // outcome
  real<lower = 0> beta_sigma ; // sd on regression coefficients
  real<lower = 0> y_sigma ; // sd on the outcome, y
  real<lower = 0> gamma_sigma_prior ; //sd on gamma factors
  real<lower = 0> delta_sigma_prior ; //sd on gamma factors

}

parameters{
  vector[N] group_1_betas; // non-interacted coefficients
  vector[J] group_2_betas; // non-interacted coefficients
  // matrix[N, K] gammas; // individual factors
  // matrix[J, K] deltas; // group 2 factors
  matrix[N, K] gamma_mu ; //gamma prior mean
  vector<lower=0, upper=pi()/2>[K] gamma_sigma_unif ; // reparameterized Cauchy
  cholesky_factor_corr[K] gamma_omega; // correlation matrix
  matrix[K, N] gamma_a ; // for non-centered parameterization
  matrix[J , K] delta_mu ; //delta prior mean
  vector<lower=0, upper=pi()/2>[K] delta_sigma_unif ; // reparameterized Cauchy
  cholesky_factor_corr[K] delta_omega; // correlation matrix
  matrix[K, J] delta_a ; // for non-centered parameterization
}

transformed parameters{
  real linear_predictor[(N*J)] ;
  vector<lower=0>[K] gamma_sigma ;
  matrix[N, K] gammas ;
  vector<lower=0>[K] delta_sigma;
  matrix[J, K] deltas ;
  for (k in 1:K) {
    gamma_sigma[k] = gamma_sigma_prior .* tan(gamma_sigma_unif[k]) ; // reparameterized cauchy
  }
  for (k in 1:K) {
    delta_sigma[k] = delta_sigma_prior .* tan(delta_sigma_unif[k]) ; // reparameterized cauchy
  }
  gammas = gamma_mu + (diag_pre_multiply(gamma_sigma, gamma_omega) * gamma_a)' ; //TODO: fix dimensions here
  deltas = delta_mu + (diag_pre_multiply(delta_sigma, delta_omega) * delta_a)' ;

  for(i in 1:(N*J)){
    linear_predictor[i] = 
      group_1_betas[X[i, 1]] + group_2_betas[X[i, 2]] + 
      (gammas[X[i, 1], ]  * deltas[ X[i, 2], ]');
  }
}

model{
  
  // regression coefficients
  group_1_betas ~ normal(0, beta_sigma) ;
  group_2_betas ~ normal(0, beta_sigma) ;

  gamma_omega ~ lkj_corr_cholesky(5) ;
  delta_omega ~ lkj_corr_cholesky(5) ;
  // gamma_sigma_unif ~ uniform(0, pi()/2) ;
  // delta_sigma_unif ~ uniform(0, pi()/2) ;
  to_vector(gamma_a) ~ std_normal() ;
  to_vector(delta_a) ~ std_normal() ;
  
  for(n in 1:N){
    gamma_mu[n, ] ~ normal(0, 1) ;
  }  
  for(j in 1:J){
    delta_mu[j, ] ~ normal(0, 1) ;
  }
  
  
  // outcome
  y ~ normal(linear_predictor, y_sigma) ;
}
generated quantities{
  real y_pred[(N*J)]  ;
  y_pred = normal_rng(linear_predictor, y_sigma) ;
}
