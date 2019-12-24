data{
  int<lower = 0> N ; // number of groups obs
  int<lower = 0> J ; // number of group 2 obs
  int<lower = 0> K ; // number of latent observations
  int X[(N*J), 2]  ; // covariate matrix
  vector[(N*J)] y ; // outcome
  real<lower = 0> beta_sigma ; // sd on regression coefficients
  real<lower = 0> y_sigma ; // sd on the outcome, y
  real<lower = 0> a_hyperprior_1 ; //ARD hyperprior
  real<lower = 0> a_hyperprior_2 ; //ARD hyperprior
  real<lower = 0> b_hyperprior_1 ; //ARD hyperprior
  real<lower = 0> b_hyperprior_2 ; //ARD hyperprior
}
transformed data {
  int<lower=1> K_choose_2;
  K_choose_2 = (K * (K - 1)) / 2;
}
parameters{
  vector[N] group_1_betas; // non-interacted coefficients
  vector[J] group_2_betas; // non-interacted coefficients
  matrix[N, K] gammas; // individual factors
  matrix[J, K] deltas; // group 2 factors
  // vector<lower = 0>[K] gamma_psi; //gamma diag
  // real<lower = 0> a; // gamma hyperprior a
  // real<lower = 0> b; // gamma hyperprior b
  vector[K_choose_2] gamma_cov_lower_tri; // fix upper tri for parameter identification
  // vector[K] gamma_mean[N]; // for non-centered parameterization
  // vector[K] z[N]; // for non-centered parameterization
  vector[K_choose_2] delta_cov_lower_tri; // fix upper tri for parameter identification
  
}
transformed parameters{
  vector[(N*J)] linear_predictor ;
  cholesky_factor_cov[K] gamma_L;
  // matrix[N, K] gammas; // individual factors
  cholesky_factor_cov[K] delta_L;
  for (k in 1:K)
    gamma_L[k, k] = 1;
  {
    int i = 1;
    for (m in 2:K) {
      for (n in 1:(m - 1)) {
        gamma_L[m, n] = gamma_cov_lower_tri[i];
        gamma_L[n, m] = 0;
        i += 1;
      }
    }
  }
  
  for (k in 1:K)
    delta_L[k, k] = 1;
  {
    int i = 1;
    for (m in 2:K) {
      for (n in 1:(m - 1)) {
        delta_L[m, n] = delta_cov_lower_tri[i];
        delta_L[n, m] = 0;
        i += 1;
      }
    }
  }

  

  
  for(i in 1:(N*J)){
    linear_predictor[i] = group_1_betas[X[i, 1]] + group_2_betas[X[i, 2]] + 
    (gammas[X[i, 1], ]  * deltas[ X[i, 2], ]');
  }
}
model{
  group_1_betas ~ normal(0, beta_sigma) ;
  group_2_betas ~ normal(0, beta_sigma) ;
  // for(n in 1:N){
  //   z[n] ~ normal(0, .1) ;
  //   gamma_mean[n] ~ normal(0, .1) ;
  // }
  delta_cov_lower_tri ~ normal(0, 1) ;
  gamma_cov_lower_tri ~ normal(0, 1) ;

  // a ~ gamma(a_hyperprior_1, a_hyperprior_2) ;
  // b ~ gamma(b_hyperprior_1, b_hyperprior_2) ;
  // gamma_psi ~ gamma(a, b) ;
  for(n in 1:N){
    gammas[n, ] ~  multi_normal_cholesky(rep_vector(0, K), gamma_L) ;
  }
  for(j in 1:J){
    deltas[j, ] ~ multi_normal_cholesky(rep_vector(0, K), gamma_L) ;
  }
  y ~ normal(linear_predictor, y_sigma) ;
}
generated quantities{
  real y_pred[(N*J)]  ;
  y_pred = normal_rng(linear_predictor, y_sigma) ;
}
