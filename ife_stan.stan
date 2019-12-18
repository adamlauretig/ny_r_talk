data{
  int<lower = 0> N ; // number of groups obs
  int<lower = 0> J ; // number of group 2 obs
  int<lower = 0> K ; // number of latent observations
  int X[(N*J), 2]  ; // covariate matrix
  vector[(N*J)] y ; // outcome
}
transformed data {
  int<lower=1> K_choose_2;
  K_choose_2 = (K * (K - 1)) / 2;
}
parameters{
  vector[N] indiv_betas; // non-interacted coefficients
  vector[J] group_2_betas; // non-interacted coefficients
  matrix[N, K] gammas; // individual factors
  matrix[J, K] deltas; // group 2 factors
  // positive_ordered[K] gamma_sd; //gamma prior sd
  // positive_ordered[K] gamma_psi; //gamma diag
  // positive_ordered[K] delta_psi; //delta diag
  // real<lower = 0> a; // gamma hyperprior a
  // real<lower = 0> b; // gamma hyperprior b
  // real<lower = 0> c; // gamma hyperprior c
  // real<lower = 0> d; // gamma hyperprior d
  vector[K_choose_2] gamma_cov_lower_tri; // fix upper tri for parameter identification
  vector[K_choose_2] delta_cov_lower_tri; // fix upper tri for parameter identification
  
}
transformed parameters{
  vector[(N*J)] linear_predictor ;
  cholesky_factor_cov[K] gamma_L;
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
    linear_predictor[i] = indiv_betas[X[i, 1]] + group_2_betas[X[i, 2]] + (gammas[X[i, 1], ]  * deltas[ X[i, 2], ]');
  }
}
model{
  indiv_betas ~ normal(0, 5) ;
  group_2_betas ~ normal(0, 5) ;
  
  delta_cov_lower_tri ~ normal(0, 1) ;
  gamma_cov_lower_tri ~ normal(0, 1) ;

  // a ~ gamma(2, 2) ;
  // b ~ gamma(2, 2) ;
  // c ~ gamma(2, 2) ;
  // d ~ gamma(2, 2) ;
  // gamma_psi ~ gamma(a, b) ;
  // delta_psi ~ gamma(c, d) ;
  for(n in 1:N){
    gammas[n, ] ~ multi_normal_cholesky(rep_vector(0, K), gamma_L) ;
  }
  for(j in 1:J){
    deltas[j, ] ~ multi_normal_cholesky(rep_vector(0, K), delta_L) ;
  }
  y ~ normal(linear_predictor, 1) ;
}

