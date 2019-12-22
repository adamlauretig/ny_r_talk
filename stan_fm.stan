data{
  int<lower = 0> N ; // number of group 1 observations
  int<lower = 0> J ; // number of group 2 observations
  int<lower = 0> K ; // number of latent dimensions
  int X[(N*J), 2]  ; // covariate matrix
  vector[(N*J)] y ; // outcome
  real<lower = 0> beta_sigma ; // sd on regression coefficients
  real<lower = 0> y_sigma ; // sd on the outcome, y
  real<lower = 0> a_hyperprior_1 ; //ARD hyperprior
  real<lower = 0> a_hyperprior_2 ; //ARD hyperprior
  real<lower = 0> b_hyperprior_1 ; //ARD hyperprior
  real<lower = 0> b_hyperprior_2 ; //ARD hyperprior
}

parameters{
  vector[N] group_1_betas; // non-interacted coefficients
  vector[J] group_2_betas; // non-interacted coefficients
  matrix[N, K] gammas; // individual factors
  matrix[J, K] deltas; // group 2 factors
  positive_ordered[K] gamma_sd; //gamma prior sd
  real<lower = 0> a; // gamma hyperprior a
  real<lower = 0> b; // gamma hyperprior b
}

transformed parameters{
  vector[(N*J)] linear_predictor ;
  for(i in 1:(N*J)){
    linear_predictor[i] = group_1_betas[X[i, 1]] + group_2_betas[X[i, 2]] + (gammas[X[i, 1], ]  * deltas[ X[i, 2], ]');
  }
}

model{
  
  // regression coefficients
  group_1_betas ~ normal(0, beta_sigma) ;
  group_2_betas ~ normal(0, beta_sigma) ;
  
  // ARD prior
  a ~ gamma(a_hyperprior_1, a_hyperprior_2) ;
  b ~ gamma(b_hyperprior_1, b_hyperprior_2) ;
  gamma_sd ~ gamma( a, b) ;

  // latent factors
  for(n in 1:N){
    gammas[n, ] ~ normal(rep_vector(0, K), gamma_sd) ;
  }
  
    for(j in 1:J){
    deltas[j, ] ~ normal(rep_vector(0, K), 1) ;
  }
  
  // outcome
  y ~ normal(linear_predictor, y_sigma) ;
}

