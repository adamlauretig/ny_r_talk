// Bayesian Factorization Machine, N(0, 1) priors
// NY R Users Meetup on January 9, 2020, presented by Adam Lauretig
data{
  int<lower = 0> N ; // number of group 1 observations
  int<lower = 0> J ; // number of group 2 observations
  int<lower = 0> K ; // number of latent dimensions
  int X[(N*J), 2]  ; // covariate matrix
  vector [(N*J)] y ; // outcome
  real<lower = 0> beta_sigma ; // sd on regression coefficients
  real<lower = 0> y_sigma ; // sd on the outcome, y
}

parameters{
  vector[N] group_1_betas; // non-interacted coefficients
  vector[J] group_2_betas; // non-interacted coefficients
  matrix[N, K] gammas; // individual factors
  matrix[J, K] deltas; // group 2 factors
}

transformed parameters{
  real linear_predictor[(N*J)] ;
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

  // latent factors
  for(n in 1:N){
    gammas[n, ] ~ normal(rep_vector(0, K), 1) ;
  }
  
  for(j in 1:J){
    deltas[j, ] ~ normal(rep_vector(0, K), 1) ;
  }
  
  // outcome
  to_vector(y) ~ normal(linear_predictor, y_sigma) ;
}
generated quantities{
  real y_pred[(N*J)]  ;
  // real mse ;
  
  y_pred = normal_rng(linear_predictor, y_sigma) ;
  // mse = ((y - y_pred)**2)/(N*J)
}
