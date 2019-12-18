// Maxdiff model, non-centered parameterization, 
// using multivariate student-t distribution for utilities
// following https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html
// Here, different slopes for best/worst choices
data {
  int<lower=2> C; // Number of alternatives (choices) in each scenario
  int<lower=1> K; // Number of total possible choices
  int<lower=1> R; // Number of respondents
  int<lower=1> S; // Number of scenarios per respondent
  int<lower=1> YB[R, S]; // best choices
  matrix[C, K] X[R, S]; // matrix of attributes for each individual, in each set
  
  real<lower = 0> L_sigma_prior ;
  int<lower = 1> L_Omega_prior ;
  real<lower = 1> Theta_raw_prior ;
  real<lower = 1> nu ; // DoF parameter. Usually set to 3
}

parameters {
  vector[K - 1] Theta_raw[R]; // for identification
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0, upper=pi()/2>[K] L_sigma_unif; // reparameterized cauchy
  vector[K] alpha ; // for non-centered parameterization 
  real<lower=0> u; // for student-t reparameterization

}

transformed parameters {
  vector[K] Beta[R];
  vector<lower=0>[K] L_sigma;

  // matrix[K, K] L_Sigma;
  vector[C] XB[R, S];
  vector[K] Theta[R];

  for (k in 1:K) {
    L_sigma[k] = L_sigma_prior * tan(L_sigma_unif[k]) ; // reparameterized cauchy
  }
  for (r in 1:R){
    Theta[r] = append_row( Theta_raw[r], 0) ; // add 0 for identification

  }
  for(r in 1:R){ // multivariate student-t
  // via https://discourse.mc-stan.org/t/non-centered-parameterization-of-the-multivariate-student-t-distribution/4176/8
      Beta[r] = Theta[r] + sqrt(nu /u) * L_sigma .* (L_Omega * alpha) ;
  }
  for (r in 1:R) {
    for (s in 1:S) {
      XB[r,s] = X[r,s] * Beta[r];
    }
  }
}

model {
  //priors
  for (r in 1:R){
    Theta_raw[r] ~ normal(0, Theta_raw_prior) ;
  }  
  L_Omega ~ lkj_corr_cholesky(L_Omega_prior) ;
  alpha ~ std_normal() ;
  u ~ chi_square(nu) ;

  for (r in 1:R) {
    for (s in 1:S) {
      if(YB[r,s] != 9){
        YB[r,s] ~ categorical_logit(XB[r,s]) ;
      } 
    }
  }
}
// generated quantities{
//   simplex[C] choice_probs[R, S];
// 
//   for(r in 1:R){
//     for (s in 1:S) {
//       choice_probs[r, s] = softmax(XB[r,s]) ;
//     }
//   }
// }

