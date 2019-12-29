rm(list = ls())
library(data.table)
library(Matrix)
library(foreign)
library(rstan)
fl_civil_war <- read.dta("repdata.dta")
m2 <- stan_model("stan_fm_3.stan")

# /* Model #1 */
# logit onset warl gdpenl lpopl lmtnest ncontig Oil nwstate instab polity2l ethfrac relfrac ,nolog
# 
m1 <- glm(onset~warl+gdpenl+lpop+lmtnest+ncontig+Oil+nwstate+instab+polity2l+ethfrac+relfrac, data = fl_civil_war)
m1_preds <- predict(m1, type = "response")

data_subset <- fl_civil_war[
  c("warl", "gdpenl", "lmtnest", "polity2l", "ethfrac", "relfrac", "cname", "year", "onset")]
data_subset <- na.omit(data_subset)
groups <- data.frame(cname = factor(data_subset$cname), 
           cnumber = as.numeric(factor(data_subset$cname)), 
           year = factor(data_subset$year), 
           cyear = as.numeric(factor(data_subset$year)))
X <- as.matrix(groups[, c(2, 4)])
X2 <- as.matrix(data_subset[, 1:6])
y <- as.numeric(data_subset$onset)
data <- list(X = X, X2 = X2, y = y, N = max(X[, 1]), J = max(X[, 2]), K = 5,
     beta_sigma = 2,
     gamma_sigma_prior = 1, 
     delta_sigma_prior = 1, 
     gamma_omega_prior = 5, 
     delta_omega_prior = 5)
# todo: figure out how to deal w/imbalananced panel

m2_fit <- vb(object = m2, data = data)
