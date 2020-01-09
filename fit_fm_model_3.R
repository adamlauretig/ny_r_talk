# NY R Users Meetup on January 9, 2020, presented by Adam Lauretig
library(data.table)
library(Matrix)
library(foreign)
library(rstan)
library(bayesplot)
library(ggplot2)
library(umap)
library(MASS)

# fearon and laitin data
fl_civil_war <- read.dta("repdata.dta")
m2 <- stan_model("stan_fm_3.stan")

# cleaning up strange data
fl_civil_war <-fl_civil_war[ fl_civil_war$onset < 2, ]
# only getting years with onsets (not including years where there is already ongoing war)
fl_civil_war <- fl_civil_war[ !(fl_civil_war$war == 1 & fl_civil_war$onset == 0), ]

# basic logistic regression
m1 <- glm(onset~warl+gdpenl+lmtnest+polity2l+ethfrac+relfrac, data = fl_civil_war, family = binomial())
m1_preds <- predict(m1, type = "response")

# prepping data for stan
# getting covariates
data_subset <- fl_civil_war[
  c("warl", "gdpenl", "lpopl1", "lmtnest", "polity2l", "ethfrac", "relfrac", 
    "ncontig",  "Oil", "nwstate", "instab",  "cname", "year", "onset")]
data_subset <- na.omit(data_subset)
data_subset <- data_subset[ data_subset$onset < 2, ]
# creating groups (factors for FM)
groups <- data.frame(cname = factor(data_subset$cname), 
           cnumber = as.numeric(factor(data_subset$cname)), 
           year = factor(data_subset$year), 
           cyear = as.numeric(factor(data_subset$year)))
X <- as.matrix(groups[, c(2, 4)])
# additional covariates
X2 <- as.matrix(data_subset[, 1:11])
y <- as.numeric(data_subset$onset)
data <- list(n_obs = nrow(X), X = X, X2 = X2, y = y, N = max(X[, 1]), 
  J = max(X[, 2]), K = 5,
     beta_sigma = 1,
     gamma_sigma_prior = 1, 
     delta_sigma_prior = 1, 
     gamma_omega_prior = 2, 
     delta_omega_prior = 2)

# fitting model
# m2_fit <- vb(object = m2, data = data, tol_rel_obj = 5e-3, seed = 123)
m2_fit <- sampling(object = m2, data = data, seed = 123, chains = 4, cores = 4, 
                   control = list(adapt_delta = .95))
save(m2_fit, file = "fm_model_3_fit.rdata")

# extracting predictions, comparing to logistic regression
m2_posterior <- extract(m2_fit)
hist(colMeans(m2_posterior$y_pred))
predictions <- data.frame(y = y, y_preds = colMeans(m2_posterior$y_pred), m1_preds = m1_preds)

# log loss to assess model fit
log_loss_binary <- function(actual, predicted, eps = 1e-15) {
  predicted = pmin(pmax(predicted, eps), 1-eps)
  - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
  }

log_loss_binary(actual = predictions$y, predicted = predictions$y_preds)
log_loss_binary(actual = predictions$y, predicted = predictions$m1_preds)

ggplot(data = predictions, aes(x = y, y = y_preds)) + 
  geom_jitter(width = .1, alpha = .2, size = .5) + 
  geom_jitter(data = predictions, aes(x = y, y = m1_preds), 
              width = .1, size = .5, alpha = .1, color = "red")

# now, plotting latent factors
country_means <- apply(m2_posterior$gammas, c(2:3), FUN = mean)
countries <- unique(groups[, 1:2])
countries[order(countries$cnumber),]
rownames(country_means) <- countries$cname
# sammon mapping
country_sim <- sammon(dist(country_means))
plot(country_sim$points, type = "n")
text(country_sim$points, labels = as.character(rownames(country_means)))

# using umap instead of sammon mapping
countries_umap <- umap(country_means)
plot(x = countries_umap$layout[, 1],
     y = countries_umap$layout[, 2],
     type = "n")
text(
  x = countries_umap$layout[, 1],
  y = countries_umap$layout[, 2],
  labels = as.character(rownames(country_means))
)


