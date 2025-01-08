library(rstanarm)
library(loo)

set.seed(123)
Z_data = matrix(0,n,q)
for(i in 1:n){
    Z_data[i,] = Z[t_design[i],county_design[i],]
}

y_glm = y[!is.na(y)]
Z_glm = Z_data[!is.na(y),]
fit_bayes <- stan_glm(y_glm ~ Z_glm, family = binomial, prior = normal(0, 10), 
                      prior_intercept = normal(0, 10), chains = 1, iter = 2000)

log_lik <- log_lik(fit_bayes) # Calculate WAIC using the loo package 
waic_result <- waic(log_lik) # Print the WAIC result 