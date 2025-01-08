set.seed(1)
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)

design_group = as.factor(h_design)
age_design = as.vector(DATA$age)

code <- nimbleCode({
  for(i in 1:n){
      y[i] ~ dbern(prob[i])
      prob[i] <- ilogit( beta_h[h_design[i]] + age_design[i]*beta_age[1] + 
                        inprod(Z[t_design[i],county_design[i],1:q], alpha[1:q]) +
                     inprod(MI[t_design[i],county_design[i],1:L], eta[t_design[i], 1:L])+
                       delta[county_design[i]])
  }
    
    
  # Age AR
  for (h in 1:H) {
    for(i_J in 1:J)  xi[h, 1, i_J] ~ dnorm(mu_xi[i_J], var = kappasq_1)
    for (t in 2:T){ for(i_J in 1:J){
          xi[h, t, i_J] ~ dnorm( inprod(B[i_J, 1:J], xi[h, t - 1, 1:J]), var = sigmasq_xi)
    }}
  }
  for(i_J in 1:J) mu_xi[i_J] ~ dnorm(0, var = kappasq_0)  
  for(i_J in 1:J){
      B[i_J, i_J] ~ dnorm(1, var = kappasq_B)
      for(ii_J in (seq(from = 1, to = J)[-i_J])){  B[i_J, ii_J] ~ dnorm(0, var = kappasq_B) }
  }
  sigmasq_xi ~ dinvgamma(a_xi, b_xi)


  # Spatial AR
  for (t in 2:T){for(i_L in 1:L){
    eta[t, i_L] ~ dnorm( inprod(A[i_L, 1:L], eta[t - 1, 1:L]), var = sigmasq_eta)
  }}
  for(i_L in 1:L) eta[1, i_L] ~ dnorm(0, var = kappasq_2)
  for(i_L in 1:L){
      A[i_L, i_L] ~ dnorm(1, var = kappasq_A)
      for(ii_L in (seq(from = 1, to = L)[-i_L])){  A[i_L, ii_L] ~ dnorm(0, var = kappasq_A) }
  }
  sigmasq_eta ~ dinvgamma(a_eta, b_eta)

  # beta
  for (i_h in 1:H) {
    beta_h[i_h] ~ dnorm(0, var = kappasq_beta_h)
  }
  beta_age[1] ~ dnorm(0, var = kappasq_beta_age)

  # alpha S-and-S
  for (i_q in 1:q) {
    include[i_q] ~ dbern(0.5)  # Spike-and-slab prior
    alpha_latent[i_q] ~ dnorm(0, var = kappasq_alpha)
    alpha[i_q] <- alpha_latent[i_q] * include[i_q]
  }
    
  # prior for delta
  for(i in 1:M){
      delta[i] ~ dnorm(mu_delta, var = sigmasq_delta)
  }
  mu_delta ~ dnorm(0, var = 100)
  sigmasq_delta ~ dinvgamma(1, 1)
    
  
})

kappasq_beta_h = 100
kappasq_beta_age = 100
kappasq_alpha = 100
kappasq_0 = 100
kappasq_1 = 100
kappasq_B = 100
a_xi = 1
b_xi = 1
kappasq_2 = 100
kappasq_A = 100
a_eta = 1
b_eta = 1

constants <- list(
    n = n,
    T = T,
    M = M,
    q = q,
    J = J,
    L = L,
    H = H,
    h_design = h_design,
    age_design = age_design,
    t_design = t_design,
    county_design = county_design,
    Z = Z,
    MI = MI,
    kappasq_beta_h = kappasq_beta_h, kappasq_beta_age = kappasq_beta_age, kappasq_alpha = kappasq_alpha, 
    kappasq_0 = kappasq_0, kappasq_1 = kappasq_1, kappasq_B = kappasq_B, a_xi = a_xi, b_xi = b_xi,   
    kappasq_2 = kappasq_2, kappasq_A = kappasq_A, a_eta = a_eta, b_eta = b_eta
)   

inits <- list(
    beta_h = rep(0,H), beta_age = c(0), alpha_latent = rep(0, q), include = rep(1,q), 
    xi = array(0, dim = c(H, T, J)), mu_xi = rep(0, J), sigmasq_xi = 1, B = diag(J),  
    eta = matrix(0, T, L), sigmasq_eta = 1,  A = diag(L),  
    mu_delta = 0, sigmasq_delta = 1, delta = rep(0, M)
)
data = list(
    y = y
)
model <- nimbleModel(code, constants = constants, data = data, inits = inits)


initInfo <- model$initializeInfo()
uninitializedNodes <- initInfo$uninitializedNodes
warnings <- initInfo$warnings
errors <- initInfo$errors
cat("\nWarnings:\n")
print(warnings)
cat("\nErrors:\n")
print(errors)



# Set up the MCMC configuration
mcmcConf <- configureMCMC(model)
mcmcConf$setMonitors('beta_h', 'beta_age', 'alpha_latent', 'alpha', 'xi', 'eta','include', 'delta', 'mu_delta', 'sigmasq_delta')
# Build and compile the MCMC
mcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(mcmc, project = model)

n_batches <- 200
batch_size <- 200
n_batches_save = 40
num_thin <- 20
mcmc.out = runMCMC(Cmcmc, niter=batch_size*n_batches, 
        nburnin=batch_size*(n_batches - n_batches_save),
        thin=num_thin, nchains=1, setSeed = TRUE, progressBar = TRUE)

samples_save = as.matrix(mcmc.out)
# Calculate WAIC
waic_results = calculateWAIC(Cmcmc)

## Prediction
## -----------## -----------## -----------## -----------## -----------## -----------
num_save = batch_size * n_batches_save / num_thin
age_range = c(10,100)

# Designs for prediction
age_grid = seq(from = 20, to = 90, by = 5)
age_grid = (age_grid - age_range[1]) / (diff(age_range))
num_G = length(age_grid)
phi_grid = Bspline.Basis(age_grid, r = J, knots = internal_nodes, degree = degree)

# Posterior prediction
posterior_params = as.data.frame(samples_save)

# Extracting beta vectors
beta_h <- as.matrix(posterior_params[,grep("^beta_h", colnames(posterior_params))])
beta_age <- as.matrix(posterior_params[,grep("^beta_age", colnames(posterior_params))])

# Extracting alpha vector
alpha <- as.matrix(posterior_params[,grep("^alpha", colnames(posterior_params))])

# # Extracting 3-d array xi
# xi <- array(0, dim = c(num_save, H, T, J))
# for(h in 1:H) {
#   for(t in 1:T) {
#     for(j in 1:J) {
#       xi[, h, t, j] <- as.matrix(posterior_params[,paste0("xi[", h,", ", t, ", ", j, "]")])
#     }
#   }
# }

# Extracting matrix eta
eta <- array(0, dim = c(num_save, T, L))
for(t in 1:T) {
  for(l in 1:L) {
    eta[, t, l] <- as.matrix(posterior_params[,paste0("eta[", t, ", ", l, "]")])
  }
}

# Extracting delta
delta <- as.matrix(posterior_params[,grep("^delta", colnames(posterior_params))])

# county-level prediction
prob_summary_M_H_T_Age_MCMC = array(0, dim = c(num_save, M,H,T,num_G))
prob_summary_M_H_T_Age = array(0, dim = c(M,H,T,num_G))
prob_summary_M_H_T_Age_upper = array(0, dim = c(M,H,T,num_G))
prob_summary_M_H_T_Age_lower = array(0, dim = c(M,H,T,num_G))

# Make prediction
for(i in 1:M){
    temp_i = delta[, i, drop = FALSE]
    for (t in 1:T){
        temp_it =  temp_i + alpha[, 1:q]%*% as.vector(Z[t,i,1:q]) + eta[, t, 1:L]%*%as.vector(MI[t,i,1:L]) 
        for(h in 1:H){
            temp_ith = temp_it + beta_h[, h]
            for(i_age in 1:num_G){
            prob_summary_M_H_T_Age_MCMC[,i,h,t,i_age] = as.vector(temp_ith + beta_age[,1]*i_age )
        }}
    }
}
prob_summary_M_H_T_Age_MCMC = ilogit( prob_summary_M_H_T_Age_MCMC )  
prob_summary_M_H_T_Age <- apply(prob_summary_M_H_T_Age_MCMC, c(2, 3, 4, 5), mean)
prob_summary_M_H_T_Age_upper <- apply(prob_summary_M_H_T_Age_MCMC, c(2, 3, 4, 5), quantile, probs = c(0.95))
prob_summary_M_H_T_Age_lower <- apply(prob_summary_M_H_T_Age_MCMC, c(2, 3, 4, 5), quantile, probs = c(0.05))
prob_summary_M_H_T_Age_average_mcmc <- apply(prob_summary_M_H_T_Age_MCMC, 1, mean)

# Make prediction to all individuals
Phi_Data = Bspline.Basis(DATA$age, r = J, knots = internal_nodes, degree = degree) 
prob_summary_individuals_MCMC = matrix(0, num_save, n)
prob_summary_individuals = rep(0, n)
prob_summary_individuals_upper = rep(0, n)
prob_summary_individuals_lower = rep(0, n)
for(i_individual in 1:n){
    i = county_design[i_individual]
    t = t_design[i_individual]
    h = h_design[i_individual]
    prob_summary_individuals_MCMC[,i_individual] = delta[, i, drop = FALSE] + 
                alpha[, 1:q]%*% as.vector(Z[t,i,1:q]) + eta[, t, 1:L]%*%as.vector(MI[t,i,1:L]) +
                                               beta_h[, h] + beta_age[,1]*age_design[i_individual]
}

prob_summary_individuals_MCMC = ilogit( prob_summary_individuals_MCMC )
prob_summary_individuals <- apply(prob_summary_individuals_MCMC, 2, mean)
prob_summary_individuals_upper <- apply(prob_summary_individuals_MCMC, 2, quantile, probs = c(0.95))
prob_summary_individuals_lower <- apply(prob_summary_individuals_MCMC, 2, quantile, probs = c(0.05))
prob_summary_individuals_average_mcmc = apply(prob_summary_individuals_MCMC, 1, mean)
