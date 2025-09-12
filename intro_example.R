n=500


setwd("/Users/cgmeixide/Dropbox/ivhal")

hid_conf=T

out = function(A,W,U) return(2 * A + W - U)

# Unified: True h(z=1|w)
make_h <- function(w, h0_0, h0_1) {
  ifelse(w == 0, h0_0, h0_1)
}

# Unified: True g(a=1|z,w) matrix style
make_g_matrix <- function(g00, g01, g10, g11) {
  matrix(c(g00, g01, g10, g11), nrow = 2, byrow = TRUE)
}



# Unified: g(a=1|z,w) function from matrix
g_fun_from_matrix <- function(Gmat) {
  function(z, w) Gmat[cbind(z+1, w+1)]
}

# Simulate data for one run
simulate_data <- function(n, h0_fun, g_fun) {
  W <- rbinom(n, 1, 0.3)
  U = runif(n,0,1)
  Z <- rbinom(n, 1, h0_fun(W))
  A <- as.numeric(U < g_fun(Z, W))
  Y <- rnorm(n, mean = out(A,W,hid_conf*U), sd = noise_sd)
  list(W = W, Z = Z, A = A, Y = Y)
}


# Estimate h(z|w) and g(a|z,w) via logistic regressions
estimate_h <- function(Z, W) {
  glm(Z ~ W, family = binomial)$fitted.values
}

estimate_g <- function(A, Z, W) {
  glm(A ~ Z + W, family = binomial)$fitted.values
}

# Q model (E[Y|A,Z,W])
estimate_Q <- function(Y, Z, W) {
  lm(Y ~ Z + W)
}

h_natural = function(z,w) z*make_h(w, h0_0, h0_1) + (1-z)*(1-make_h(w, h0_0, h0_1)) 
h_star = function(z,w) z*make_h(w, h_star_0, h_star_1) + (1-z)*(1-make_h(w, h_star_0, h_star_1)) 


# Calculate plugin and clever covariates
calc_plugin_eif_terms <- function(Q_model, Y, Z, W, h_star, h_natural) {
  Q_hat <- Q_model(Z,W)
  clever_term <- (h_star(Z,W)/ h_natural(Z,W)) * (Y - Q_hat)
  Q0_hat <- Q_model(rep(0,n),W)
  Q1_hat <- Q_model(rep(1,n),W)
  plugin_term <- (Q0_hat * h_star(0,W) + Q1_hat * h_star(1,W))
  list(clever_term = clever_term, plugin_term = plugin_term)
}

# Simulate intervention (counterfactual) data
simulate_intervention <- function(n, h_star_fun, g_fun) {
  W <- rbinom(n, 1, 0.3)
  Z_int <- rbinom(n, 1, h_star_fun(W))
  U = runif(n,0,1)
  A_int <- as.numeric(U < g_fun(Z_int, W))
  Y_int <- rnorm(n, mean = out(A_int,W,hid_conf*U), sd = noise_sd)
  list(Z_int = Z_int, A_int = A_int, Y_int = Y_int)
}

noise_sd=0.05
h0_0 = 0.3
h0_1 = 0.8

h_star_0 = 0.7
h_star_1 = 0.4

gmat =  make_g_matrix(0.3, 0.8, 0.7, 0.5)
g_fun <- g_fun_from_matrix(gmat)


h_star_grid_0=seq(0,1,length.out=100)
h_star_grid_1=seq(0,1,length.out=100)

implied_0=g_fun(0,0)*(1-h_star_grid_0) + g_fun(1,0)*h_star_grid_0
implied_1=g_fun(0,1)*(1-h_star_grid_1) + g_fun(1,1)*h_star_grid_1



h0_fun <- function(w) make_h(w, h0_0, h0_1)
hstar_fun <- function(w) make_h(w, h_star_0, h_star_1)



EYstar=mean(simulate_intervention(1e6, function(w) make_h(w, h_star_0 , h_star_1 ), g_fun)$Y_int)
EYstar

B=100

psinontmle=numeric(B)
coverage=numeric(B)
psi=numeric(B)
sdbound=numeric(B)
for (b in 1:B) {
  
  
  data=simulate_data(n, h0_fun, g_fun) 
  
  Y=data$Y
  A=data$A
  W=data$W
  Z=data$Z
  
  
  estimate_Q2 <- function(Y, A, W) {
    lm(Y ~ A + W)
  }
  
  summary(estimate_Q2(Y,A,W))
  aggregate(A ~ Z + W,mean,data=data)
  
  
  Qmodel=estimate_Q(Y, Z, W)
  Qn.0<- function(Z, W) coef(Qmodel)[1] + coef(Qmodel)[2]*Z + coef(Qmodel)[3]*W
  
  
  
  G_comp_z = function(Qbar,h_star_0,h_star_1) {
    Z0_df <- data.frame(Z = 0, W = W)
    Z1_df <- data.frame(Z = 1, W = W)
    Q0_hat <- Qbar(rep(0,n),W)
    Q1_hat <- Qbar(rep(1,n),W)
    hstar_fun <- function(w) make_h(w, h_star_0, h_star_1)
    return(mean(Q0_hat*(1-hstar_fun(W)) + Q1_hat*hstar_fun(W)))
  }
  
  ###
  
  psinontmle[b]=G_comp_z(Qn.0,h_star_0 ,h_star_1)
  
  
  # Initial estimate
  
  
  
  # Weights = clever covariate
  weights <- h_star(Z,W) / h_natural(Z,W)
  
  # Fluctuation model (weighted OLS)
  fluctuation_model <- lm(Y  ~offset(Qn.0(Z, W)) + 1, weights=weights)
  summary(fluctuation_model)
  epsilon <- coef(fluctuation_model)
  Qn.0old <- Qn.0  # save current version
  
  # define new Qn.0 by wrapping the fixed old version
  Qn.0 <- (function(Qold, eps) {
    function(Z, W) Qold(Z, W) + eps
  })(Qn.0old, epsilon)
  
  ###
  
  EYstar
  psi[b]=G_comp_z(Qn.0,h_star_0 ,h_star_1)
  
  #eifold=calc_plugin_eif_terms(Qn.0old, Y, Z, W, h_star, h_natural)
  
  
  #mean(eifold[[1]]) + mean(eifold[[2]])- EYstar # non zero bc of tmle
  
  #
  terms=calc_plugin_eif_terms(Qn.0, Y, Z, W,  h_star, h_natural)
  #mean(terms[[1]]) + mean(terms[[2]])- EYstar
  
  eif=terms[[1]] + terms[[2]]- psi[b]
  
  #tmle makes the first term zero
  # in the old eif, we needed to know EYstar. in the new one, we do not
  # both eif have mean zero though
  
  sdbound[b]=sqrt(var(eif)/n)
  
  
  
  
}


# Combine the vectors into a data frame
results <- data.frame(
  psinontmle = psinontmle,
  psi = psi,
  sdbound = sdbound
)

# Create directory if it doesn't exist
dir.create("sims/part1", recursive = TRUE, showWarnings = FALSE)

# Define the file name with the value of n
filename <- paste0("sims/part1/", n, ".txt")

# Write the data frame to the file with headers
write.table(results, file = filename, row.names = FALSE, col.names = TRUE, sep = "\t")


