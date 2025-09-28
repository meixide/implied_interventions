# ---- Modularized Simulation to Understand h*, g*, and EIF ----

#set.seed(123)

n=500

pw=0.3

hid_conf=T

out = function(A,W,U) return(2 * A + W - U)

# Utility: Project a vector onto the probability simplex (used for some optimization routines)
project_onto_simplex <- function(v, z = 1) {
  if (z <= 0) stop("z must be positive")
  if (!is.numeric(v)) stop("v must be a numeric vector")
  n <- length(v)
  if (n == 1) return(z)
  if (all(v <= 0)) return(rep(0, n))
  mu <- sort(v, decreasing = TRUE)
  cumsum_mu <- cumsum(mu)
  indices <- 1:n
  rho <- max(which(mu - (cumsum_mu - z)/indices > 0))
  theta <- (cumsum_mu[rho] - z)/rho
  w <- pmax(v - theta, 0)
  w * (z / sum(w))
}

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
  W <- rbinom(n, 1, pw)
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
  W <- rbinom(n, 1, pw)
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
h0_fun <- function(w) make_h(w, h0_0, h0_1)
hstar_fun <- function(w) make_h(w, h_star_0, h_star_1)


g_fun <- g_fun_from_matrix(gmat)