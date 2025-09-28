# =============================================
# EM-HAL for KL Projection of g* onto {g(h)}
# =============================================

set.seed(1)

# --- Libraries
library(hal9001)
library(glmnet)

# --- Switches
use_simulated_data <- FALSE   # TRUE = synthetic; FALSE = real (needs data.R / medicare_loader)
max_iters          <- 10
tol                <- 1e-4    # convergence tolerance on log-likelihood

# --- Desired intervention g*(A|W): only parameter(s) you set
mustar <- 1      # if A is continuous, mean for N(mustar, sdstar^2)
sdstar <- 0.5

# =============================================
# Helpers for HAL design on W
# =============================================
make_hal_design <- function(W) {
  basis_list <- enumerate_basis(
    x = W,
    max_degree = NULL,
    smoothness_orders = rep(0, ncol(W)),
    include_zero_order = TRUE,
    include_lower_order = TRUE
  )
  H <- make_design_matrix(W, basis_list)
  list(H = as.matrix(H), basis = basis_list)
}

make_newx_from_basis <- function(basis_list, W_new) {
  make_design_matrix(as.matrix(W_new), basis_list)
}

# =============================================
# Load data & define \hat p(A|Z,W)
# =============================================

# --- Synthetic data
n <- 1000
p <- 2
gamma <- 2
sdnoise <- 0.2

# Covariates
W <- matrix(runif(n * p, -2, 2), nrow = n, ncol = p)

# Instrument Z ~ Bernoulli(plogis(ezw(W)))
ezw <- function(W) W[,1] * sqrt(abs(W[,2])) * sign(W[,2])
Z <- rbinom(n, 1, plogis(ezw(W)))

  # Continuous A | Z,W ~ N(mu(Z,W), sdnoise^2)
  eazw <- function(z, W) gamma * z + sin(W[,1]) * log(1 + W[,2]^2)
  mu <- eazw(Z, W)
  A <- rnorm(n, mean = mu, sd = sdnoise)
  
  pawz <- function(a, w, z) {           # \hat p(A=a|Z=z,W=w)
    m <- gamma * z + sin(w[1]) * log(1 + w[2]^2)
    dnorm(a, mean = m, sd = sdnoise)
  }
  # Desired g*: N(mustar, sdstar^2) independent of W
  Astar <- rnorm(n, mean = mustar, sd = sdstar)
  gstar <- function(a, w) dnorm(a, mean = mustar, sd = sdstar)



# =============================================
# Build HAL design for h(Z=1|W) (logistic link)
# =============================================
hal_W <- make_hal_design(W)
X <- hal_W$H
basis_W <- hal_W$basis

# =============================================
# EM Initialization: h^{(0)}(1|W) (start at 0.5)
# =============================================
h_curr <- rep(0.5, n)

# Utility: observed-data log-likelihood for current h
loglik_obs <- function(h_vec) {
  # \ell(h) = sum_i log( p(A_i^*|0,W_i)*(1-h_i) + p(A_i^*|1,W_i)*h_i )
  denom <- mapply(function(ai, wi, h) {
    p0 <- pawz(ai, wi, z = 0)
    p1 <- pawz(ai, wi, z = 1)
    p0 * (1 - h) + p1 * h
  }, Astar, split(W, row(W)), h_vec)
  sum(log(pmax(denom, .Machine$double.eps)))
}

# Duplication trick to fit the M-step with glmnet
duplication_for_glmnet <- function(W, tau, basis_W) {
  n <- nrow(W)
  X <- make_newx_from_basis(basis_W, W)
  X2 <- rbind(X, X)
  y2 <- c(rep(1L, n), rep(0L, n))
  w2 <- c(tau, 1 - tau)
  list(X = X2, y = y2, w = w2)
}

# =============================================
# EM Loop
# =============================================
ll_prev <- -Inf
for (k in 0:max_iters) {
  print(k)
  # --- E-step: tau_i = P_k(Z=1 | A_i^*, W_i)
  tau <- mapply(function(ai, wi, h) {
    num <- pawz(ai, wi, z = 1) * h
    den <- pawz(ai, wi, z = 0) * (1 - h) + pawz(ai, wi, z = 1) * h
    num / pmax(den, .Machine$double.eps)
  }, Astar, split(W, row(W)), h_curr)
  
  # --- M-step: maximize Q via HAL/logistic with duplication trick
  dup <- duplication_for_glmnet(W, tau, basis_W)
  cv_fit <- cv.glmnet(dup$X, dup$y, family = "binomial",
                      weights = dup$w, intercept = FALSE, standardize = FALSE)
  fit <- glmnet(dup$X, dup$y, family = "binomial",
                lambda = cv_fit$lambda.1se, weights = dup$w,
                intercept = FALSE, standardize = FALSE)
  
  # Update h^{(k+1)}(1|W_i)
  XW <- make_newx_from_basis(basis_W, W)
  h_next <- as.numeric(predict(fit, newx = XW, type = "response"))
  
  # Check convergence using observed-data log-likelihood
  ll_curr <- loglik_obs(h_next)
  if (k > 0 && abs(ll_curr - ll_prev) < tol) {
    message(sprintf("Converged at iteration %d: Δll=%.3e", k, ll_curr - ll_prev))
    h_curr <- h_next
    break
  }
  
  h_curr <- h_next
  ll_prev <- ll_curr
  if (k == max_iters) {
    message("Reached max_iters without hitting tol; returning last iterate.")
  }
}

# =============================================
# Outputs
# =============================================
# h_hat(W) = P(Z=1 | W) under the KL-projection g(h^\dagger)
h_hat <- h_curr

# A compact list you can save/return:
em_hal_result <- list(
  h_hat = h_hat,
  basis_W = basis_W,
  # helper to predict h_hat at new W:
  predict_h = function(Wnew) {
    Xnew <- make_newx_from_basis(basis_W, as.matrix(Wnew))
    as.numeric(predict(fit, newx = Xnew, type = "response"))
  },
  loglik = ll_prev,
  use_simulated_data = use_simulated_data
)


# =============================================
# Diagnostics & Visualization
# =============================================


  # Simulate A_tilde from implied mixture g(h)(A|W) = sum_z p(A|z,W) h_hat(z|W)
  # For binary A: sample Z~Bern(h_hat), then A~Bern(p1(z,W))
  # For continuous A: sample Z~Bern(h_hat), then A~N(mu(z,W), sdnoise^2)
  Z_tilde <- rbinom(n, 1, em_hal_result$h_hat)
 
    # Continuous case: need mu(z,W) and sdnoise from simulation setup
    mu_tilde <- gamma * Z_tilde + sin(W[,1]) * log(1 + W[,2]^2)
    A_tilde <- rnorm(n, mean = mu_tilde, sd = sdnoise)
    
    # Density overlays
    rng <- range(A, Astar, A_tilde)
    
    pdf("/Users/cgmeixide/Dropbox/Aplicaciones/Overleaf/causal_inference_implied_interventions/figures/dens1.pdf", width = 3, height = 5.5, useDingbats = FALSE)
    
    plot(density(A), xlim = rng, ylim = c(0, max(density(A)$y, density(Astar)$y, density(A_tilde)$y)),
         main = "mu*=1", xlab = "A")
    lines(density(Astar), lty = 2)
    lines(density(A_tilde),col='red')
    #legend("topright", legend = c("Empirical A", "Target A*", "Projection A*"),
    #       lty = c(1, 2, 1), col = c("black", "black", "red"), bty = "n")
    dev.off()


 
    library(plotly)
    W1_seq <- seq(min(W[,1]), max(W[,1]), length.out = 80)
    W2_seq <- seq(min(W[,2]), max(W[,2]), length.out = 80)
    grd <- expand.grid(W1 = W1_seq, W2 = W2_seq)
    h_grid <- em_hal_result$predict_h(as.matrix(grd))
    H <- matrix(h_grid, nrow = length(W1_seq), byrow = FALSE)
    plot_ly(x = W1_seq, y = W2_seq, z = ~H) |>
      add_surface(showscale = FALSE) |>
      layout(scene = list(
        xaxis = list(title = "W1"),
        yaxis = list(title = "W2"),
        zaxis = list(title = "h*"),
        aspectratio = list(x = 1, y = 1, z = 0.6)
      ))
  





