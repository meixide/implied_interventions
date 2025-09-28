set.seed(123)

setwd("/Users/cgmeixide/Projects/implied_interventions")
source("data.R")
data <- medicare_loader(4)

Wall <- data[, c("W.zip_msa_list", "W.birthyear_list", "W.female_list", "W.dep_dx_0m")]
W    <- data$W.dep_dx_0m
Z    <- data$treatment
A    <- data$approved_app
Yprim <- data[, 1]
Y    <- ifelse(Yprim == 0, 0, 1)

## --- HAL helpers ------------------------------------------------------------
# Minimal reimplementation to build design matrices for new data given basis list
make_newx_from_basis <- function(basis_list, X_new) {
  X_new <- as.matrix(X_new)
  n <- nrow(X_new)
  B <- matrix(NA_real_, nrow = n, ncol = length(basis_list))
  for (j in seq_along(basis_list)) {
    info <- basis_list[[j]]
    ind <- rep(1, n)
    for (k in seq_along(info$cols)) {
      ind <- ind * (X_new[, info$cols[k]] >= info$cutoffs[k])
    }
    B[, j] <- ind
  }
  B
}

## --- Prob helpers -----------------------------------------------------------
make_h <- function(w, h0_0, h0_1) ifelse(w == 0, h0_0, h0_1)
hfun   <- function(z, f) z * f + (1 - z) * (1 - f)  # P(Z=z | W) in 1-dim form

## --- Libraries --------------------------------------------------------------
library(hal9001)
library(glmnet)

## --- Bases & fits -----------------------------------------------------------

Win <- as.matrix(Wall)
n   <- nrow(Win)

# h(z|w)
hal_w <- enumerate_basis(x = Win, max_degree = 2)
HW    <- make_design_matrix(Win, hal_w, p_reserve = 1)
fitcvw <- cv.glmnet(HW, Z, family = "binomial", intercept = FALSE, standardize = FALSE)
fitw   <- glmnet(HW, Z, family = "binomial", lambda = fitcvw$lambda.min,
                 intercept = FALSE, standardize = FALSE)
hw <- as.numeric(predict(fitw, newx = HW, type = "response"))

# Q(z,w) = E[Y|Z,W]
hal_zw <- enumerate_basis(x = cbind(Z, Win), max_degree = 2)
HZW    <- make_design_matrix(cbind(Z, Win), hal_zw, p_reserve = 1)
fitcvy <- cv.glmnet(HZW, Y, family = "binomial", intercept = FALSE, standardize = FALSE)
fity   <- glmnet(HZW, Y, family = "binomial", lambda = fitcvy$lambda.min,
                 intercept = FALSE, standardize = FALSE)

# Wrapper producing Q(Z,W) on demand (compatible with offset fluctuation)
Q_from_fit <- function(Znew, Win_new) {
  Xnew <- cbind(Znew, Win_new)
  Hnew <- make_newx_from_basis(hal_zw, Xnew)
  as.numeric(predict(fity, newx = Hnew, type = "response"))
}

## --- EIF components ---------------------------------------------------------
calc_plugin_eif_terms <- function(Q_model, Y, Z, Win, W, h_star_fun, hw_vec) {
  Q_hat   <- Q_model(Z, Win)
  clever  <- (h_star_fun(Z, W) / hfun(Z, hw_vec)) * (Y - Q_hat)
  
  Q0_hat  <- Q_model(rep(0, length(Z)), Win)
  Q1_hat  <- Q_model(rep(1, length(Z)), Win)
  plugin  <- Q0_hat * (1 - h_star_fun(W = W, z = 1)) + Q1_hat * h_star_fun(W = W, z = 1)
  
  list(clever_term = clever, plugin_term = plugin)
}

## --- Grid of target h* ------------------------------------------------------
lgrid <- 10
grid_matrix <- as.matrix(expand.grid(seq(0, 1, length.out = lgrid),
                                     seq(0, 1, length.out = lgrid)))
colnames(grid_matrix) <- c("h_star_0", "h_star_1")

psi     <- numeric(lgrid^2)
sdbound <- numeric(lgrid^2)

## --- Loop over target h* ----------------------------------------------------
for (i in seq_len(lgrid^2)) {
  print(i)
  h_star_0 <- grid_matrix[i, 1]
  h_star_1 <- grid_matrix[i, 2]
  
  # h*(z|w) as a function in (z,w)
  h_star_fun <- function(z, W) z * make_h(W, h_star_0, h_star_1) +
    (1 - z) * (1 - make_h(W, h_star_0, h_star_1))
  
  # Initial Q
  Qn0 <- function(Znew, Win_new) Q_from_fit(Znew, Win_new)
  
  # Clever covariate weights
  weights <- h_star_fun(Z, W) / hfun(Z, hw)
  
  # One-dimensional fluctuation (TMLE)
  offset_eta <- qlogis(Qn0(Z, Win))
  fluctuation <- glm(Y ~ offset(offset_eta) + 1, weights = weights,
                     family = binomial)
  eps <- as.numeric(coef(fluctuation))
  
  # Updated Q
  Qn <- (function(Qold, eps_) {
    function(Znew, Win_new) plogis(qlogis(Qold(Znew, Win_new)) + eps_)
  })(Qn0, eps)
  
  # Parameter map G(Q; h*)
  G_comp_z <- function(Q_model, h0, h1) {
    hW <- make_h(W, h0, h1)
    Q0 <- Q_model(rep(0, n), Win)
    Q1 <- Q_model(rep(1, n), Win)
    mean(Q0 * (1 - hW) + Q1 * hW)
  }
  
  psi[i] <- G_comp_z(Qn, h_star_0, h_star_1)
  
  # EIF & std bound
  terms <- calc_plugin_eif_terms(Qn, Y, Z, Win, W, h_star_fun, hw)
  eif   <- terms$clever_term + terms$plugin_term - psi[i]
  sdbound[i] <- sqrt(var(eif) / n)
}

grid_matrix[13,] # Policy 1
grid_matrix[88,] # Policy 2
grid_matrix[1,] # Intervention h* such that g(h*)=0

psi[13]*100
sdbound[13]*100

psi[88]*100
sdbound[88]*100

psi[1]*100
sdbound[1]*100

# Implied intervention 1
gmat[2,1]*grid_matrix[13,1]*100
gmat[2,2]*grid_matrix[13,2]*100

# Implied intervention 2
gmat[2,1]*grid_matrix[88,1]*100
gmat[2,2]*grid_matrix[88,2]*100





