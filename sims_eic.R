source('load_example.R')

EYstar=mean(simulate_intervention(1e6, function(w) make_h(w, h_star_0 , h_star_1 ), g_fun)$Y_int)
EYstar

B=1000

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


