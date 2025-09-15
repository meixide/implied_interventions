set.seed(1)
setwd("~/Projects/implied_interventions")
source('load_example.R')

# install.packages(c("ggplot2", "patchwork", "colorspace"))  # if needed
library(ggplot2)
library(patchwork)
library(colorspace)

gstar0=0.6
gstar1=0.7


map_inverse <- function(im0, im1) {
  # corner constants
  g00 <- g_fun(0,0)
  g10 <- g_fun(1,0)
  g01 <- g_fun(0,1)
  g11 <- g_fun(1,1)
  
  # invert the linear relations
  h0 <- (im0 - g00) / (g10 - g00)
  h1 <- (im1 - g01) / (g11 - g01)
  
  c(h0, h1)
}
hstar0=map_inverse(gstar0,gstar1)[1]
hstar1=map_inverse(gstar0,gstar1)[2]




data=simulate_data(1e6, h0_fun, g_fun) 

Y=data$Y
A=data$A
W=data$W
Z=data$Z

mean(Y)


econd= function(a,w) {
  print(sd(data$Y[data$A==a & data$W==w]))
  return(  mean(data$Y[data$A==a & data$W==w])
  )
}

econdz= function(z,w)
  mean(Y[Z==z & W==w])


econd(1,0)




# g_fun: function(z, w) returning P(A=1 | Z=z, W=w)
econdz <- function(z, w) {
  2 * g_fun(z, w) + w - 0.5
}
gcomp_z=(econdz(1,0)*hstar0 + econdz(0,0)*(1-hstar0))*(1-pw) + (econdz(1,1)*hstar1 + econdz(0,1)*(1-hstar1))*pw


gcomp_z
g_star_fun <- function(w) ifelse(w == 0, gstar0, gstar1)
h_star_fun=function(w) ifelse(w == 0, hstar0, hstar1)

g_fun(1,0)*hstar0 + 
  g_fun(0,0)*(1-hstar0)  

g_fun(1,1)*hstar1 + 
  g_fun(0,1)*(1-hstar1)  

# Simulate intervention (counterfactual) data
simulate_treatment_intervention <- function(n, h_star_fun,g_star_fun) {
  W <- rbinom(n, 1, pw)
  U = runif(n,0,1)
  Z_int <- rbinom(n, 1, h_star_fun(W))
  U = runif(n,0,1)
  A_hstar <- as.numeric(U < g_fun(Z_int, W))
  A_gstar <- as.numeric(U < g_star_fun(W))
  Y_hstar <- rnorm(n, mean = out(A_hstar,W,hid_conf*U), sd = noise_sd)
  Y_gstar <- rnorm(n, mean = out(A_gstar,W,hid_conf*U), sd = noise_sd)
  
  list( W=W,A_hstar = A_hstar, Y_hstar = Y_hstar,A_gstar = A_gstar, Y_gstar = Y_gstar)
}


# --- Map from h* = (h0, h1) to implied = (im0, im1) ---
map_to_implied <- function(h0, h1) {
  im0 <- g_fun(0,0) * (1 - h0) + g_fun(1,0) * h0  # depends on h0 only
  im1 <- g_fun(0,1) * (1 - h1) + g_fun(1,1) * h1  # depends on h1 only
  c(im0, im1)
}

map_to_implied(h0_0,h0_1)



inter=simulate_treatment_intervention(1e6, h_star_fun,g_star_fun)
mean(inter$Y_hstar)
mean(inter$Y_gstar)

mean(Y)

mean(inter$A_gstar[inter$W==0])
mean(inter$A_gstar[inter$W==1])

mean(inter$A_hstar[inter$W==0])
mean(inter$A_hstar[inter$W==1])



# --- High-resolution grid (SURFACES) ---
n <- 401  # dense: true “surface” look
grid <- expand.grid(h0 = seq(0, 1, length.out = n),
                    h1 = seq(0, 1, length.out = n))

mapped <- t(apply(grid, 1, \(r) map_to_implied(r["h0"], r["h1"])))
grid$im0 <- mapped[, 1]
grid$im1 <- mapped[, 2]

# --- Bivariate, non-periodic colors that pop (HCL sweep) ---
grid$col <- hcl(
  h = 200 + 120*grid$h1,   # hue varies with h0 (no wraparound)
  c = 30  + 40*grid$h0,    # chroma increases with h1
  l = 35  + 55*grid$h0     # luminance increases with h1
)

# --- Panel A: Filled "surface" of the source square ---
p_src <- ggplot(grid, aes(h0, h1)) +
  geom_raster(aes(fill = col)) +
  scale_fill_identity() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(title = "Instrumental interventions", x = "h*(W=0)", y = "h*(W=1)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank()) + 
geom_point(aes(x = h0_0, y = h0_1), 
           color = "green", 
           size = 3) + 
  geom_point(aes(x = hstar0, y = hstar1), 
             color = "green", shape=17,
             size = 3)



# --- Panel B: Filled "surface" of the image square ---
p_img <- ggplot(grid, aes(im0, im1)) +
  geom_raster(aes(fill = col)) +
  scale_fill_identity() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(title = "Implied interventions", x = "g(h*)(W=0)", y = "g(h*)(W=1)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank()) + 
  geom_point(aes(x = map_to_implied(h0_0,h0_1)[1], y = map_to_implied(h0_0,h0_1)[2]), 
             color = "green", 
             size = 3) +
geom_point(aes(x = gstar0, y = gstar1), shape=17,
           color = "green", 
           size = 3)

final_plot=(p_src | p_img)

final_plot

# Save as PNG for LaTeX
ggsave(
  filename = "/Users/cgmeixide/Dropbox/Aplicaciones/Overleaf/causal_inference_implied_interventions/gmap.png",
  plot     = final_plot,
  width    = 6.5,     # inches (suitable for a 2-column article)
  height   = 5.0,     # adjust proportionally
  units    = "in",
  dpi      = 600,     # high resolution
  device   = "png"
)


mout=lm(Y ~ A + Z + W)


estar=function(a,w)
  predict(mout,newdata=data.frame(A=a,W=w,Z=1))*g_fun(1,w)*h_star_fun(w) + predict(mout,newdata=data.frame(A=a,W=w,Z=0))*g_fun(0,w)*(1-h_star_fun(w))


outerg=function(w) estar(1,w) + estar(0,w)

outerg(1)*pw + outerg(0)*(1-pw)

gcomp
EYstarSCM


## h_fun:   w -> P(Z=1 | W=w)
## g_fun: (z,w) -> P(A=1 | Z=z, W=w)
## Returns E[Y | A=a, W=w] (vectorized in a, w)
EY_given_AW <- function(a, w, h_fun, g_fun) {
  h  <- h_fun(w)
  g1 <- g_fun(1, w)
  g0 <- g_fun(0, w)
  
  # E[U | A=1, W=w]
  EU_A1 <- (h * g1^2 + (1 - h) * g0^2) /
    (2 * (h * g1 + (1 - h) * g0))
  
  # E[U | A=0, W=w]
  EU_A0 <- (h * (1 - g1^2) + (1 - h) * (1 - g0^2)) /
    (2 * (h * (1 - g1) + (1 - h) * (1 - g0)))
  
  # E[Y | A=a, W=w] = 2a + w - E[U | A=a, W=w]
  ifelse(a == 1, 2 + w - EU_A1, w - EU_A0)
}

## Convenience wrappers for each arm
EY_A1_given_W <- function(w, h_fun, g_fun) {
  h  <- h_fun(w); g1 <- g_fun(1, w); g0 <- g_fun(0, w)
  EU <- (h * g1^2 + (1 - h) * g0^2) / (2 * (h * g1 + (1 - h) * g0))
  2 + w - EU
}

EY_A0_given_W <- function(w, h_fun, g_fun) {
  h  <- h_fun(w); g1 <- g_fun(1, w); g0 <- g_fun(0, w)
  EU <- (h * (1 - g1^2) + (1 - h) * (1 - g0^2)) /
    (2 * (h * (1 - g1) + (1 - h) * (1 - g0)))
  w - EU
}



EY_A1_given_W(0, h0_fun, g_fun)  # E[Y | A=1, W=0]
EY_A0_given_W(1, h0_fun, g_fun)  # E[Y | A=0, W=1]

# Vectorized over (A,W):
A <- c(0,1,1,0); W <- c(0,0,1,1)
econd=function(A,W) 
  EY_given_AW(A, W, h0_fun, g_fun)

econd(0,1)
econd(1,1)
econd(1,0)
econd(0,0)



