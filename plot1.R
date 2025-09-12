# install.packages(c("ggplot2", "patchwork", "colorspace"))  # if needed
library(ggplot2)
library(patchwork)
library(colorspace)



# --- Map from h* = (h0, h1) to implied = (im0, im1) ---
map_to_implied <- function(h0, h1) {
  im0 <- g_fun(0,0) * (1 - h0) + g_fun(1,0) * h0  # depends on h0 only
  im1 <- g_fun(0,1) * (1 - h1) + g_fun(1,1) * h1  # depends on h1 only
  c(im0, im1)
}

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
  theme(panel.grid = element_blank())

# --- Panel B: Filled "surface" of the image square ---
p_img <- ggplot(grid, aes(im0, im1)) +
  geom_raster(aes(fill = col)) +
  scale_fill_identity() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(title = "Implied interventions", x = "g(h*)(W=0)", y = "g(h*)(W=0)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

final_plot=(p_src | p_img)

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
