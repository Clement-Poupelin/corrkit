# ====================================================================
# data-raw/logo.R -- hex logo of corrkit
# ====================================================================
# Lower-triangular mini heatmap (blue < 0 < red, same palette as the package)
# above the package name, drawn with ggplot2 and assembled with hexSticker.
# Run from the root of the project:   source("data-raw/logo.R")
# Writes man/figures/logo.png (used by the README and the pkgdown site).
# Packages: ggplot2, hexSticker (install.packages("hexSticker")).
# ====================================================================

library(ggplot2)
library(hexSticker)

pal <- c(low = "#2166AC", high = "#B2182B")
navy <- "#1B3A5C"

# 4 x 4 lower triangle of correlations (diagonal = 1)
r <- c(1,
       0.62, 1,
      -0.35, 0.18, 1,
       0.81, -0.58, 0.27, 1)
cells <- data.frame(i = rep(1:4, times = 1:4),
                    j = unlist(lapply(1:4, seq_len)),
                    r = r)

heat <- ggplot(cells, aes(x = j, y = -i, fill = r)) +
  geom_tile(width = 0.88, height = 0.88) +
  scale_fill_gradient2(low = pal[["low"]], mid = "white", high = pal[["high"]],
                       midpoint = 0, limits = c(-1, 1), guide = "none") +
  coord_equal() +
  theme_void() +
  theme_transparent()

sticker(heat,
        package  = "corrkit", p_size = 22, p_y = 0.55, p_color = navy, p_fontface = "bold",
        s_x = 1, s_y = 1.22, s_width = 0.95, s_height = 0.95,
        h_fill = "white", h_color = navy, h_size = 1.6,
        url = "r · rₛ · τ", u_size = 5, u_color = "#5A6B7D", u_x = 1, u_y = 0.3, u_angle = 0,
        filename = "man/figures/logo.png", dpi = 300)

