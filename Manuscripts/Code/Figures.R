# ------------------------------------------------------------
# L-shaped function f(x) and its curvature kappa(x)
#
# f(x) is a smooth "softplus" curve that behaves like -x for
# x << 0 (steep/vertical branch) and flattens to 0 for x >> 0
# (flat/horizontal branch) — i.e. an L-shape.
#
# Curvature of a curve y = f(x) is:
#     kappa(x) = |f''(x)| / (1 + f'(x)^2)^(3/2)
# ------------------------------------------------------------


rm(list=ls())

#### Set up  ####
source("Simulation/Code/SimHelpers.r")
source("Simulation/Code/SimData.R")
source("Simulation/Code/SimFit.R")

library(ggplot2)
theme_set(theme_minimal())

#### Example L curve ####

k <- 1.5   # sharpness of the corner (larger k = sharper bend)

f  <- function(x) -(1/k) * log(1 + exp(k * x))       # the L-shaped function
f1 <- function(x) -1 / (1 + exp(k * x))               # first derivative f'(x)
f2 <- function(x) k * exp(k * x) / (1 + exp(k * x))^2 # second derivative f''(x)

kappa <- function(x) abs(f2(x)) / (1 + f1(x)^2)^1.5   # curvature formula

#  Find x that maximizes curvature 
opt <- optimize(kappa, interval = c(-10, 10), maximum = TRUE)
x_star <- opt$maximum
k_star <- opt$objective

cat("Maximum curvature at x =", round(x_star, 4),
    " with kappa =", round(k_star, 4), "\n")
cat("f(x*) =", round(f(x_star), 4), "\n")

# Build data frame for plotting 
x <- seq(-5, 5, length.out = 1000)
df <- data.frame(x = x, y = f(x), kap = kappa(x))

# Top panel: the L-shaped function 
pdf("Manuscripts/Figures/example_lcurve.pdf", height = 4, width = 4)
ggplot(df, aes(x, y)) +
  geom_line(linewidth = 1.2) +
  geom_point(data = data.frame(x = x_star, y = f(x_star)),
             aes(x, y), color = "red", size = 3) +
  annotate("text", x = x_star, y = f(x_star),
           label = "Corner",
           color = "red", fontface = "bold", hjust = -0.1, vjust = 1.1,
           size = 3.4) +
  annotate("text", x = -4, y = f(-4)-0.5,
           label = "Unfair trade-off",
           color = "steelblue", hjust = 0, size = 3.2) +
  annotate("text", x = 3.5, y = f(3.5) + 0.15,
           label = "Unfair trade-off",
           color = "Steelblue", hjust = 0.5, size = 3.2) +
  labs(title = " ", x = "Regularization", y = "Stress") +
  theme(axis.text = element_blank())
dev.off()

#### Example from simulation ####
library(tidyverse)

output_i <- read_rds("Simulation/Output/out_smooth_p10_T10.rds")

best_lam <- c(
  menger = menger_corner(output_i, plot = F)$lambda_star,
  perp_dist = lcurve_corner_dist(output_i)$lambda_star)


output_i <- output_i %>%
  mutate_at(2:4, log)

# KEY CHANGE: rescale movement and stress to [0, 1]
# so that both axes are on the same numerical scale, and a
# perpendicular line computed via dot-product projection will
# actually render as perpendicular under default (non-fixed) coords
range01 <- function(x) (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))

output_i <- output_i %>%
  mutate(movement = range01(movement),
         stress   = range01(stress))


# point coordinates
max_dist_point <- output_i[output_i$lambda == best_lam["perp_dist"],
                           c("movement", "stress")]
menger_point <- output_i[output_i$lambda == best_lam["menger"],
                         c("movement", "stress")]
points <- rbind(max_dist_point, menger_point)
points$method <- c("Max dist", "Menger")

# get first and last points (defines the blue reference line)
first_pt <- output_i %>% slice(1)
last_pt  <- output_i %>% slice(n())

x1 <- first_pt$movement; y1 <- first_pt$stress
x2 <- last_pt$movement;  y2 <- last_pt$stress

# max dist point
xm <- points$movement[1]
ym <- points$stress[1]

# compute the perpendicular foot (projection of max dist point onto the line)
dx <- x2 - x1;  dy <- y2 - y1
t  <- ((xm - x1) * dx + (ym - y1) * dy) / (dx^2 + dy^2)
xf <- x1 + t * dx
yf <- y1 + t * dy


# plot
pdf("Manuscripts/Figures/simulation_lcurve.pdf", height = 4, width = 4)
output_i %>%
  ggplot(aes(x = movement, y = stress)) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_point(data = points, aes(x = movement, y = stress),
             size = 2, col = "red") +
  geom_label(data = points, aes(x = movement, y = stress, label = method),
             nudge_y = 0.05,
             nudge_x = 0.1,
             color   = "red",
             size    = 3) +
  # 1. straight line from first to last point
  annotate("segment",
           x = x1, y = y1, xend = x2, yend = y2,
           color = "blue", linetype = "dashed") +
  # 2. perpendicular line from max dist point to the blue line
  annotate("segment",
           x = xm, y = ym, xend = xf, yend = yf,
           color = "darkgreen", linetype = "dotted") +
  # force equal visual scaling on both axes so the right angle
  # is actually rendered as a right angle
  coord_fixed(ratio = 1)+
  labs(x = "Regularization (normalized)", y = "Stress (normalized)")
dev.off()

#### IFED network plot ####


