################################################################################
# Copyright (C) 2025  Michiel van Boven
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or any later version.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
################################################################################

# packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)
library(grid)
library(ggrastr)

# sizes
FIG_W_CM   <- 8.7        # PNAS single-column width
FIG_ASPECT <- 1.0        # default aspect
BASE_PT    <- 9          # base font size (pt)
BASE_FAM   <- "CMUSerif" # or other serif
cm_to_in <- function(x) x / 2.54

# simple scientific format
sci_fmt <- function(x) {
  fmt <- format(x, scientific = TRUE)
  gsub("e\\+0?", "e", fmt)
}

# global theme
theme_pub <- theme_minimal(base_size = BASE_PT, base_family = BASE_FAM) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    axis.ticks       = element_line(linewidth = 0.25),
    axis.title       = element_text(margin = margin(t = 2, r = 2, b = 2, l = 2)),
    plot.title       = element_text(face = "bold"),
    plot.margin      = margin(4, 4, 4, 4),
    aspect.ratio     = FIG_ASPECT
  )
theme_set(theme_pub)

# geom defaults
update_geom_defaults("line",    list(linewidth = 0.4))
update_geom_defaults("segment", list(linewidth = 0.4))
update_geom_defaults("point",   list(size = 1.4))

# pdf saver (PNAS style)
save_pdf <- function(plot, file_pdf,
                     width_cm = FIG_W_CM,
                     aspect   = FIG_ASPECT,
                     dpi      = 300) {
  w_in <- cm_to_in(width_cm)
  h_in <- w_in * aspect
  ggplot2::ggsave(
    filename  = file_pdf,
    plot      = plot,
    width     = w_in,
    height    = h_in,
    units     = "in",
    dpi       = dpi,
    device    = cairo_pdf,
    limitsize = FALSE
  )
  message("Wrote: ", file_pdf)
}

#################################################
# Figure 4. Time plots for influenza A scenario #
#################################################

#df_ts <- read_csv("data/time_series_influenza_2025_20251117.csv")
#df_ta <- read_csv("data/time_age_influenza_2025_20251117.csv")
df_ts <- read_csv("data/time_series_influenza_2025_20251222.csv")
df_ta <- read_csv("data/time_age_influenza_2025_20251222.csv")

t_ref <- 150

df_ts <- df_ts %>%
  mutate(t_rel = t - t_ref)

# colours
col_low    <- "#56B4E9"
col_mid    <- "#E69F00"
col_high   <- "red"
col_low_w  <- "white" # "#c7e4fb"

# time-series theme
theme_time <- theme_bw(base_size = BASE_PT, base_family = BASE_FAM) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.25),
    axis.title       = element_text(margin = margin(t = 2, r = 2, b = 2, l = 2)),
    plot.margin      = margin(2, 10, 4, 6)
  )

remove_x_all <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

# top panel: compartments over time
df_top <- df_ts %>%
  select(t_rel, S, I, R, V) %>%
  pivot_longer(-t_rel, names_to = "compartment", values_to = "value") %>%
  mutate(
    compartment = factor(compartment, levels = c("S", "I", "R", "V"))
  )

p_top <- ggplot() +
  geom_line(
    data = dplyr::filter(df_top, compartment != "I"),
    aes(x = t_rel, y = value, colour = compartment),
    linewidth = 0.4
  ) +
  geom_line(
    data = dplyr::filter(df_top, compartment == "I"),
    aes(x = t_rel, y = value, colour = compartment),
    linewidth = 0.4
  ) +
  scale_colour_manual(
    values = c(
      S = col_low,
      I = col_mid,
      R = "grey20",
      V = "grey70"
    ),
    breaks = c("S", "I", "R", "V"),  # <- force legend order
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    breaks = seq(0, 100, 20)
  ) +
  scale_y_continuous(labels = sci_fmt) +
  labs(
    x = NULL,
    y = "Number"
  ) +
  theme_time +
  remove_x_all +
  theme(
    legend.position      = c(0.02, 0.97),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = "grey80"),
    legend.text          = element_text(size = BASE_PT - 1),
    legend.margin        = margin(t = 1, r = 5, b = 2, l = 2),
    legend.key.width     = unit(0.6, "lines"),
    legend.key.height    = unit(0.6, "lines")
  )

# middle panel: age-time incidence heatmap
df_mid <- df_ta %>%
  mutate(t_rel = t - t_ref) %>%
  pivot_longer(
    cols = starts_with("age_"),
    names_to = "age_label",
    values_to = "I"
  ) %>%
  mutate(
    age = as.numeric(sub("age_", "", age_label))
  ) %>%
  arrange(t_rel, age)

t_vals   <- sort(unique(df_mid$t_rel))
dt_vals  <- diff(t_vals)
dt       <- max(dt_vals)

age_vals <- sort(unique(df_mid$age))
da_vals  <- diff(age_vals)
da       <- max(da_vals)

I_min <- min(df_mid$I, na.rm = TRUE)
I_max <- max(df_mid$I, na.rm = TRUE)

# warning: large file. try rasterisation
p_mid <- ggplot(df_mid, aes(x = t_rel, y = age, fill = I)) +
#  ggrastr::rasterise(
#    geom_tile(width = dt * 1.02, height = da * 1.02),
#    dpi = 150   
#  ) +
  geom_tile(width = dt * 1.02, height = da * 1.02) +
  geom_text(
    inherit.aes = FALSE,
    aes(x = 98, y = 115, label = "italic(I)(t,a)"),
    parse = TRUE,
    family = BASE_FAM,
    size   = 3.2,
    hjust  = 1,
    vjust  = 1
  ) +
  scale_y_continuous(
    name   = "Age (yr)",
    expand = c(0, 0),
    limits = range(df_mid$age),
    breaks = seq(0, 120, 20)
  ) +
  scale_x_continuous(
    name   = NULL,
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_fill_gradient2(
    low      = col_low_w,
    mid      = col_mid,
    high     = col_high,
    midpoint = (I_min + I_max) / 2,
    limits   = c(I_min, I_max),
    name     = NULL,
    labels   = sci_fmt,
    guide    = guide_colorbar(
      title          = NULL,
      barheight      = unit(2.0, "cm"),
      barwidth       = unit(0.25, "cm"),
      ticks          = TRUE,
      label.position = "right"
    )
  ) +
  theme_time +
  remove_x_all +
  theme(
    legend.position      = c(0.02, 0.97),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = "grey80"),
    legend.text          = element_text(size = BASE_PT - 3)
  )

# bottom panel: mortality over time
p_bottom <- ggplot(df_ts, aes(x = t_rel, y = M_inc)) +
  geom_line(linewidth = 0.4, colour = col_high) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(0, 1e4, 2e4),
    labels = sci_fmt
  ) +
  labs(
    x = "Time (yr)",
    y = expression(plain(Mort) ~ plain("(") * yr^{-1} * plain(")"))
  ) +
  theme_time +
  theme(
    legend.position = "none",
    plot.margin     = margin(t = 6, r = 10, b = 2, l = 6)
  )

# combine panels
width_cm  <- FIG_W_CM
height_cm <- FIG_W_CM * 1.05
aspect_three <- height_cm / width_cm

figure_4 <- cowplot::plot_grid(
  p_top,
  p_mid,
  p_bottom,
  ncol           = 1,
  labels         = c("A", "B", "C"),
  label_fontface = "bold",
  label_size     = BASE_PT + 1,
  label_x        = 0.015,
  label_y        = 0.985,
  align          = "v",
  axis           = "l",
  rel_heights    = c(1, 1, 1)
)

figure_4

#library(ragg)  # optional since you used ragg::agg_png, but fine
save_png <- function(plot, file_png,
                     width_cm = FIG_W_CM,
                     aspect   = 1.05,
                     dpi      = 160) {
  w_in <- width_cm / 2.54
  h_in <- w_in * aspect
  ragg::agg_png(file_png, width = w_in, height = h_in, units = "in", res = dpi)
  print(plot)
  dev.off()
}
save_png(figure_4, "figures/figure4_overleaf.png", width_cm = FIG_W_CM, aspect = aspect_three, dpi = 160)

save_pdf(
  plot     = figure_4,
  file_pdf = "figures/figure4_150dpi.pdf",
  width_cm = width_cm,
  aspect   = aspect_three
)
