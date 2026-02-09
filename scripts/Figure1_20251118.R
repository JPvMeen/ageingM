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

# sizes
FIG_W_CM   <- 8.7        # PNAS single-column width
FIG_ASPECT <- 1.0        # default aspect
BASE_PT    <- 9          # base font size (pt)
BASE_FAM   <- "CMUSerif" # or other serif

cm_to_in <- function(x) x / 2.54

# simple scientific format (if needed)
sci_fmt <- function(x) {
  fmt <- format(x, scientific = TRUE)
  gsub("e\\+0?", "e", fmt)
}

# global theme
theme_pub <- theme_minimal(base_size = BASE_PT, base_family = BASE_FAM) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
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

# pdf saver (PNAS style, Cairo backend)
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

##########################################
# Data: AG2022 mortality (qx male/female)
##########################################

qx_female_wide <- read_csv("./data/qxfemale2022.csv")
qx_male_wide   <- read_csv("./data/qxmale2022.csv")

qx_female_long <- qx_female_wide %>%
  pivot_longer(!x, names_to = "year", values_to = "rate") %>%
  rename(age = x) %>%
  mutate(across(where(is.character), as.numeric))

qx_male_long <- qx_male_wide %>%
  pivot_longer(!x, names_to = "year", values_to = "rate") %>%
  rename(age = x) %>%
  mutate(across(where(is.character), as.numeric))

# overall mortality (simple mean male/female)
qx_overall <- (qx_male_long + qx_female_long) / 2

# log-mortality and year offset
qx_overall <- qx_overall %>%
  mutate(
    log_rate    = log10(rate),
    year_offset = year - 2020
  )

# fit log10(mu) ~ age + (year - 2020) for ages >= 15
qx_overall_filtered <- qx_overall %>%
  filter(age >= 15)

linear_model <- lm(log_rate ~ age + year_offset, data = qx_overall_filtered)

# coefficients
a0 <- coef(linear_model)["(Intercept)"]
a1 <- coef(linear_model)["year_offset"]
a2 <- coef(linear_model)["age"]

##########################################
# Predicted mortality and survival
##########################################

ages  <- 0:120
years <- 2020:2100

combined_data <- tidyr::expand_grid(age = ages, year = years) %>%
  arrange(year, age) %>%
  group_by(year) %>%
  mutate(
    log_predicted_rate = a0 + a1 * (year - 2020) + a2 * age,
    mortality_rate     = 10^log_predicted_rate,
    dx                 = c(0, diff(age)),
    cumulative_mort    = cumsum(mortality_rate * dx),
    survival_prob      = exp(-cumulative_mort)
  ) %>%
  ungroup()

##########################################
# Figure 1: panels + legend
##########################################

# shared colour scale
color_scale <- scale_colour_gradientn(
  colors = c("#E69F00", "#56B4E9"),
  breaks = c(2020, 2060, 2100),
  labels = c("2020", "2060", "2100"),
  name   = "Year",
  guide  = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(160, "pt"),
    barheight = unit(6,  "pt"),
    ticks = FALSE, label.position = "bottom"
  )
)

ylab_mort_expr <- expression(log[10](mu(t, a)) ~ yr^{-1})
ylab_surv_expr <- expression(Pr(alive))

# mortality panel
mortality_plot <- ggplot(qx_overall %>% filter(year <= 2100)) +
  geom_line(aes(age, log_rate, group = year, colour = year), alpha = 0.18) +
  geom_line(
    data = combined_data,
    aes(age, log_predicted_rate, group = year, colour = year),
    linewidth = 0.3, alpha = 0.9
  ) +
  color_scale +
  scale_x_continuous(
    "Age (years)",
    limits = c(0, 120), expand = c(0, 0),
    breaks = c(0, 25, 50, 75, 100, 120)
  ) +
  scale_y_continuous(
    limits = c(-5, 0), expand = c(0, 0),
    breaks = -5:0
  ) +
  labs(y = ylab_mort_expr, title = "Mortality") +
  theme_bw(base_size = BASE_PT) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    legend.position  = "none",
    plot.margin      = margin(8, 8, 6, 8)
  )

# survival panel
survival_plot <- ggplot(combined_data, aes(age, survival_prob, group = year, colour = year)) +
  geom_line(alpha = 0.18) +
  color_scale +
  scale_x_continuous(
    "Age (years)",
    limits = c(0, 120), expand = c(0, 0),
    breaks = c(0, 25, 50, 75, 100, 120)
  ) +
  scale_y_continuous(
    limits = c(0, 1), expand = c(0, 0),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(y = ylab_surv_expr, title = "Survival") +
  theme_bw(base_size = BASE_PT) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    legend.position  = "none",
    plot.margin      = margin(8, 8, 6, 8)
  )

# legend helper
get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}

# legend grob
legend <- get_legend(
  ggplot(combined_data, aes(age, log_predicted_rate, colour = year)) +
    geom_line(alpha = 0.2) +
    color_scale +
    theme_minimal(base_size = BASE_PT) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_text(margin = margin(b = 3)),
      legend.key.height = unit(6, "pt"),
      legend.key.width  = unit(18, "pt"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

# spacer between panels
spacer <- ggplot() + theme_void()

# horizontal row with A/B labels and small gap
plots_row <- cowplot::plot_grid(
  mortality_plot + theme(plot.margin = margin(8, 2, 6, 2)),
  spacer,
  survival_plot  + theme(plot.margin = margin(8, 2, 6, 2)),
  ncol        = 3,
  rel_widths  = c(1.05, 0.08, 1.00),
  align       = "v",
  axis        = "lr",
  labels      = c("A", "", "B"),
  label_size  = BASE_PT + 1,
  label_fontface = "bold",
  label_x     = c(0.01, NA, 0.01),
  label_y     = c(0.98, NA, 0.98)
)

final_plot <- cowplot::plot_grid(
  plots_row, legend,
  ncol        = 1,
  rel_heights = c(1, 0.22)
)

figure_1 <- cowplot::ggdraw(final_plot) +
  theme(plot.margin = margin(8, 10, 12, 10))

figure_1

# save 
save_pdf(
  plot     = figure_1,
  file_pdf = "figures/figure1.pdf",
  width_cm = FIG_W_CM,
  aspect   = 0.75
)
