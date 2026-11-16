# General PLOS formatting and export for all figure types.
# Source from the project root: source("function/format_PLOS.R")
# Keep plot-specific fonts, guides, aspect ratios and panel assembly in a
# separate layout file or in each figure script.

required <- c("ggplot2", "grid", "sysfonts", "showtext", "gtable")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing packages first: ", paste(missing, collapse = ", "))
}
invisible(lapply(required, library, character.only = TRUE))

FONT_PT <- 8:12
BASE_FAM <- "Arial"
FIG_W_CM <- 19.05
# Default export shape only; override aspect in save_plos() for each figure.
FIG_ASPECT <- 1.0
EXPORT_DPI <- 300

showtext::showtext_opts(dpi = EXPORT_DPI)
showtext::showtext_auto(TRUE)
cm_to_in <- function(x) x / 2.54

# Shared panel theme: all text sizes are set here, including legend text.
# Use theme_pub in each panel instead of overriding it with theme_bw().
theme_pub <- theme_bw(base_size = FONT_PT[4], base_family = BASE_FAM) +
  theme(
    text             = element_text(size = FONT_PT[4]),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    panel.border     = element_rect(fill = NA, colour = "black", linewidth = 2*0.25),
    axis.ticks       = element_line(linewidth = 0.25),
    axis.text        = element_text(size = FONT_PT[3]),
    axis.title       = element_text(size = FONT_PT[4]),
    axis.title.x     = element_text(margin = margin(t = 4)),
    axis.title.y     = element_text(margin = margin(t = 2, r = 4, b = 2, l = 2)),
    plot.title       = element_text(
      size = FONT_PT[5], face = "bold", hjust = 0,
      margin = margin(t = 3, r = 4, b = 2)
    ),
    plot.subtitle    = element_text(size = FONT_PT[4]),
    plot.caption     = element_text(size = FONT_PT[2]),
    plot.tag         = element_text(size = FONT_PT[5]),
    strip.text       = element_text(size = FONT_PT[3]),
    legend.background = element_blank(),
    legend.title     = element_text(size = FONT_PT[3]),  # 10 pt
    legend.text      = element_text(size = FONT_PT[2]),  # 9 pt
    plot.margin      = margin(4, 4, 4, 4)
  )
theme_set(theme_pub)

# geom defaults
update_geom_defaults("line",    list(linewidth = 0.4))
update_geom_defaults("segment", list(linewidth = 0.4))
update_geom_defaults("point",   list(size = 1.4))


# Apply once, to the final assembled figure; use the original figure for each save.
plos_outer_border <- function(plot) {
  g <- if (inherits(plot, "gtable")) plot else ggplot2::ggplotGrob(plot)
  g <- gtable::gtable_add_rows(g, unit(2, "pt"), pos = 0)
  g <- gtable::gtable_add_rows(g, unit(2, "pt"), pos = -1)
  g <- gtable::gtable_add_cols(g, unit(2, "pt"), pos = 0)
  gtable::gtable_add_cols(g, unit(2, "pt"), pos = -1)
}

# Export both formats with identical dimensions and 2-point outer whitespace.
# filename is a stem WITHOUT extension, e.g. "figures/figure6".
# PDF is the companion preview; TIFF is the submission format.
# Choose the final submission filename to match the journal/caption numbering.
save_plos <- function(plot, filename,
                      width_cm = FIG_W_CM, aspect = FIG_ASPECT,
                      dpi = EXPORT_DPI) {
  height_cm <- width_cm * aspect
  stopifnot(length(filename) == 1L, nzchar(filename),
            is.finite(width_cm), is.finite(height_cm),
            width_cm >= 6.68, width_cm <= 19.05,
            height_cm > 0, height_cm <= 22.23,
            dpi >= 300, dpi <= 600)
  if (!capabilities("cairo")) stop("Cairo support is required for these exports.")
  if (grepl("\\.(pdf|tif|tiff)$", filename, ignore.case = TRUE)) {
    stop("Supply a filename stem without .pdf or .tif.")
  }
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  plot_export <- plos_outer_border(plot)
  pdf_file <- paste0(filename, ".pdf")
  tif_file <- paste0(filename, ".tif")
  showtext::showtext_opts(dpi = dpi)

  ggplot2::ggsave(pdf_file, plot = plot_export, device = grDevices::cairo_pdf,
                 width = width_cm, height = height_cm, units = "cm",
                 bg = "white", limitsize = FALSE)
  ggplot2::ggsave(tif_file, plot = plot_export, device = "tiff", type = "cairo",
                 width = width_cm, height = height_cm, units = "cm", dpi = dpi,
                 compression = "lzw", bg = "white", limitsize = FALSE)
  size_mb <- file.info(tif_file)$size / 1e6
  message("Wrote: ", pdf_file)
  message(sprintf("Wrote: %s (%.2f MB)", tif_file, size_mb))
  if (size_mb >= 10) stop("TIFF is at or above 10 MB; reduce it before submission.")
  invisible(c(pdf = pdf_file, tiff = tif_file))
}
