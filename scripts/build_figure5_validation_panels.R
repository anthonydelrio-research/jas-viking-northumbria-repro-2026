#!/usr/bin/env Rscript

# Build Figure 5 (validation diagnostics panels) from JAS outputs.
# No external R packages are required (base R only).
# The script works with either:
#  - public package layout: data/derived/*
#  - author full-run layout: outputs_submission_* / outputs_750_*

args <- commandArgs(trailingOnly = TRUE)
jas_root <- if (length(args) >= 1) args[[1]] else normalizePath(".", mustWork = TRUE)

first_existing <- function(paths, required = TRUE) {
  for (p in paths) {
    if (file.exists(p) || dir.exists(p)) return(p)
  }
  if (required) {
    stop(
      paste("Could not resolve required path from candidates:", paste(paths, collapse = "; "))
    )
  }
  return(NA_character_)
}

path_primary_tables <- first_existing(c(
  file.path(jas_root, "data", "derived", "primary", "tables"),
  file.path(jas_root, "outputs_750_primary_pasall_v1", "tables"),
  file.path(jas_root, "outputs_submission_2026-03-05_pasall_context50", "tables")
))
path_sensitivity_tables <- first_existing(c(
  file.path(jas_root, "data", "derived", "sensitivity", "tables"),
  file.path(jas_root, "outputs_750_sensitivity_pasall_v2", "tables"),
  file.path(jas_root, "outputs_submission_2026-03-05_temporal_weighted", "tables")
))
path_holdout_stress <- first_existing(c(
  file.path(jas_root, "data", "derived", "holdout", "HOLDOUT_BREADTH_STRESS_TEST.tsv"),
  file.path(jas_root, "outputs_submission_2026-03-05_blocker_resolution_final", "HOLDOUT_BREADTH_STRESS_TEST.tsv")
))
path_holdout_points <- first_existing(c(
  file.path(jas_root, "data", "derived", "holdout", "validation_holdout_sites_used_with_top15.csv"),
  file.path(jas_root, "figures", "generated", "validation_holdout_sites_used_with_top15.csv")
), required = FALSE)
path_output_dir <- file.path(jas_root, "outputs", "figures")

dir.create(path_output_dir, showWarnings = FALSE, recursive = TRUE)

req_files <- c(
  file.path(path_primary_tables, "decile_lift_table.tsv"),
  file.path(path_primary_tables, "auc_metrics.tsv"),
  file.path(path_sensitivity_tables, "auc_metrics.tsv"),
  path_holdout_stress
)

missing_files <- req_files[!file.exists(req_files)]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required input file(s):",
      paste(missing_files, collapse = "; ")
    )
  )
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

decile <- read_tsv(file.path(path_primary_tables, "decile_lift_table.tsv"))
auc_primary <- read_tsv(file.path(path_primary_tables, "auc_metrics.tsv"))
auc_sensitivity_tbl <- read_tsv(file.path(path_sensitivity_tables, "auc_metrics.tsv"))
holdout_stress <- read_tsv(path_holdout_stress)
holdout_pts <- NULL
if (!is.na(path_holdout_points) && file.exists(path_holdout_points)) {
  holdout_pts <- read.csv(path_holdout_points, stringsAsFactors = FALSE)
}

auc_sensitivity_stress <- holdout_stress$auc_uniform[
  holdout_stress$holdout_set == "primary_plus_sensitivity"
]
if (length(auc_sensitivity_stress) != 1) {
  stop("Could not extract unique primary_plus_sensitivity AUC from holdout stress table.")
}
auc_sensitivity <- as.numeric(auc_sensitivity_tbl$auc_mean[1])
auc_sensitivity_w <- as.numeric(auc_sensitivity_tbl$auc_temporal_weighted[1])

if (!is.null(holdout_pts)) {
  holdout_pts$site_type <- tolower(holdout_pts$site_type)
  holdout_pts$capture_status <- ifelse(holdout_pts$in_top15 == 1, "Captured (Top 15%)", "Not captured")
}

plot_panel_a <- function() {
  x <- decile$top_area_pct
  y_model <- decile$capture_rate * 100
  y_random <- decile$expected_random_rate * 100
  y_opp <- decile$opportunity_capture_rate * 100

  plot(
    x, y_model,
    type = "b",
    lwd = 2.2,
    pch = 19,
    col = "#1B9E77",
    xlim = c(10, 100),
    ylim = c(0, 100),
    xlab = "Top ranked area (%)",
    ylab = "Captured holdouts (%)",
    main = "Hit@Area curves"
  )
  lines(x, y_random, lwd = 2, lty = 2, col = "#4D4D4D")
  lines(x, y_opp, lwd = 2, lty = 3, col = "#D95F02")
  points(x, y_opp, pch = 17, cex = 0.8, col = "#D95F02")
  legend(
    "topleft",
    bty = "n",
    cex = 0.82,
    legend = c("SICI model", "CSR baseline", "Opportunity baseline"),
    col = c("#1B9E77", "#4D4D4D", "#D95F02"),
    lty = c(1, 2, 3),
    pch = c(19, NA, 17),
    lwd = c(2.2, 2, 2)
  )
}

plot_panel_b <- function() {
  labels <- c(
    "Primary\nAUC",
    "Primary\nAUCw",
    "Sensitivity\nAUC",
    "Sensitivity\nAUCw",
    "Stress test\nAUC"
  )
  vals <- c(
    as.numeric(auc_primary$auc_mean[1]),
    as.numeric(auc_primary$auc_temporal_weighted[1]),
    as.numeric(auc_sensitivity),
    as.numeric(auc_sensitivity_w),
    as.numeric(auc_sensitivity_stress[1])
  )
  cols <- c("#2C7BB6", "#2C7BB6", "#ABD9E9", "#ABD9E9", "#FDAE61")

  bp <- barplot(
    vals,
    names.arg = labels,
    col = cols,
    border = "#333333",
    ylim = c(0.5, 0.75),
    ylab = "AUC",
    main = "Temporal and sensitivity diagnostics",
    las = 1
  )
  abline(h = 0.5, lty = 2, col = "#7F7F7F")
  text(
    x = bp,
    y = vals + 0.006,
    labels = sprintf("%.3f", vals),
    cex = 0.8
  )
  legend(
    "topleft",
    bty = "n",
    cex = 0.8,
    legend = c("Primary run", "Sensitivity run", "Sensitivity stress test"),
    fill = c("#2C7BB6", "#ABD9E9", "#FDAE61")
  )
}

plot_panel_c <- function() {
  if (is.null(holdout_pts)) {
    plot.new()
    title("Primary holdout site distribution")
    text(
      x = 0.5, y = 0.56,
      labels = "Holdout point CSV not included in this package.",
      cex = 0.95
    )
    text(
      x = 0.5, y = 0.46,
      labels = "Panels A-B are fully reproducible from redistributed tables.",
      cex = 0.88
    )
    return(invisible(NULL))
  }

  x_col <- if ("X" %in% names(holdout_pts)) "X" else if ("x" %in% names(holdout_pts)) "x" else NA_character_
  y_col <- if ("Y" %in% names(holdout_pts)) "Y" else if ("y" %in% names(holdout_pts)) "y" else NA_character_
  if (is.na(x_col) || is.na(y_col)) {
    stop("Holdout point file present but missing X/Y (or x/y) columns.")
  }

  x <- holdout_pts[[x_col]]
  y <- holdout_pts[[y_col]]
  cols <- ifelse(holdout_pts$capture_status == "Captured (Top 15%)", "#1B9E77", "#BDBDBD")
  pchs <- ifelse(holdout_pts$site_type == "settlement", 24, 21)

  plot(
    x, y,
    type = "n",
    asp = 1,
    xlab = "Easting (m, EPSG:27700)",
    ylab = "Northing (m, EPSG:27700)",
    main = "Primary holdout site distribution"
  )
  grid(col = "#ECECEC", lty = 1)
  points(
    x, y,
    pch = pchs,
    bg = cols,
    col = "#1F1F1F",
    cex = 1.05
  )

  legend(
    "bottomleft",
    bty = "n",
    cex = 0.8,
    legend = c("Captured (Top 15%)", "Not captured", "Burial", "Settlement"),
    pch = c(21, 21, 21, 24),
    pt.bg = c("#1B9E77", "#BDBDBD", "#FFFFFF", "#FFFFFF"),
    col = c("#1F1F1F", "#1F1F1F", "#1F1F1F", "#1F1F1F")
  )
}

draw_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3.4, 1.3), oma = c(0.2, 0.2, 0.2, 0.2))
  plot_panel_a()
  mtext("A", side = 3, line = 1.0, adj = -0.12, font = 2, cex = 1.1)
  plot_panel_b()
  mtext("B", side = 3, line = 1.0, adj = -0.12, font = 2, cex = 1.1)
  plot_panel_c()
  mtext("C", side = 3, line = 1.0, adj = -0.12, font = 2, cex = 1.1)
}

out_pdf <- file.path(path_output_dir, "FIG05_validation_diagnostics_public.pdf")
out_tif <- file.path(path_output_dir, "FIG05_validation_diagnostics_public.tif")

pdf(out_pdf, width = 13, height = 4.8, useDingbats = FALSE)
draw_figure()
dev.off()

tiff(
  out_tif,
  width = 13,
  height = 4.8,
  units = "in",
  res = 600
)
draw_figure()
dev.off()

cat("Wrote:\n")
cat(out_pdf, "\n")
cat(out_tif, "\n")
