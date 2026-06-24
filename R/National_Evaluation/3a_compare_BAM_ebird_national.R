# ==============================================================================
# compare_BAM_ebird_national.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1 – Libraries
# ------------------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggrepel)
library(patchwork)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(units)

# ------------------------------------------------------------------------------
# Section 2 – Configuration
# ------------------------------------------------------------------------------
source("R/config.R")
source("R/functions.R")

study_area <- rnaturalearth::ne_states(country = "Canada") %>%
  sf::st_transform(3978) %>%
  sf::st_make_valid()

bcr <- sf::st_read(BCR_GPKG_PATH, quiet = TRUE) %>%
  sf::st_transform(sf::st_crs(study_area)) %>%
  sf::st_make_valid() %>%
  sf::st_intersection(study_area) %>%
  dplyr::filter(!(bcr_label %in% BCR_EXCLUDE))

region_national <- sf::st_union(bcr) %>%
  sf::st_as_sf() %>%
  remove_multipolygon_holes()   # defined in functions_optimized.R



OUT_COMPARISON <- file.path(OUT_ROOT, "National", "comparison")
dir.create(OUT_COMPARISON, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Section 3 – Load results
# ------------------------------------------------------------------------------
bam_stats_path   <- file.path(OUT_NATIONAL_BAM_STATS,
                              "BAM_preds_vs_BAM_data.rds")
ebird_stats_path <- file.path(OUT_NATIONAL_EBIRD_STATS,
                              "eBird_preds_vs_BAM_data.rds")

if (!file.exists(bam_stats_path))
  stop("BAM stats not found at: ", bam_stats_path,
       "\nRun assess_BAM_national.R first.")
if (!file.exists(ebird_stats_path))
  stop("eBird stats not found at: ", ebird_stats_path,
       "\nRun assess_ebird_national.R first.")

bam_stats   <- readRDS(bam_stats_path)
ebird_stats <- readRDS(ebird_stats_path)

message("BAM stats:   ", nrow(bam_stats),   " species")
message("eBird stats: ", nrow(ebird_stats), " species")

# ------------------------------------------------------------------------------
# Section 4 – Join and derive comparison metrics
# ------------------------------------------------------------------------------
combined <- bam_stats %>%
  dplyr::select(
    spcd_bam, sp_english,
    rho_bam          = spearman_obs_pred_above_absence,
    n_surveys_bam    = n_surveys_total,
    n_hex_surv_bam   = n_hex_surveyed,
    n_absent_bam     = n_flag_pred_absent_detected,
    n_highpred_bam   = n_flag_high_pred_no_detection,
    n_flag_any_bam   = n_flag_any,
    prop_flagged_bam = prop_surveyed_flagged
  ) %>%
  dplyr::inner_join(
    ebird_stats %>%
      dplyr::select(
        spcd_bam,
        rho_ebird          = spearman_obs_pred_above_absence,
        n_surveys_ebird    = n_surveys_total,
        n_hex_surv_ebird   = n_hex_surveyed,
        n_absent_ebird     = n_flag_pred_absent_detected,
        n_highpred_ebird   = n_flag_high_pred_no_detection,
        n_flag_any_ebird   = n_flag_any,
        prop_flagged_ebird = prop_surveyed_flagged,
        sp_english_ebird   = sp_english
      ),
    by = "spcd_bam"
  ) %>%
  dplyr::mutate(
    sp_english = dplyr::coalesce(sp_english, sp_english_ebird)
  ) %>%
  dplyr::select(-sp_english_ebird) %>%
  dplyr::mutate(
    # Signed differences (BAM - eBird); positive = BAM higher
    delta_rho       = rho_bam        - rho_ebird,
    delta_absent    = n_absent_bam   - n_absent_ebird,
    delta_highpred  = n_highpred_bam - n_highpred_ebird,
    delta_flag_any  = n_flag_any_bam - n_flag_any_ebird,
    delta_prop_flag = prop_flagged_bam - prop_flagged_ebird,
    
    # z-scores of absolute delta — used only for legacy compatibility
    zdelta_rho      = abs(delta_rho)      / sd(delta_rho,      na.rm = TRUE),
    zdelta_absent   = abs(delta_absent)   / sd(delta_absent,   na.rm = TRUE),
    zdelta_highpred = abs(delta_highpred) / sd(delta_highpred, na.rm = TRUE),
    zdelta_flag_any = abs(delta_flag_any) / sd(delta_flag_any, na.rm = TRUE)
  )

message("Species in both outputs: ", nrow(combined))

# ------------------------------------------------------------------------------
# Section 5 – Win scoring
#
# BAM wins a metric if:
#   rho:       delta_rho       > 0  (higher rho is better)
#   absent:    delta_absent    < 0  (fewer flags is better)
#   highpred:  delta_highpred  < 0
#   prop_flag: delta_prop_flag < 0
# ------------------------------------------------------------------------------
combined <- combined %>%
  dplyr::mutate(
    win_rho = dplyr::case_when(
      delta_rho > 0 ~  1L,
      delta_rho < 0 ~ -1L,
      TRUE          ~  0L
    ),
    win_absent = dplyr::case_when(
      delta_absent < 0 ~  1L,
      delta_absent > 0 ~ -1L,
      TRUE             ~  0L
    ),
    win_highpred = dplyr::case_when(
      delta_highpred < 0 ~  1L,
      delta_highpred > 0 ~ -1L,
      TRUE               ~  0L
    ),
    win_propflag = dplyr::case_when(
      delta_prop_flag < 0 ~  1L,
      delta_prop_flag > 0 ~ -1L,
      TRUE                ~  0L
    ),
    overall_score = win_rho + win_absent + win_highpred
  )

# ------------------------------------------------------------------------------
# Section 6 – Shared theme and helpers
# ------------------------------------------------------------------------------
theme_comparison <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 10),
      plot.caption     = ggplot2::element_text(colour = "grey50", size = 8,
                                               hjust = 0)
    )
}

LABEL_COLOUR   <- "#C0392B"
REF_LINE_COLOUR <- "grey60"

# Build a multi-line summary string shown in the plot corner.
# win_col: name of an integer win column (+1 BAM, -1 eBird, 0 tied).
# n_plotted: number of species actually shown in this panel (after filtering).
make_summary_text <- function(data, win_col, n_plotted = NULL) {
  wins <- data[[win_col]]
  n_bam   <- sum(wins >  0, na.rm = TRUE)
  n_ebird <- sum(wins <  0, na.rm = TRUE)
  n_tied  <- sum(wins == 0, na.rm = TRUE)
  n_shown <- if (!is.null(n_plotted)) n_plotted else sum(!is.na(wins))
  paste0(
    "n shown = ", n_shown, "\n",
    "BAM better:   ", n_bam,   "\n",
    "eBird better: ", n_ebird, "\n",
    "Tied:         ", n_tied
  )
}

# Scatterplot with 1:1 line, discrepant-point labels, and corner summary text.
# summary_text: pre-built character string (use make_summary_text()).
scatter_with_labels <- function(data, x, y, x_label, y_label,
                                zdelta_col, title, subtitle = NULL,
                                summary_text = NULL,
                                log_axes = FALSE,
                                n_label = 25) {
  
  data <- data %>%
    dplyr::mutate(
      .x      = .data[[x]],
      .y      = .data[[y]],
      .zdelta = .data[[zdelta_col]],
      .delta  = .y - .x,
      .rank_above = dplyr::if_else(.delta > 0,
                                   rank(-.delta, ties.method = "first"),
                                   NA_integer_),
      .rank_below = dplyr::if_else(.delta < 0,
                                   rank( .delta, ties.method = "first"),
                                   NA_integer_),
      .label = dplyr::if_else(
        (!is.na(.rank_above) & .rank_above <= n_label) |
          (!is.na(.rank_below) & .rank_below <= n_label),
        sp_english, NA_character_
      )
    )
  
  # Break generator safe across zero: 0, 1, 10, 100, ... up to data max
  pseudo_log_breaks <- function(limits) {
    pos_max <- max(limits[is.finite(limits) & limits > 0], na.rm = TRUE)
    if (length(pos_max) == 0 || !is.finite(pos_max) || pos_max <= 0) {
      return(0)
    }
    pow <- seq(0, floor(log10(pos_max)), by = 1)
    sort(unique(c(0, 10^pow)))
  }
  
  # Compute axis limits in data space; for log panels expand slightly above
  # zero so the 0-vs-0 mass is not clipped by coord_equal
  raw_range <- range(c(data$.x, data$.y), na.rm = TRUE)
  
  if (log_axes) {
    # Ensure lower limit is exactly 0 (never negative)
    axis_lim <- c(0, raw_range[2])
  } else {
    axis_lim <- raw_range
  }
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .x, y = .y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = REF_LINE_COLOUR, linetype = "dashed") +
    ggplot2::geom_point(
      ggplot2::aes(colour = !is.na(.label)),
      size = 2.2, alpha = 0.8
    ) +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = .label),
      size           = 2,
      colour         = LABEL_COLOUR,
      fill           = "white",
      max.overlaps   = 20,
      na.rm          = TRUE,
      box.padding    = 0.4,
      segment.colour = LABEL_COLOUR,
      segment.size   = 0.3,
      alpha          = 0.5
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "steelblue", "TRUE" = LABEL_COLOUR),
      guide  = "none"
    ) +
    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = title,
      subtitle = subtitle,
      caption  = paste0("Red = top ", n_label,
                        " most extreme species above/below 1:1 line;",
                        " dashed line = 1:1")
    ) +
    theme_comparison()
  
  if (log_axes) {
    p <- p +
      ggplot2::scale_x_continuous(
        trans  = scales::pseudo_log_trans(base = 10),
        breaks = pseudo_log_breaks,
        limits = axis_lim
      ) +
      ggplot2::scale_y_continuous(
        trans  = scales::pseudo_log_trans(base = 10),
        breaks = pseudo_log_breaks,
        limits = axis_lim
      ) +
      ggplot2::coord_equal()
  } else {
    p <- p +
      ggplot2::coord_equal(xlim = axis_lim, ylim = axis_lim)
  }
  
  if (!is.null(summary_text)) {
    # Annotation coordinates are always in transformed space for the label
    # anchor; use -Inf/Inf to stay safely in the plot corners regardless of
    # axis scaling
    p <- p +
      ggplot2::annotate(
        "label",
        x         = -Inf,
        y         =  Inf,
        label     = summary_text,
        hjust     = -0.05,
        vjust     =  1.05,
        size      = 3,
        colour    = "grey20",
        fill      = "white",
        linewidth = 0.3,
        family    = "mono"
      )
  }
  
  p
}

# ------------------------------------------------------------------------------
# Section 7 – Build and save individual scatterplots
# ------------------------------------------------------------------------------

# Filtered data for each panel (mirrors original filter logic)
dat_rho      <- combined %>%
  dplyr::filter(is.finite(rho_bam), is.finite(rho_ebird))

dat_absent   <- combined %>%
  dplyr::filter(is.finite(n_absent_bam), is.finite(n_absent_ebird))

dat_highpred <- combined %>%
  dplyr::filter(is.finite(n_highpred_bam), is.finite(n_highpred_ebird))

dat_propflag <- combined %>%
  dplyr::filter(is.finite(prop_flagged_bam), is.finite(prop_flagged_ebird))

# -- 7a. Spearman correlation --------------------------------------------------
p_rho <- scatter_with_labels(
  data         = dat_rho,
  x            = "rho_bam",
  y            = "rho_ebird",
  x_label      = "BAM v5  —  Spearman rho",
  y_label      = "eBird  —  Spearman rho",
  zdelta_col   = "zdelta_rho",
  title        = "Spearman correlation: BAM v5 vs eBird",
  subtitle     = paste("Predicted abundance vs observed mean count per effort",
                       "(hexagons above absence threshold only)", sep = "\n"),
  summary_text = make_summary_text(dat_rho, "win_rho", nrow(dat_rho))
)

ggplot2::ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_rho.pdf"),
  p_rho, width = 9, height = 8
)
message("Saved: BAM_vs_eBird_rho.pdf")

# -- 7b. Predicted-absent-but-detected flags ----------------------------------
p_absent <- scatter_with_labels(
  data         = dat_absent,
  x            = "n_absent_bam",
  y            = "n_absent_ebird",
  x_label      = "BAM v5  —  n hexagons flagged (pred absent, detected)",
  y_label      = "eBird  —  n hexagons flagged (pred absent, detected)",
  zdelta_col   = "zdelta_absent",
  title        = "Flag count: predicted absent, detected",
  subtitle     = "Higher = more hexagons where model predicts absence but birds were recorded",
  log_axes     = TRUE,
  summary_text = make_summary_text(dat_absent, "win_absent", nrow(dat_absent))
)

ggplot2::ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_absent_flags.pdf"),
  p_absent, width = 9, height = 8
)
message("Saved: BAM_vs_eBird_absent_flags.pdf")

# -- 7c. High-predicted-no-detection flags ------------------------------------
p_highpred <- scatter_with_labels(
  data         = dat_highpred,
  x            = "n_highpred_bam",
  y            = "n_highpred_ebird",
  x_label      = "BAM v5  —  n hexagons flagged (high pred, not detected)",
  y_label      = "eBird  —  n hexagons flagged (high pred, not detected)",
  zdelta_col   = "zdelta_highpred",
  title        = "Flag count: high predicted abundance, not detected",
  subtitle     = "Higher = more hexagons in the top predicted-abundance tail with zero detections",
  log_axes     = TRUE,
  summary_text = make_summary_text(dat_highpred, "win_highpred",
                                   nrow(dat_highpred))
)

ggplot2::ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_highpred_flags.pdf"),
  p_highpred, width = 9, height = 8
)
message("Saved: BAM_vs_eBird_highpred_flags.pdf")

# -- 7d. Proportion of surveyed hexagons flagged ------------------------------
p_propflag <- scatter_with_labels(
  data         = dat_propflag,
  x            = "prop_flagged_bam",
  y            = "prop_flagged_ebird",
  x_label      = "BAM v5  —  proportion of surveyed hexagons flagged",
  y_label      = "eBird  —  proportion of surveyed hexagons flagged",
  zdelta_col   = "zdelta_flag_any",
  title        = "Flag rate: proportion of surveyed hexagons with any flag",
  subtitle     = "Combines both flag types; higher = more inconsistency between model and data",
  summary_text = make_summary_text(dat_propflag, "win_propflag",
                                   nrow(dat_propflag))
)

ggplot2::ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_prop_flagged.pdf"),
  p_propflag, width = 9, height = 8
)
message("Saved: BAM_vs_eBird_prop_flagged.pdf")

# -- Combined 2×2 panel (unchanged layout) ------------------------------------
p_scatter_combined <- (p_rho | p_propflag) / (p_absent | p_highpred) +
  patchwork::plot_annotation(
    title   = "BAM v5 vs eBird — national model assessment comparison",
    caption = paste("Survey data: BAM point counts.",
                    " Both models assessed against identical filtered",
                    " survey dataset."),
    theme   = ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = 15),
      plot.caption = ggplot2::element_text(colour = "grey50", size = 9,
                                           hjust = 0)
    )
  )

ggplot2::ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_scatterplots.pdf"),
  p_scatter_combined, width = 18, height = 16
)
message("Saved: BAM_vs_eBird_scatterplots.pdf")

# -- 7e. Delta-rho histogram --------------------------------------------------
n_bam_rho   <- sum(combined$win_rho >  0, na.rm = TRUE)
n_ebird_rho <- sum(combined$win_rho <  0, na.rm = TRUE)
n_tied_rho  <- sum(combined$win_rho == 0, na.rm = TRUE)

hist_summary <- paste0(
  "BAM better:   ", n_bam_rho,   "\n",
  "eBird better: ", n_ebird_rho, "\n",
  "Tied:         ", n_tied_rho
)

p_delta_rho <- combined %>%
  dplyr::filter(is.finite(delta_rho)) %>%
  dplyr::mutate(favours = dplyr::if_else(delta_rho > 0, "BAM v5", "eBird")) %>%
  ggplot2::ggplot(ggplot2::aes(x = delta_rho, fill = favours)) +
  ggplot2::geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                      colour = REF_LINE_COLOUR) +
  ggplot2::annotate(
    "label",
    x = Inf, y = Inf,
    label      = hist_summary,
    hjust      = 1.05,
    vjust      = 1.05,
    size       = 3,
    colour     = "grey20",
    fill       = "white",
    label.size = 0.3,
    family     = "mono"
  ) +
  ggplot2::scale_fill_manual(
    values = c("BAM v5" = "#2471A3", "eBird" = "#E67E22"),
    name   = "Higher rho"
  ) +
  ggplot2::labs(
    x        = "Spearman rho (BAM) \u2212 Spearman rho (eBird)",
    y        = "Number of species",
    title    = "Distribution of rho differences (BAM \u2212 eBird)",
    subtitle = paste("Positive = BAM has higher rank correlation;",
                     "negative = eBird has higher rank correlation")
  ) +
  theme_comparison()

ggplot2::ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_rho_delta_histogram.pdf"),
  p_delta_rho, width = 8, height = 5
)
message("Saved: BAM_vs_eBird_rho_delta_histogram.pdf")

# ------------------------------------------------------------------------------
# Section 8 – Summary table
# ------------------------------------------------------------------------------

summarise_winner <- function(delta_vec, model_a = "BAM v5",
                             model_b = "eBird", metric_label) {
  tibble::tibble(
    metric       = metric_label,
    n_total      = sum(!is.na(delta_vec)),
    n_model_a    = sum(delta_vec > 0, na.rm = TRUE),
    n_model_b    = sum(delta_vec < 0, na.rm = TRUE),
    n_tied       = sum(delta_vec == 0, na.rm = TRUE),
    median_delta = median(delta_vec, na.rm = TRUE),
    mean_delta   = mean(delta_vec,   na.rm = TRUE)
  ) %>%
    dplyr::rename(
      "{model_a}_better" := n_model_a,
      "{model_b}_better" := n_model_b
    )
}

summary_table <- dplyr::bind_rows(
  summarise_winner(combined$delta_rho,
                   metric_label = "Spearman rho (higher = better)"),
  summarise_winner(-combined$delta_absent,
                   metric_label = "Pred-absent flags (lower = better)"),
  summarise_winner(-combined$delta_highpred,
                   metric_label = "High-pred flags (lower = better)"),
  summarise_winner(-combined$delta_prop_flag,
                   metric_label = "Proportion flagged (lower = better)")
)

print(summary_table, n = Inf)

readr::write_csv(
  summary_table,
  file.path(OUT_COMPARISON, "BAM_vs_eBird_summary_table.csv")
)
message("Saved: BAM_vs_eBird_summary_table.csv")

# ------------------------------------------------------------------------------
# Section 9 – Full comparison table
# ------------------------------------------------------------------------------

readr::write_csv(
  combined %>%
    dplyr::select(
      spcd_bam, sp_english,
      rho_bam, rho_ebird, delta_rho,
      n_absent_bam, n_absent_ebird, delta_absent,
      n_highpred_bam, n_highpred_ebird, delta_highpred,
      prop_flagged_bam, prop_flagged_ebird, delta_prop_flag,
      n_hex_surv_bam,
      overall_score
    ) %>%
    dplyr::arrange(spcd_bam),
  file.path(OUT_COMPARISON, "BAM_vs_eBird_full_comparison.csv")
)
message("Saved: BAM_vs_eBird_full_comparison.csv")

# ------------------------------------------------------------------------------
# Comparison between national hotspots
# ------------------------------------------------------------------------------
species_lookup <- combined %>%
  distinct(spcd_bam, sp_english)

species_to_plot <- combined$spcd_bam %>% unique()

categories_to_plot <- c(
  "Effectively Absent",
  "Very High"
)

hex_summary_dir <- "output/National/hex_summaries"
category_output_dir <- "output/National/comparison/"

overlap_plot_output_dir <- paste0(category_output_dir, "category_overlap_plots/")
side_by_side_output_dir <- paste0(category_output_dir, "side_by_side_category_maps/")

dir.create(overlap_plot_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(side_by_side_output_dir, recursive = TRUE, showWarnings = FALSE)


# ── Compare BAM vs eBird category assignment ─────────────────────────────────

compare_hex_category <- function(hex_summary_bam,
                                 hex_summary_ebird,
                                 target_category) {
  
  ebird_dat <- hex_summary_ebird %>%
    st_drop_geometry() %>%
    select(
      hex_id,
      category_ebird = category
    )
  
  hex_summary_bam %>%
    select(
      hex_id,
      category_bam = category,
      geometry
    ) %>%
    left_join(ebird_dat, by = "hex_id") %>%
    mutate(
      bam_in_category = coalesce(category_bam == target_category, FALSE),
      ebird_in_category = coalesce(category_ebird == target_category, FALSE),
      category_agreement = case_when(
        bam_in_category & ebird_in_category ~ "Both",
        bam_in_category & !ebird_in_category ~ "BAM only",
        !bam_in_category & ebird_in_category ~ "eBird only",
        TRUE ~ NA_character_
      ),
      category_agreement = factor(
        category_agreement,
        levels = c("Both", "BAM only", "eBird only")
      )
    )
}


# ── Collapse categories for side-by-side BAM/eBird maps ──────────────────────

prepare_three_class_map_data <- function(hex_summary) {
  
  hex_summary %>%
    mutate(
      category_3class = case_when(
        category == "Effectively Absent" ~ "Effectively absent",
        category %in% c("Low", "Moderate", "High") ~ "Low / moderate / high",
        category == "Very High" ~ "Very high",
        TRUE ~ NA_character_
      ),
      category_3class = factor(
        category_3class,
        levels = c(
          "Effectively absent",
          "Low / moderate / high",
          "Very high"
        )
      )
    )
}


# ── Plot overlap for one category for one species ────────────────────────────

plot_category_comparison <- function(hex_compare,
                                     region_national,
                                     target_category) {
  
  ggplot() +
    geom_sf(data = region_national, fill = "gray95", color = NA) +
    geom_sf(
      data = hex_compare %>% filter(!is.na(category_agreement)),
      aes(fill = category_agreement),
      color = NA
    ) +
    scale_fill_manual(
      values = c(
        "Both" = "grey20",
        "BAM only" = "#4C78A8",
        "eBird only" = "#E07B73"
      ),
      name = paste0("Hexagons categorized\nas '", target_category, "'"),
      drop = FALSE
    ) +
    ggtitle(target_category) +
    theme_bw() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.8),
        color = "grey70"
      ),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm")
    )
}


# ── Plot one BAM/eBird three-class map ───────────────────────────────────────

plot_three_class_map <- function(hex_summary,
                                 region_national,
                                 plot_title) {
  
  hex_plot <- prepare_three_class_map_data(hex_summary)
  
  ggplot() +
    geom_sf(data = region_national, fill = "gray80", color = "gray60") +
    geom_sf(
      data = hex_plot %>% filter(!is.na(category_3class)),
      aes(fill = category_3class),
      color = NA
    ) +
    scale_fill_manual(
      values = c(
        "Effectively absent" = "white",
        "Low / moderate / high" = "#8FAE8B",
        "Very high" = "black"
      ),
      name = "Predicted category",
      drop = FALSE
    ) +
    geom_sf(data = region_national, fill = NA, color = "gray60", lwd = 0.2) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      panel.background = element_rect(
        fill = "grey98",
        colour = NA
      ),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.8),
        color = "grey70"
      ),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm")
    )
}


# ── Create side-by-side BAM/eBird map for one species ────────────────────────

plot_bam_ebird_side_by_side <- function(hex_summary_bam,
                                        hex_summary_ebird,
                                        region_national,
                                        sp_english) {
  
  p_bam <- plot_three_class_map(
    hex_summary = hex_summary_bam,
    region_national = region_national,
    plot_title = "BAM v5"
  )
  
  p_ebird <- plot_three_class_map(
    hex_summary = hex_summary_ebird,
    region_national = region_national,
    plot_title = "eBird"
  )
  
  p_bam + p_ebird + plot_annotation(title = sp_english)
    # plot_layout(guides = "collect") +
    # theme(
    #   legend.position = "bottom"
    # )
}


# ── Summarize one category for one species ───────────────────────────────────

summarize_category_comparison <- function(hex_compare,
                                          spcd_bam,
                                          target_category) {
  
  hex_compare_area <- hex_compare %>%
    mutate(
      area_km2 = as.numeric(set_units(st_area(.), "km^2")),
      both_in_category = bam_in_category & ebird_in_category,
      disagreement = xor(bam_in_category, ebird_in_category)
    )
  
  total_area_bam <- sum(hex_compare_area$area_km2[hex_compare_area$bam_in_category], na.rm = TRUE)
  total_area_ebird <- sum(hex_compare_area$area_km2[hex_compare_area$ebird_in_category], na.rm = TRUE)
  overlap_area <- sum(hex_compare_area$area_km2[hex_compare_area$both_in_category], na.rm = TRUE)
  disagreement_area <- sum(hex_compare_area$area_km2[hex_compare_area$disagreement], na.rm = TRUE)
  
  tibble(
    species = spcd_bam,
    category = target_category,
    total_area_bam_km2 = total_area_bam,
    total_area_ebird_km2 = total_area_ebird,
    overlap_area_km2 = overlap_area,
    disagreement_area_km2 = disagreement_area,
    overlap_to_disagreement_ratio = if_else(
      disagreement_area > 0,
      overlap_area / disagreement_area,
      NA_real_
    )
  )
}


# ── Process one species ──────────────────────────────────────────────────────

process_species <- function(spcd_bam,
                            sp_english,
                            categories_to_plot,
                            region_national,
                            hex_summary_dir,
                            overlap_plot_output_dir,
                            side_by_side_output_dir) {
  
  message("Processing ", spcd_bam)
  
  hex_summary_bam <- readRDS(
    file.path(hex_summary_dir, paste0(spcd_bam, "_BAM_hex_summary.rds"))
  )
  
  hex_summary_ebird <- readRDS(
    file.path(hex_summary_dir, paste0(spcd_bam, "_eBird_hex_summary.rds"))
  )
  
  # ── Existing overlap plots and stats ───────────────────────────────────────
  
  category_results <- map(categories_to_plot, function(target_category) {
    
    hex_compare <- compare_hex_category(
      hex_summary_bam = hex_summary_bam,
      hex_summary_ebird = hex_summary_ebird,
      target_category = target_category
    )
    
    p <- plot_category_comparison(
      hex_compare = hex_compare,
      region_national = region_national,
      target_category = target_category
    )
    
    stats <- summarize_category_comparison(
      hex_compare = hex_compare,
      spcd_bam = spcd_bam,
      target_category = target_category
    )
    
    list(
      plot = p,
      stats = stats
    )
  })
  
  overlap_plot <- wrap_plots(
    map(category_results, "plot"),
    ncol = length(categories_to_plot)
  ) +
    plot_annotation(title = paste0(spcd_bam, " — category overlap"))
  
  ggsave(
    filename = file.path(
      overlap_plot_output_dir,
      paste0(spcd_bam, "_category_overlap.png")
    ),
    plot = overlap_plot,
    width = 14,
    height = 7,
    dpi = 300
  )
  
  # ── New side-by-side BAM/eBird maps ────────────────────────────────────────
  
  side_by_side_plot <- plot_bam_ebird_side_by_side(
    hex_summary_bam = hex_summary_bam,
    hex_summary_ebird = hex_summary_ebird,
    region_national = region_national,
    sp_english = sp_english
  )
  
  ggsave(
    filename = file.path(
      side_by_side_output_dir,
      paste0(spcd_bam, "_BAM_eBird_side_by_side_categories.png")
    ),
    plot = side_by_side_plot,
    width = 14,
    height = 7,
    dpi = 300
  )
  
  bind_rows(map(category_results, "stats"))
}


# ── Run all species ──────────────────────────────────────────────────────────

comparison_stats_all <- map_dfr(
  species_to_plot,
  function(spcd_bam) {
    
    sp_english <- species_lookup %>%
      filter(spcd_bam == !!spcd_bam) %>%
      pull(sp_english) %>%
      first()
    
    process_species(
      spcd_bam = spcd_bam,
      sp_english = sp_english,
      categories_to_plot = categories_to_plot,
      region_national = region_national,
      hex_summary_dir = hex_summary_dir,
      overlap_plot_output_dir = overlap_plot_output_dir,
      side_by_side_output_dir = side_by_side_output_dir
    )
  }
)

comparison_stats_all

write_csv(
  comparison_stats_all,
  file.path(category_output_dir, "category_comparison_summary_stats.csv")
)


# ── Summarize across species ─────────────────────────────────────────────────

comparison_stats_summary <- comparison_stats_all %>%
  group_by(category) %>%
  summarise(
    n_species = n(),
    median_overlap_to_disagreement = median(
      overlap_to_disagreement_ratio,
      na.rm = TRUE
    ),
    prop_bam_larger = mean(
      total_area_bam_km2 > total_area_ebird_km2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

comparison_stats_summary

# ------------
# Compare total area occupied
summarize_non_absent_area <- function(spcd_bam, hex_summary_dir) {
  
  hex_summary_bam <- readRDS(
    file.path(hex_summary_dir, paste0(spcd_bam, "_BAM_hex_summary.rds"))
  )
  
  hex_summary_ebird <- readRDS(
    file.path(hex_summary_dir, paste0(spcd_bam, "_eBird_hex_summary.rds"))
  )
  
  bam_area <- hex_summary_bam %>%
    filter(!is.na(category), category != "Effectively Absent") %>%
    mutate(area_km2 = as.numeric(set_units(st_area(.), "km^2"))) %>%
    summarise(total_area_bam_km2 = sum(area_km2, na.rm = TRUE)) %>%
    pull(total_area_bam_km2)
  
  ebird_area <- hex_summary_ebird %>%
    filter(!is.na(category), category != "Effectively Absent") %>%
    mutate(area_km2 = as.numeric(set_units(st_area(.), "km^2"))) %>%
    summarise(total_area_ebird_km2 = sum(area_km2, na.rm = TRUE)) %>%
    pull(total_area_ebird_km2)
  
  tibble(
    species = spcd_bam,
    total_area_bam_km2 = bam_area,
    total_area_ebird_km2 = ebird_area,
    bam_larger = total_area_bam_km2 > total_area_ebird_km2
  )
}

total_area_occupied <- map_dfr(
  species_to_plot,
  summarize_non_absent_area,
  hex_summary_dir = hex_summary_dir
)

total_area_occupied

# ------------------------------------------------------------------------------
# Section X – Species-level summed predicted abundance from each program
# ------------------------------------------------------------------------------

summarise_pred_total <- function(spcd_bam, model, hex_summary_dir) {
  
  path <- file.path(
    hex_summary_dir,
    paste0(spcd_bam, "_", model, "_hex_summary.rds")
  )
  
  if (!file.exists(path)) {
    warning("Missing hex summary: ", path)
    return(tibble(
      spcd_bam = spcd_bam,
      model = model,
      pred_sum = NA_real_,
      n_hex_valid = NA_integer_
    ))
  }
  
  hex_summary <- readRDS(path)
  
  tibble(
    spcd_bam = spcd_bam,
    model = model,
    pred_sum = sum(hex_summary$pred_population_index, na.rm = TRUE),
    n_hex_valid = sum(is.finite(hex_summary$pred_population_index))
  )
}

pred_totals <- tidyr::expand_grid(
  spcd_bam = combined$spcd_bam,
  model = c("BAM", "eBird")
) %>%
  purrr::pmap_dfr(
    ~ summarise_pred_total(
      spcd_bam = ..1,
      model = ..2,
      hex_summary_dir = hex_summary_dir
    )
  ) %>%
  tidyr::pivot_wider(
    names_from = model,
    values_from = c(pred_sum, n_hex_valid),
    names_glue = "{.value}_{model}"
  )


# ------------------------------------------------------------------------------
# Summed relative abundance: BAM vs eBird (how correlated are they?)
# ------------------------------------------------------------------------------

dat_pred_sum <- combined %>%
  filter(
    is.finite(pred_sum_BAM),
    is.finite(pred_sum_eBird)
  ) %>%
  mutate(
    x_log = log10(pred_sum_BAM + 1),
    y_log = log10(pred_sum_eBird + 1)
  )

mod <- lm(y_log ~ x_log, data = dat_pred_sum)

dat_pred_sum <- dat_pred_sum %>%
  mutate(
    residual = residuals(mod),
    abs_residual = abs(residual),
    label = if_else(
      rank(-abs_residual, ties.method = "first") <= 20,
      sp_english,
      NA_character_
    )
  )

p_pred_sum <- ggplot(
  dat_pred_sum,
  aes(x = pred_sum_BAM, y = pred_sum_eBird)
) +
  geom_point(alpha = 0.8, size = 2.2) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_label_repel(
    aes(label = label),
    size = 2.5,
    na.rm = TRUE,
    max.overlaps = Inf
  ) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(
    x = "BAM v5 — summed predicted relative abundance",
    y = "eBird — summed predicted relative abundance",
    title = "Species-level summed predicted abundance: BAM v5 vs eBird",
    subtitle = "Each point is one species; axes use pseudo-log scaling"
  ) +
  theme_comparison()

ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_summed_relative_abundance.pdf"),
  p_pred_sum, width = 8, height = 7
)

# ------------------------------------------------------------------------------
# Summed relative abundance: BAM vs eBird (how correlated are they?)
# ------------------------------------------------------------------------------

dat_pred_sum <- combined %>%
  filter(
    is.finite(pred_sum_BAM),
    is.finite(pred_sum_eBird)
  ) %>%
  mutate(
    x_log = log10(pred_sum_BAM + 1),
    y_log = log10(pred_sum_eBird + 1)
  )

mod <- lm(y_log ~ x_log, data = dat_pred_sum)

dat_pred_sum <- dat_pred_sum %>%
  mutate(
    residual = residuals(mod),
    abs_residual = abs(residual),
    label = if_else(
      rank(-abs_residual, ties.method = "first") <= 20,
      sp_english,
      NA_character_
    )
  )

p_pred_sum <- ggplot(
  dat_pred_sum,
  aes(x = pred_sum_BAM, y = pred_sum_eBird)
) +
  geom_point(alpha = 0.8, size = 2.2) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_label_repel(
    aes(label = label),
    size = 2.5,
    na.rm = TRUE,
    max.overlaps = Inf
  ) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(
    x = "BAM v5 — summed predicted relative abundance",
    y = "eBird — summed predicted relative abundance",
    title = "Species-level summed predicted abundance: BAM v5 vs eBird",
    subtitle = "Each point is one species; axes use pseudo-log scaling"
  ) +
  theme_comparison()

ggsave(
  file.path(OUT_COMPARISON, "BAM_vs_eBird_summed_relative_abundance.pdf"),
  p_pred_sum, width = 8, height = 7
)

# ------------------------------------------------------------------------------
# Are flag counts related to species rarity?
# ------------------------------------------------------------------------------

flag_rarity_long <- combined %>%
  select(
    spcd_bam, sp_english,
    pred_sum_BAM,
    n_absent_bam, n_highpred_bam, total_flags_bam,
    n_absent_ebird, n_highpred_ebird, total_flags_ebird
  ) %>%
  pivot_longer(
    cols = c(
      n_absent_bam, n_highpred_bam, total_flags_bam,
      n_absent_ebird, n_highpred_ebird, total_flags_ebird
    ),
    names_to = "flag_metric",
    values_to = "n_flags"
  ) %>%
  mutate(
    model = case_when(
      stringr::str_ends(flag_metric, "_bam") ~ "BAM v5",
      stringr::str_ends(flag_metric, "_ebird") ~ "eBird"
    ),
    flag_type = case_when(
      stringr::str_detect(flag_metric, "absent") ~ "Predicted absent, detected",
      stringr::str_detect(flag_metric, "highpred") ~ "High predicted, not detected",
      stringr::str_detect(flag_metric, "total") ~ "Total flags"
    ),
    flag_type = factor(
      flag_type,
      levels = c(
        "Predicted absent, detected",
        "High predicted, not detected",
        "Total flags"
      )
    )
  ) %>%
  filter(is.finite(pred_sum_BAM), is.finite(n_flags))

p_flags_vs_rarity <- ggplot(
  flag_rarity_long,
  aes(x = pred_sum_BAM, y = n_flags)
) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), se = TRUE) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  facet_grid(flag_type ~ model, scales = "free_y") +
  labs(
    x = "BAM v5 summed predicted relative abundance",
    y = "Number of flagged hexagons",
    title = "Are flag counts related to species rarity?",
    subtitle = "Lower BAM summed abundance = rarer species in the BAM prediction surface"
  ) +
  theme_comparison()
p_flags_vs_rarity


ggsave(
  file.path(OUT_COMPARISON, "flags_vs_BAM_summed_relative_abundance.pdf"),
  p_flags_vs_rarity, width = 11, height = 9
)