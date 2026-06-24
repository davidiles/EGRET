# ==============================================================================
# assess_BAM_national.R
#
# National-scale assessment of BAM v5 model predictions against BAM
# point-count survey observations.
#
# For each species with at least MIN_DETECTIONS detections in the BAM dataset,
# this script:
#   1. Loads the species-specific BAM v5 raster (mean layer).
#   2. Derives a breeding-season date window from eBird seasonal metadata
#      (±DATE_BUFFER_DAYS around the published breeding start/end dates).
#   3. Filters the BAM survey data to that date window.
#   4. Summarises observations and predictions to the national hex grid.
#   5. Computes species-specific cumulative-population thresholds via
#      get_thresholds_cumpop() and attaches a five-category "category" column
#      to the hex summary via apply_thresholds():
#        - Effectively Absent: hexagons summing to the lowest CUMPOP_BREAKS$absent_upper
#          fraction of the national predicted population.
#        - Low / Moderate / High / Very High: successive population quartiles.
#      The absence threshold (absent_upper) and high-abundance threshold
#      (high_upper) drive the diagnostic flags in assess_region().
#   6. Saves a multi-panel PDF figure, the region-wide statistics, and
#      the full hex summary sf object (for use in region-specific analyses
#      without re-running this script).
#
# All shared constants and file paths live in config.R.
# All analysis functions live in functions_optimized.R.
# ==============================================================================


# ------------------------------------------------------------------------------
# Section 1 – Libraries
# ------------------------------------------------------------------------------

library(BAMexploreR)
library(tidyverse)
library(terra)
library(viridis)
library(rnaturalearth)
library(patchwork)
library(sf)
library(lwgeom)
library(ggnewscale)
library(lubridate)
library(exactextractr)
library(ebirdst)


# ------------------------------------------------------------------------------
# Section 2 – Configuration and functions
# ------------------------------------------------------------------------------

source("R/config.R")
source("R/functions.R")

# ------------------------------------------------------------------------------
# Section 3 – Output directories
# ------------------------------------------------------------------------------

dir.create(OUT_NATIONAL_BAM_FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_NATIONAL_BAM_STATS,   recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_HEX_SUMMARIES,        recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_CACHE,                recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# Section 4 – Study area: Canadian provinces intersected with BCRs
#
# The national study area is Canada's provincial boundaries dissolved and
# clipped to Bird Conservation Regions, excluding BCR_EXCLUDE units (defined
# in config.R) that fall outside the boreal/temperate zone modelled by BAM.
# Interior holes and very small fragments are removed for a clean boundary.
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# Section 5 – Load and filter BAM survey data
#
# Surveys outside the acceptable duration and distance ranges are dropped here,
# once, before any species-specific processing.  The date-window filter is
# applied per species inside the loop (Section 9) because each species has a
# different breeding season.
# ------------------------------------------------------------------------------

load(BAM_DATA_PATH)   # loads object 'dat'

bam_sf <- dat %>%
  dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
  dplyr::mutate(
    survey_datetime = lubridate::ymd_hms(date_time, tz = "UTC"),
    day_of_year     = lubridate::yday(survey_datetime),
    year            = lubridate::year(survey_datetime)
  ) %>%
  dplyr::filter(
    duration >= SURVEY_DURATION_MIN_S,
    duration <= SURVEY_DURATION_MAX_S,
    distance >  SURVEY_DISTANCE_MIN_M
  ) %>%
  sf::st_as_sf(
    coords = c("longitude", "latitude"),
    crs    = 4326,
    remove = FALSE
  ) %>%
  sf::st_transform(sf::st_crs(region_national)) %>%
  sf::st_filter(region_national, .predicate = sf::st_intersects)

rm(dat)   # free memory; bam_sf is the canonical survey object from here on


# ------------------------------------------------------------------------------
# Section 6 – Species list
#
# Retain only species columns (ABBE:YTWA) with at least MIN_DETECTIONS
# detections across all surveys in the study area.  MIN_DETECTIONS is defined
# in config.R and applied consistently across both assessment scripts.
# ------------------------------------------------------------------------------

n_detections <- bam_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(dplyr::across(ABBE:YTWA, ~ sum(.x > 0, na.rm = TRUE)))

species_vec   <- names(n_detections)[n_detections >= MIN_DETECTIONS]
species_table <- BAMexploreR::spp_tbl
ebirdst_runs  <- ebirdst::ebirdst_runs


# ------------------------------------------------------------------------------
# Section 7 – Hex grid (create or load from cache)
#
# The hex grid is shared with assess_ebird_national.R.  If that script has
# already run, both cache files will already exist and will be loaded directly.
# The cache is keyed only on existence; the grid geometry is stable as long as
# HEX_WIDTH_KM and region_national do not change.
# ------------------------------------------------------------------------------

hex_grid_path <- file.path(OUT_CACHE, "hex_grid_national_10_km.rds")

hex_grid_national <- load_or_compute(
  cache_path = hex_grid_path,
  compute_fn = function() make_hex_grid(region_national,
                                        hex_width_km = HEX_WIDTH_KM)
)


# ------------------------------------------------------------------------------
# Section 8 – Precomputed hex metadata and survey–hex spatial join
#
# PRECOMPUTATION 1: hex area, inner diameter, and study-area coverage fraction.
# PRECOMPUTATION 2: assign each survey its hexagon ID.
#
# Both are expensive to compute over 80,000+ hexagons and are therefore cached.
# The survey–hex cache is invalidated if BAM_DATA_PATH is newer than the cache
# file, ensuring stale spatial joins are not silently reused after a data update.
# ------------------------------------------------------------------------------

hex_meta_path <- file.path(OUT_CACHE, "hex_grid_national_10_km_meta.rds")

hex_grid_national <- load_or_compute(
  cache_path = hex_meta_path,
  compute_fn = function() precompute_hex_metadata(
    hex_grid   = hex_grid_national,
    study_area = region_national
  )
)

bam_hex_path <- file.path(OUT_CACHE, "bam_sf_with_hex_id.rds")

bam_sf <- load_or_compute(
  cache_path  = bam_hex_path,
  compute_fn  = function() precompute_survey_hex_id(
    survey_sf = bam_sf,
    hex_grid  = hex_grid_national
  ),
  source_path = BAM_DATA_PATH   # invalidate cache if source data is updated
)


# ------------------------------------------------------------------------------
# Section 9 – Resume logic
#
# Results are appended to an on-disk RDS after each species.  If the script is
# interrupted and restarted, already-completed species are skipped.
# ------------------------------------------------------------------------------

stats_path <- file.path(
  OUT_NATIONAL_BAM_STATS, "BAM_preds_vs_BAM_data.rds"
)

national_stats_df <- if (file.exists(stats_path)) {
  readRDS(stats_path)
} else {
  tibble::tibble()
}

completed_species <- if (nrow(national_stats_df) > 0) {
  national_stats_df %>%
    dplyr::filter(region_id == "national", model == "BAM v5") %>%
    dplyr::pull(spcd_bam) %>%
    unique()
} else {
  character()
}


# ------------------------------------------------------------------------------
# Section 10 – Main species loop
# ------------------------------------------------------------------------------

results_list <- vector("list", length(species_vec))
names(results_list) <- species_vec

for (spcd in species_vec) {

  # -- 10a. Skip if already complete ------------------------------------------
  if (spcd %in% completed_species) {
    message("Skipping ", spcd, " — already in national_stats_df")
    next
  }

  # -- 10b. Resolve common name -----------------------------------------------
  sp_english <- species_table %>%
    dplyr::filter(speciesCode == spcd) %>%
    dplyr::pull(commonName) %>%
    dplyr::first()

  if (length(sp_english) == 0 || is.na(sp_english)) {
    message("Skipping ", spcd, " — no common name found")
    next
  }

  # -- 10c. Check raster exists -----------------------------------------------
  bam_raster_path_sp <- file.path(
    BAM_RASTER_PATH,
    paste0(spcd, "_Canada_2020.tif")
  )

  if (!file.exists(bam_raster_path_sp)) {
    message("Skipping ", spcd, " — raster not found at ", bam_raster_path_sp)
    next
  }

  # -- 10d. Breeding-season date window ---------------------------------------
  # get_species_date_window() is defined in functions_optimized.R and is called
  # identically in assess_ebird_national.R, ensuring consistent survey filtering
  # across both model assessments.
  date_window <- get_species_date_window(
    sp_english     = sp_english,
    ebirdst_runs   = ebirdst_runs,
    buffer_days    = DATE_BUFFER_DAYS,
    fallback_start = FALLBACK_START_DATE,
    fallback_end   = FALLBACK_END_DATE
  )

  message("Processing ", spcd, " — ", sp_english,
          " (", date_window$date_label, ")")
  start_time <- Sys.time()

  # -- 10e. Species assessment ------------------------------------------------
  result <- tryCatch({

    # Load BAM raster (mean layer only)
    bam_rast <- terra::rast(bam_raster_path_sp)[["mean"]]

    # Filter surveys to the species-specific date window
    sp_dat <- make_sp_dat(bam_sf, spcd) %>%
      dplyr::mutate(
        day_of_year = lubridate::yday(survey_datetime)
      ) %>%
      dplyr::filter(
        !is.na(day_of_year),
        day_of_year >= date_window$start_doy,
        day_of_year <= date_window$end_doy,
        !is.na(count),
        !is.na(count_per_effort)
      )

    # Compute species-specific cumulative-population thresholds.
    # summarize_hex() is run first so pred_mean is available; the resulting
    # hex_summary is reused by assess_region() via hex_summary_precomputed to
    # avoid a second raster extraction.
    hex_summary <- summarize_hex(
      dat      = sp_dat,
      hex_grid = hex_grid_national,
      rast     = bam_rast,
      min_coverage = 0.25
    )

    # Derive five-category thresholds from the predicted-abundance distribution.
    # CUMPOP_BREAKS is defined in config.R and controls all five breakpoints.
    thresholds <- get_thresholds_cumpop(
      hex_summary   = hex_summary,
      cumpop_breaks = CUMPOP_BREAKS
    )

    # Attach the category column to the hex summary so it is saved to disk
    # and available to downstream region-specific scripts.
    hex_summary <- apply_thresholds(hex_summary, thresholds)

    # Run the assessment, passing pre-computed absolute thresholds.
    #   flag_absence_threshold  <- thresholds$absent_upper (Effectively Absent)
    #   flag_high_pred_threshold <- thresholds$high_upper  (top of "High" tier)
    # flag_high_pred_cumprop and flag_absence_fraction_of_high are set to NULL
    # to disable the internal cumulative calculations inside assess_region().
    assess_region(
      flag_high_pred_cumprop        = NULL,
      flag_high_pred_threshold      = thresholds$mod_upper,
      flag_absence_fraction_of_high = NULL,
      flag_absence_threshold        = thresholds$absent_upper,
      flag_high_pred_min_surveys    = FLAG_HIGH_PRED_MIN_SURVEYS,
      flag_absence_min_surveys      = FLAG_ABSENCE_MIN_SURVEYS,

      region         = region_national,
      sp_dat         = sp_dat,
      rast           = bam_rast,
      hex_grid       = hex_grid_national,
      hex_summary_precomputed = hex_summary,
      water          = NULL,
      transform      = TRANSFORM,
      flag_linewidth = FLAG_LINEWIDTH,
      title          = paste(
        sp_english,
        "- 2020 National scale\nData = BAM | Model = BAM v5\nData filtered:",
        date_window$date_label
      )
    )

  }, error = function(e) {
    message("Failed ", spcd, " — ", sp_english, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(result)) next

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  message(sprintf("  Done in %.1f min", as.numeric(elapsed)))

  # -- 10f. Save figure -------------------------------------------------------
  pdf_path <- file.path(
    OUT_NATIONAL_BAM_FIGURES,
    paste0(spcd, "_BAM_national.pdf")
  )
  grDevices::pdf(pdf_path, width = 30, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(result$plot_combined)
  grDevices::dev.off()
  on.exit(NULL)   # clear the on.exit handler after successful dev.off()

  # -- 10g. Save hex summary sf object ----------------------------------------
  # The full hex summary (geometry + all columns) is saved per species so that
  # region-specific scripts can subset it spatially without re-running this
  # national analysis.
  hex_summary_path <- file.path(
    OUT_HEX_SUMMARIES,
    paste0(spcd, "_BAM_hex_summary.rds")
  )
  saveRDS(result$hex_summary, hex_summary_path)

  # -- 10h. Accumulate statistics ---------------------------------------------
  new_stats <- result$region_stats %>%
    dplyr::mutate(
      spcd_bam   = spcd,
      sp_english = sp_english,
      model      = "BAM v5",
      data       = "BAM v5",
      region_id  = "national"
    )

  results_list[[spcd]] <- new_stats

  # Append to the on-disk file after each species so partial results survive
  # an interruption
  national_stats_df <- dplyr::bind_rows(
    national_stats_df,
    new_stats
  )
  saveRDS(national_stats_df, stats_path)
  completed_species <- c(completed_species, spcd)
}


# ------------------------------------------------------------------------------
# Section 11 – Summary tables
# ------------------------------------------------------------------------------

national_stats_df %>%
  dplyr::relocate(sp_english) %>%
  dplyr::arrange(spearman_obs_pred_above_absence)

national_stats_df %>%
  dplyr::relocate(sp_english) %>%
  dplyr::arrange(dplyr::desc(n_flag_pred_absent_detected))

national_stats_df %>%
  dplyr::relocate(sp_english) %>%
  dplyr::arrange(dplyr::desc(n_flag_high_pred_no_detection))
