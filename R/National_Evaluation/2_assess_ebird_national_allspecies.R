# ==============================================================================
# assess_ebird_national.R
#
# National-scale assessment of eBird Status & Trends seasonal abundance
# predictions against BAM point-count survey observations.
#
# For each species available in both the BAM dataset (≥ MIN_DETECTIONS) and
# eBird Status & Trends, this script:
#   1. Loads the species-specific eBird seasonal raster (breeding or resident).
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
# All shared constants and file paths live in config_new.R.
# All analysis functions live in functions_optimized.R.
# ==============================================================================


# ------------------------------------------------------------------------------
# Section 1 – Libraries
# ------------------------------------------------------------------------------

library(BAMexploreR)
library(tidyverse)
library(terra)
library(viridis)
library(ebirdst)
library(rnaturalearth)
library(patchwork)
library(sf)
library(lwgeom)
library(ggnewscale)
library(lubridate)
library(exactextractr)
library(readxl)


# ------------------------------------------------------------------------------
# Section 2 – Configuration and functions
# ------------------------------------------------------------------------------

source("R/config.R")
source("R/functions.R")


# ------------------------------------------------------------------------------
# Section 3 – Output directories
# ------------------------------------------------------------------------------

dir.create(OUT_NATIONAL_EBIRD_FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_NATIONAL_EBIRD_STATS,   recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_HEX_SUMMARIES,          recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_CACHE,                  recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Section 4 – Study area: Canadian provinces intersected with BCRs
#
# Identical to assess_BAM_national.R.  The BCR exclusions and hole-removal
# step are controlled by BCR_EXCLUDE (config.R) and remove_multipolygon_holes()
# (functions_optimized.R) so any change propagates to both scripts.
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
# Section 6 – Species list
#
# Consistent with assess_BAM_national.R: retain species with at least
# MIN_DETECTIONS detections in the BAM dataset (defined in config.R).
# The eBird availability check happens per species inside the loop.
# ------------------------------------------------------------------------------

# load avian core
avian_core <- read_xlsx(AVIAN_CORE_PATH) %>%
  filter(Taxon == "Aves",
         Full_Species__Espèce_complète == "Yes - Oui",
         CDN_Status__Statut_CDN %in% c("BRE","BRE_OCC","RNB"))

ebirdst_runs  <- ebirdst::ebirdst_runs %>%
  filter((common_name %in% avian_core$English_Name__Nom_Anglais |
            scientific_name %in% avian_core$Scientific_Name__Nom_Scientifique) &
           species_code != "yebsap-example") %>%
  mutate(ebird_common_name = common_name,
         ebird_species_code = species_code)

# ------------------------------------------------------------------------------
# Section 7 – Hex grid (create or load from cache)
#
# Shared cache with assess_BAM_national.R.  If BAM has already run the grid
# will be loaded directly.
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
# Shared cache with assess_BAM_national.R.  The survey–hex cache is
# invalidated automatically if BAM_DATA_PATH is newer than the cache.
# ------------------------------------------------------------------------------

hex_meta_path <- file.path(OUT_CACHE, "hex_grid_national_10_km_meta.rds")

hex_grid_national <- load_or_compute(
  cache_path = hex_meta_path,
  compute_fn = function() precompute_hex_metadata(
    hex_grid   = hex_grid_national,
    study_area = region_national
  )
)

# bam_hex_path <- file.path(OUT_CACHE, "bam_sf_with_hex_id.rds")
# 
# bam_sf <- load_or_compute(
#   cache_path  = bam_hex_path,
#   compute_fn  = function() precompute_survey_hex_id(
#     survey_sf = bam_sf,
#     hex_grid  = hex_grid_national
#   ),
#   source_path = BAM_DATA_PATH   # invalidate if source data is updated
# )


# ------------------------------------------------------------------------------
# Section 9 – Resume logic
# ------------------------------------------------------------------------------

stats_path <- file.path(
  OUT_NATIONAL_EBIRD_STATS, "eBird_preds.rds"
)

national_stats_df <- if (file.exists(stats_path)) {
  readRDS(stats_path)
} else {
  tibble::tibble()
}

completed_species <- if (nrow(national_stats_df) > 0) {
  national_stats_df %>%
    dplyr::filter(region_id == "national", model == "eBird") %>%
    dplyr::pull(spcd_bam) %>%
    unique()
} else {
  character()
}


# ------------------------------------------------------------------------------
# Section 10 – Main species loop
# ------------------------------------------------------------------------------

species_table <- avian_core %>%
  dplyr::select(English_Name__Nom_Anglais,Scientific_Name__Nom_Scientifique,Alpha_Code_BBL__Code_Alpha_BBL) %>%
  dplyr::rename(english_name = English_Name__Nom_Anglais,
                scientific_name = Scientific_Name__Nom_Scientifique,
                bbl_code = Alpha_Code_BBL__Code_Alpha_BBL) %>%
  left_join(ebirdst_runs) %>%
  
  # Removes 8 species
  filter(!is.na(ebird_species_code))

species_vec <- species_table$ebird_species_code

results_list <- vector("list", length(species_vec))
names(results_list) <- species_vec

for (spcd_ebird in species_vec) {
  
  # -- 10a. Skip if already complete ------------------------------------------
  if (spcd_ebird %in% completed_species) {
    message("Skipping ", spcd_ebird, " — already in national_stats_df")
    next
  }
  
  # -- 10b. Resolve common name -----------------------------------------------
  sp_english <- species_table %>%
    dplyr::filter(ebird_species_code == spcd_ebird) %>%
    dplyr::pull(english_name) %>%
    dplyr::first()
  
  if (length(sp_english) == 0 || is.na(sp_english)) {
    message("Skipping ", spcd_ebird, " — no common name found")
    next
  }
  
  # -- 10d. Breeding-season date window ---------------------------------------
  # get_species_date_window() is defined in functions_optimized.R and is called
  # identically in assess_BAM_national.R, ensuring consistent survey filtering
  # across both model assessments.
  date_window <- get_species_date_window(
    sp_english     = sp_english,
    ebirdst_runs   = ebirdst_runs,
    buffer_days    = DATE_BUFFER_DAYS,
    fallback_start = FALLBACK_START_DATE,
    fallback_end   = FALLBACK_END_DATE
  )
  
  message("Processing ", sp_english,
          " — eBird: ", spcd_ebird,
          " (", date_window$date_label, ")")
  start_time <- Sys.time()
  
  # -- 10e. Species assessment ------------------------------------------------
  result <- tryCatch({
    
    # Load eBird seasonal raster
    rast_ebird <- ebirdst::load_raster(
      spcd_ebird,
      product    = "abundance",
      period     = "seasonal",
      resolution = "3km"
    )
    
    # Select the appropriate seasonal layer.
    # Resident species use a combined resident layer; migratory species use
    # the breeding layer.  If neither is present, skip the species.
    if ("resident" %in% names(rast_ebird)) {
      rast_ebird <- rast_ebird[["resident"]]
    } else if ("breeding" %in% names(rast_ebird)) {
      rast_ebird <- rast_ebird[["breeding"]]
    } else {
      stop("No resident or breeding layer found in eBird raster")
    }
    
    # # Filter surveys to the species-specific date window
    # sp_dat <- make_sp_dat(bam_sf, spcd_ebird) %>%
    #   dplyr::mutate(
    #     day_of_year = lubridate::yday(survey_datetime)
    #   ) %>%
    #   dplyr::filter(
    #     !is.na(day_of_year),
    #     day_of_year >= date_window$start_doy,
    #     day_of_year <= date_window$end_doy,
    #     !is.na(count),
    #     !is.na(count_per_effort)
    #   )
    
    # Compute species-specific cumulative-population thresholds.
    # summarize_hex() is run first so pred_mean is available; the resulting
    # hex_summary is reused by assess_region() via hex_summary_precomputed to
    # avoid a second raster extraction.
    #
    # eBird-specific subtlety preserved:
    #   Thresholds are computed from the selected seasonal layer (resident if
    #   available, otherwise breeding), not from all eBird seasonal layers.
    
    hex_summary <- extract_hex_predictions(
      hex_grid_national,
      rast_ebird,
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
    #   flag_absence_threshold   <- thresholds$absent_upper (Effectively Absent)
    #   flag_high_pred_threshold <- thresholds$high_upper   (top of "High" tier)
    #   flag_high_pred_cumprop and flag_absence_fraction_of_high are set to NULL
    #   to disable the internal cumulative calculations inside assess_region().
    
    # assess_region(
    #   flag_high_pred_cumprop        = NULL,
    #   flag_high_pred_threshold      = thresholds$mod_upper,
    #   flag_absence_fraction_of_high = NULL,
    #   flag_absence_threshold        = thresholds$absent_upper,
    #   flag_high_pred_min_surveys    = FLAG_HIGH_PRED_MIN_SURVEYS,
    #   flag_absence_min_surveys      = FLAG_ABSENCE_MIN_SURVEYS,
    #   
    #   region         = region_national,
    #   sp_dat         = sp_dat,
    #   rast           = rast_ebird,
    #   hex_grid       = hex_grid_national,
    #   hex_summary_precomputed = hex_summary,
    #   water          = NULL,
    #   transform      = TRANSFORM,
    #   flag_linewidth = FLAG_LINEWIDTH,
    #   title          = paste(
    #     sp_english,
    #     "- 2020 National scale\nData = BAM | Model = eBird\nData filtered:",
    #     date_window$date_label
    #   )
    # )
    
  }, error = function(e) {
    message("Failed ", spcd_ebird, " — ", sp_english, ": ", conditionMessage(e))
    NULL
  })
  
  if (is.null(result)) next
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  message(sprintf("  Done in %.1f min", as.numeric(elapsed)))
  
  # -- 10f. Save figure -------------------------------------------------------
  # pdf_path <- file.path(
  #   OUT_NATIONAL_EBIRD_FIGURES,
  #   paste0(spcd_ebird, "_eBird_national.pdf")
  # )
  # grDevices::pdf(pdf_path, width = 30, height = 10)
  # on.exit(grDevices::dev.off(), add = TRUE)
  # print(result$plot_combined)
  # grDevices::dev.off()
  # on.exit(NULL)   # clear the on.exit handler after successful dev.off()
  # 
  # -- 10g. Save hex summary sf object ----------------------------------------
  # The full hex summary (geometry + all columns) is saved per species so that
  # region-specific scripts can subset it spatially without re-running this
  # national analysis.
  hex_summary_path <- file.path(
    OUT_HEX_SUMMARIES,
    paste0(spcd_ebird, "_eBird_hex_summary.rds")
  )
  saveRDS(result$hex_summary, hex_summary_path)
  
  # -- 10h. Accumulate statistics ---------------------------------------------
  new_stats <- result$region_stats %>%
    dplyr::mutate(
      spcd_bam   = spcd_ebird,
      spcd_ebird = spcd_ebird,
      sp_english = sp_english,
      model      = "eBird",
      data       = "BAM v5",
      region_id  = "national"
    )
  
  results_list[[spcd_ebird]] <- new_stats
  
  national_stats_df <- dplyr::bind_rows(
    national_stats_df,
    new_stats
  )
  saveRDS(national_stats_df, stats_path)
  completed_species <- c(completed_species, spcd_ebird)
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
