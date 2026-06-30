# ==============================================================================
# Combine BAM and eBird hex-level abundance categories,
# with raw data summaries and BirdLife range information
# ==============================================================================

library(tidyverse)
library(sf)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_CACHE, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

species_lookup_path <- file.path(
  OUT_CACHE,
  "species_lookup_ebird_bam.rds"
)

out_path_rds <- file.path(
  OUT_CACHE,
  "hex_species_abundance_categories_ebird_bam.rds"
)

out_path_csv <- file.path(
  OUT_CACHE,
  "hex_species_abundance_categories_ebird_bam.csv"
)

raw_hex_species_summary_path <- file.path(
  OUT_CACHE,
  "raw_hex_species_summary_bam.rds"
)

# ------------------------------------------------------------------------------
# Species lookup
# ------------------------------------------------------------------------------

species_lookup <- readRDS(species_lookup_path) |>
  # filter(
  #   !is.na(ebird_species_code),
  #   !is.na(bam_species_code)
  # ) |>
  distinct(english_name, scientific_name, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Helper readers
# ------------------------------------------------------------------------------

read_ebird_hex_summary <- function(ebird_species_code) {
  
  path <- file.path(
    OUT_HEX_SUMMARIES_EBIRD,
    paste0(ebird_species_code, "_eBird_hex_summary.rds")
  )
  
  if (!file.exists(path)) return(NULL)
  
  x <- readRDS(path)
  
  x$hex_summary |>
    sf::st_drop_geometry() |>
    transmute(
      hex_id,
      ebird_predicted_abundance = pred_mean,
      ebird_category = as.character(category),
      ebird_coverage = pred_coverage,
      ebird_season = x$season
    )
}

read_bam_hex_summary <- function(bam_species_code) {
  
  path <- file.path(
    OUT_HEX_SUMMARIES_BAM,
    paste0(bam_species_code, "_BAM_hex_summary.rds")
  )
  
  if (!file.exists(path)) return(NULL)
  
  x <- readRDS(path)
  
  x$hex_summary |>
    sf::st_drop_geometry() |>
    transmute(
      hex_id,
      bam_predicted_abundance = pred_mean,
      bam_category = as.character(category),
      bam_coverage = pred_coverage
    )
}

read_birdlife_hex_summary <- function(scientific_name) {
  
  path <- file.path(
    OUT_HEX_SUMMARIES_BIRDLIFE,
    paste0(make.names(scientific_name), "_BirdLife_hex_summary.rds")
  )
  
  if (!file.exists(path)) return(NULL)
  
  x <- readRDS(path)
  
  # Only rows where the hexagon is inside the species' BirdLife range are
  # kept; everything else resolves to NA after the join below, giving the
  # "In range" / NA encoding requested.
  x$hex_summary |>
    sf::st_drop_geometry() |>
    dplyr::filter(in_range_birdlife_any == 1L) |>
    transmute(
      hex_id,
      range_BirdLife = "In range"
    )
}

# ------------------------------------------------------------------------------
# Combine species x hex summaries
# ------------------------------------------------------------------------------

combined_hex_species <- purrr::map_dfr(seq_len(nrow(species_lookup)), function(i) {
  
  sp <- species_lookup[i, ]
  
  ebird_dat    <- read_ebird_hex_summary(sp$ebird_species_code)
  bam_dat      <- read_bam_hex_summary(sp$bam_species_code)
  birdlife_dat <- read_birdlife_hex_summary(sp$scientific_name)
  
  if (is.null(ebird_dat)) {
    ebird_dat <- tibble(
      hex_id = integer(),
      ebird_predicted_abundance = numeric(),
      ebird_category = character(),
      ebird_coverage = numeric(),
      ebird_season = character()
    )
  }
  
  if (is.null(bam_dat)) {
    bam_dat <- tibble(
      hex_id = integer(),
      bam_predicted_abundance = numeric(),
      bam_category = character(),
      bam_coverage = numeric()
    )
  }
  
  if (is.null(birdlife_dat)) {
    birdlife_dat <- tibble(
      hex_id = integer(),
      range_BirdLife = character()
    )
  }
  
  full_join(ebird_dat, bam_dat, by = "hex_id") |>
    left_join(birdlife_dat, by = "hex_id") |>
    mutate(
      english_name = sp$english_name,
      scientific_name = sp$scientific_name,
      ebird_species_code = sp$ebird_species_code,
      bam_species_code = sp$bam_species_code,
      source_status = case_when(
        !is.na(ebird_category) & !is.na(bam_category) ~ "both",
        !is.na(ebird_category) &  is.na(bam_category) ~ "ebird_only",
        is.na(ebird_category) & !is.na(bam_category) ~ "bam_only",
        TRUE ~ "neither"
      )
    ) |>
    select(
      hex_id,
      english_name,
      scientific_name,
      ebird_species_code,
      bam_species_code,
      season = ebird_season,
      bam_predicted_abundance,
      bam_category,
      bam_coverage,
      ebird_predicted_abundance,
      ebird_category,
      ebird_coverage,
      range_BirdLife,
      source_status
    )
})


# ------------------------------------------------------------------------------
# Add raw BAM survey and detection summaries
# ------------------------------------------------------------------------------

if (!file.exists(raw_hex_species_summary_path)) {
  stop(
    "Raw hex-species summary does not exist: ",
    raw_hex_species_summary_path,
    "\nRun script 04 first."
  )
}

raw_hex_species_summary <- readRDS(raw_hex_species_summary_path) |>
  select(
    hex_id,
    ebird_species_code,
    bam_species_code,
    raw_species_code,
    code_source,
    date_label,
    start_doy,
    end_doy,
    n_surveys_bam_raw,
    n_surveys_detected_bam_raw,
    total_count_bam_raw,
    raw_detected_bam
  )

combined_hex_species <- combined_hex_species |>
  left_join(
    raw_hex_species_summary,
    by = c("hex_id", "ebird_species_code", "bam_species_code")
  ) |>
  mutate(
    n_surveys_bam_raw = replace_na(n_surveys_bam_raw, 0L),
    n_surveys_detected_bam_raw = replace_na(n_surveys_detected_bam_raw, 0L),
    total_count_bam_raw = replace_na(total_count_bam_raw, 0),
    raw_detected_bam = replace_na(raw_detected_bam, FALSE)
  )

# ------------------------------------------------------------------------------
# Limit to columns of interest
# ------------------------------------------------------------------------------

dat_to_save <- combined_hex_species %>%
  dplyr::select(hex_id,english_name,scientific_name,season,ebird_category,bam_category,range_BirdLife,date_label,n_surveys_bam_raw,n_surveys_detected_bam_raw) %>%
  dplyr::rename(date_range = date_label,
                n_surveys = n_surveys_bam_raw,
                n_detections = n_surveys_detected_bam_raw)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

saveRDS(dat_to_save, out_path_rds)
readr::write_csv(dat_to_save, out_path_csv)

message("Rows: ", nrow(dat_to_save))
message("Species: ", n_distinct(dat_to_save$english_name))
message("Hexagons: ", n_distinct(dat_to_save$hex_id, na.rm = TRUE))
message("Saved RDS: ", out_path_rds)
message("Saved CSV: ", out_path_csv)

# ------------------------------------------------------------------------------
# Summary features
# ------------------------------------------------------------------------------

table(dat_to_save$hex_id)       # 455 species in every hexagon
table(dat_to_save$english_name) # 83276 hexagons per species

length(unique(species_lookup$bam_species_code))
species_lookup


# Number of species with non-zero detections in at least one hexagon
species_detected <- dat_to_save %>%
  subset(n_detections>0)
length(unique(species_detected$english_name)) # 417 species detected at least once

# BirdLife range coverage check
table(dat_to_save$range_BirdLife, useNA = "always")
mean(!is.na(dat_to_save$range_BirdLife)) # proportion of hex-species rows flagged "In range"