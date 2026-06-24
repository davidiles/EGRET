# ==============================================================================
# 04_summarize_raw_data_by_hexagon.R
#
# Summarize raw BAM survey data by national 10-km hexagon and species,
# using species-specific eBird breeding-season date windows
# ==============================================================================

library(tidyverse)
library(sf)
library(lubridate)
library(ebirdst)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_CACHE, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Output paths
# ------------------------------------------------------------------------------

bam_hex_path <- file.path(OUT_CACHE, "bam_sf_with_hex_id.rds")

raw_hex_species_summary_path <- file.path(
  OUT_CACHE,
  "raw_hex_species_summary_bam.rds"
)

raw_hex_species_summary_csv_path <- file.path(
  OUT_CACHE,
  "raw_hex_species_summary_bam.csv"
)

raw_species_unmatched_path <- file.path(
  OUT_CACHE,
  "raw_species_codes_unmatched_bam.csv"
)

# ------------------------------------------------------------------------------
# National study area
# ------------------------------------------------------------------------------

region_national_path <- file.path(OUT_CACHE, "region_national.rds")

if (!file.exists(region_national_path)) {
  stop("Missing region_national file: ", region_national_path)
}

region_national <- readRDS(region_national_path)

# ------------------------------------------------------------------------------
# Hex grid
# ------------------------------------------------------------------------------

hex_grid_path <- file.path(
  OUT_CACHE,
  paste0("hex_grid_national_", HEX_WIDTH_KM, "_km.rds")
)

if (!file.exists(hex_grid_path)) {
  stop(
    "Hex grid does not exist: ", hex_grid_path, "\n",
    "Run scripts/01_create_national_hex_grid.R first."
  )
}

hex_grid_national <- readRDS(hex_grid_path)

stopifnot(inherits(hex_grid_national, "sf"))
stopifnot("hex_id" %in% names(hex_grid_national))
stopifnot("hex_area_km2" %in% names(hex_grid_national))

message("Loaded hex grid with ", nrow(hex_grid_national), " hexagons.")

# ------------------------------------------------------------------------------
# Species lookup
# ------------------------------------------------------------------------------

species_lookup_path <- file.path(OUT_CACHE, "species_lookup_ebird_bam.rds")

if (!file.exists(species_lookup_path)) {
  stop(
    "Species lookup table does not exist: ", species_lookup_path, "\n",
    "Run the species lookup creation script first."
  )
}

species_table <- readRDS(species_lookup_path)

stopifnot(
  "english_name" %in% names(species_table),
  "scientific_name" %in% names(species_table),
  "bbl_code" %in% names(species_table),
  "bam_species_code" %in% names(species_table),
  "ebird_species_code" %in% names(species_table)
)

if (!"AOS_CODE_4L" %in% names(species_table)) {
  species_table$AOS_CODE_4L <- NA_character_
}

# ------------------------------------------------------------------------------
# Load and filter BAM raw survey data
# ------------------------------------------------------------------------------

load(BAM_DATA_PATH) # loads object 'dat'

if (!exists("dat")) {
  stop("BAM_DATA_PATH did not load an object named 'dat'.")
}

bam_sf <- dat |>
  filter(!is.na(longitude), !is.na(latitude)) |>
  mutate(
    survey_datetime = lubridate::ymd_hms(date_time, tz = "UTC"),
    day_of_year     = lubridate::yday(survey_datetime),
    year            = lubridate::year(survey_datetime)
  ) |>
  filter(
    duration >= SURVEY_DURATION_MIN_S,
    duration <= SURVEY_DURATION_MAX_S,
    distance >  SURVEY_DISTANCE_MIN_M
  ) |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs    = 4326,
    remove = FALSE
  ) |>
  st_transform(st_crs(region_national)) |>
  st_filter(region_national, .predicate = st_intersects)

rm(dat)
gc()

message("Filtered BAM surveys in study area: ", nrow(bam_sf))

# ------------------------------------------------------------------------------
# Assign each BAM survey to a hexagon
# ------------------------------------------------------------------------------

bam_sf <- load_or_compute(
  cache_path = bam_hex_path,
  compute_fn = function() {
    precompute_survey_hex_id(
      survey_sf = bam_sf,
      hex_grid  = hex_grid_national
    )
  },
  source_path = BAM_DATA_PATH
)

bam_sf <- bam_sf |>
  filter(!is.na(hex_id))

message("BAM surveys assigned to hexagons: ", nrow(bam_sf))

# ------------------------------------------------------------------------------
# Identify raw species columns
# ------------------------------------------------------------------------------

if (!all(c("ABBE", "YTWA") %in% names(bam_sf))) {
  stop("Could not find expected species column range ABBE:YTWA.")
}

species_cols <- names(bam_sf)[
  which(names(bam_sf) == "ABBE"):which(names(bam_sf) == "YTWA")
]

species_cols_detected <- bam_sf |>
  st_drop_geometry() |>
  summarise(
    across(
      all_of(species_cols),
      ~ sum(.x > 0, na.rm = TRUE)
    )
  ) |>
  pivot_longer(
    cols = everything(),
    names_to = "raw_species_code",
    values_to = "n_surveys_detected_total"
  ) |>
  filter(n_surveys_detected_total > 0) |>
  pull(raw_species_code)

message("Species with >=1 raw detection: ", length(species_cols_detected))

# ------------------------------------------------------------------------------
# Build raw-code-to-species lookup
# ------------------------------------------------------------------------------

species_code_lookup <- species_table |>
  select(
    english_name,
    scientific_name,
    bbl_code,
    AOS_CODE_4L,
    bam_species_code,
    ebird_species_code
  ) |>
  mutate(
    bam_species_code_final = bam_species_code
  ) |>
  pivot_longer(
    cols = c(bbl_code, AOS_CODE_4L, bam_species_code),
    names_to = "code_source",
    values_to = "raw_species_code"
  ) |>
  filter(!is.na(raw_species_code), raw_species_code != "") |>
  arrange(
    raw_species_code,
    desc(code_source == "bam_species_code"),
    desc(code_source == "AOS_CODE_4L"),
    desc(code_source == "bbl_code")
  ) |>
  distinct(raw_species_code, .keep_all = TRUE) |>
  rename(
    bam_species_code = bam_species_code_final
  )

species_to_summarize <- species_code_lookup |>
  filter(
    raw_species_code %in% species_cols,
    !is.na(english_name),
    !is.na(ebird_species_code)
  ) |>
  mutate(
    code_priority = case_when(
      code_source == "bam_species_code" ~ 1L,
      code_source == "AOS_CODE_4L" ~ 2L,
      code_source == "bbl_code" ~ 3L,
      TRUE ~ 99L
    )
  ) |>
  arrange(
    english_name,
    scientific_name,
    ebird_species_code,
    bam_species_code,
    code_priority
  ) |>
  distinct(
    english_name,
    scientific_name,
    ebird_species_code,
    bam_species_code,
    .keep_all = TRUE
  ) |>
  select(-code_priority)

species_to_summarize |>
  count(english_name, scientific_name, ebird_species_code, bam_species_code) |>
  filter(n > 1)


message("Species to summarize with date filtering: ", nrow(species_to_summarize))

raw_species_unmatched <- tibble(
  raw_species_code = setdiff(species_cols_detected, species_to_summarize$raw_species_code)
) |>
  arrange(raw_species_code)

message("Unmatched raw species codes: ", nrow(raw_species_unmatched))

# ------------------------------------------------------------------------------
# Species-specific summary function
# ------------------------------------------------------------------------------

summarize_species_raw_by_hex <- function(sp_row, bam_sf) {
  
  sp_english <- sp_row$english_name
  raw_code   <- sp_row$raw_species_code
  
  date_window <- get_species_date_window(
    sp_english     = sp_english,
    ebirdst_runs   = ebirdst::ebirdst_runs,
    buffer_days    = DATE_BUFFER_DAYS,
    fallback_start = FALLBACK_START_DATE,
    fallback_end   = FALLBACK_END_DATE
  )
  
  start_doy <- date_window$start_doy
  end_doy   <- date_window$end_doy
  
  if (start_doy <= end_doy) {
    sp_surveys <- bam_sf |>
      filter(day_of_year >= start_doy, day_of_year <= end_doy)
  } else {
    sp_surveys <- bam_sf |>
      filter(day_of_year >= start_doy | day_of_year <= end_doy)
  }
  
  if (nrow(sp_surveys) == 0) {
    return(tibble())
  }
  
  sp_surveys |>
    st_drop_geometry() |>
    group_by(hex_id) |>
    summarise(
      n_surveys_bam_raw = n(),
      n_surveys_detected_bam_raw = sum(.data[[raw_code]] > 0, na.rm = TRUE),
      total_count_bam_raw = sum(.data[[raw_code]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      english_name = sp_row$english_name,
      scientific_name = sp_row$scientific_name,
      ebird_species_code = sp_row$ebird_species_code,
      bam_species_code = sp_row$bam_species_code,
      raw_species_code = raw_code,
      code_source = sp_row$code_source,
      date_label = date_window$date_label,
      start_doy = start_doy,
      end_doy = end_doy,
      raw_detected_bam = n_surveys_detected_bam_raw > 0
    ) |>
    select(
      hex_id,
      english_name,
      scientific_name,
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
}

# ------------------------------------------------------------------------------
# Run species-specific summaries
# ------------------------------------------------------------------------------

raw_hex_species_summary <- purrr::map_dfr(
  seq_len(nrow(species_to_summarize)),
  function(i) {
    summarize_species_raw_by_hex(
      sp_row = species_to_summarize[i, ],
      bam_sf = bam_sf
    )
  }
) |>
  arrange(english_name, hex_id)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

saveRDS(raw_hex_species_summary, raw_hex_species_summary_path)

readr::write_csv(
  raw_hex_species_summary,
  raw_hex_species_summary_csv_path
)

readr::write_csv(
  raw_species_unmatched,
  raw_species_unmatched_path
)

message("Saved raw hex-species summary: ", raw_hex_species_summary_path)
message("Rows: ", nrow(raw_hex_species_summary))
message("Unmatched raw species codes: ", nrow(raw_species_unmatched))