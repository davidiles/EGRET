# ==============================================================================
# config.R
#
# Central configuration for the BAM/eBird model assessment pipeline.
# All file paths, analysis constants, and output directories are defined here.
# Both assess_BAM_national.R and assess_ebird_national.R source this file.
#
# To run the pipeline on a different machine, only this file needs to be edited.
# ==============================================================================


# ------------------------------------------------------------------------------
# Input file paths
#
# Override any of these with environment variables for portability, e.g.:
#   Sys.setenv(BAM_RASTER_PATH = "/data/BAMv5/2020")
# If the environment variable is not set, the default Windows path is used.
# ------------------------------------------------------------------------------

BAM_RASTER_PATH <- Sys.getenv(
  "BAM_RASTER_PATH",
  unset = "C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/BAMv5_model_outputs/2020"
)

BAM_DATA_PATH <- Sys.getenv(
  "BAM_DATA_PATH",
  unset = "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/BAMDataset/04_BAMDataset_WT-2026-03-09_EBd-Jan-2026.Rdata"
)

BCR_GPKG_PATH <- Sys.getenv(
  "BCR_GPKG_PATH",
  unset = "C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/BCR_2026/bcr_2026g.gpkg"
)

AVIAN_CORE_PATH <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/Avian_Core/Avian_Core_20251124.xlsx"
AOS_PATH <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/AOS_BIRD_CODES/IBP-AOS-list25.csv"
BIRDLIFE_RANGE_PATH <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/BirdLife_Range_Boundaries/BOTW_2025.gpkg"
BIRDLIFE_SPECIES_PATH <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/BirdLife_Range_Boundaries/birdlife_species.xlsx"

# ------------------------------------------------------------------------------
# Output directories
# ------------------------------------------------------------------------------

OUT_ROOT <- "output"

# National-scale outputs (shared by both assess scripts)
# OUT_NATIONAL_BAM_FIGURES  <- file.path(OUT_ROOT, "National", "BAM",   "figures")
# OUT_NATIONAL_BAM_STATS    <- file.path(OUT_ROOT, "National", "BAM",   "stats")
# OUT_NATIONAL_EBIRD_FIGURES <- file.path(OUT_ROOT, "National", "eBird", "figures")
OUT_NATIONAL_EBIRD_STATS   <- file.path(OUT_ROOT, "National", "eBird", "stats")

# Per-species hex summaries — written by both assess scripts, read by any
# region-specific script that wants to subset without rerunning the full analysis.
# Files are named: <OUT_HEX_SUMMARIES>/<spcd>_BAM_hex_summary.rds
#                  <OUT_HEX_SUMMARIES>/<spcd>_eBird_hex_summary.rds
OUT_HEX_SUMMARIES_BAM <- file.path(OUT_ROOT, "hex_summaries/BAM")
OUT_HEX_SUMMARIES_EBIRD <- file.path(OUT_ROOT, "hex_summaries/eBird")
OUT_HEX_SUMMARIES_BIRDLIFE <- file.path(OUT_ROOT, "hex_summaries/BirdLife")

# Shared spatial cache files (hex grid and precomputed metadata)
OUT_CACHE <- file.path(OUT_ROOT, "cache")


# ------------------------------------------------------------------------------
# Survey filter parameters
#
# Applied to the BAM dataset before any species loop.
# Both the BAM and eBird scripts use the same survey filters so that
# predictions from the two models are assessed against an identical observation
# dataset.
# ------------------------------------------------------------------------------

# Minimum and maximum survey duration (seconds)
SURVEY_DURATION_MIN_S <- 1  * 60   #  1 minute
SURVEY_DURATION_MAX_S <- 10 * 60   # 10 minutes

# Minimum detection radius (metres)
SURVEY_DISTANCE_MIN_M <- 50

# Minimum number of BAM detections (across all surveys) for a species to be
# included in either analysis.  Set to 1 so that any species ever detected
# at least once is assessed.
MIN_DETECTIONS <- 1L


# ------------------------------------------------------------------------------
# Date window parameters
#
# Both scripts derive a species-specific breeding date window from eBird
# seasonal metadata (breeding_start / breeding_end).  The window is expanded
# by DATE_BUFFER_DAYS on each side.  If a species has no eBird metadata, the
# fallback dates below are used.
# ------------------------------------------------------------------------------

DATE_BUFFER_DAYS <- 7L

FALLBACK_START_DATE <- "2023-06-01"
FALLBACK_END_DATE   <- "2023-07-15"


# ------------------------------------------------------------------------------
# Hex grid parameters
# ------------------------------------------------------------------------------

HEX_WIDTH_KM <- 10


# ------------------------------------------------------------------------------
# Raster extraction and transformation constants
# ------------------------------------------------------------------------------

# Upper quantile used to cap extreme raster values before visualisation and
# extraction.  Values above this quantile are clamped to the quantile value.
RAST_MAX_Q <- 0.995

# Variance-stabilising transformation applied to predicted values in plots.
# Options: "identity", "sqrt", "log" (anything accepted by ggplot2 trans=).
TRANSFORM <- "identity"

# ------------------------------------------------------------------------------
# Cumulative-population threshold parameters
#
# get_thresholds_cumpop() classifies every hexagon into one of five abundance
# categories by walking the predicted-abundance distribution from lowest to
# highest and finding the pred_mean value at which cumulative population
# crosses each breakpoint.
#
# All five breakpoints are defined here so they can be tuned in one place.
# Pass them to get_thresholds_cumpop() via the cumpop_breaks argument (see
# functions.R).  The names must match the expected list structure:
#   absent_upper  – upper edge of "Effectively Absent"
#   low_upper     – upper edge of "Low"
#   mod_upper     – upper edge of "Moderate"
#   high_upper    – upper edge of "High"
#   (anything above high_upper is "Very High")
#
# The absence threshold also drives the flag logic:
#   flag_pred_absent_detected    uses absent_upper
#   flag_high_pred_no_detection  uses high_upper
# ------------------------------------------------------------------------------

CUMPOP_BREAKS <- list(
  absent_upper = 0.001,   # hexagons summing to the lowest 0.1% of population
  low_upper    = 0.25,    # 0.1% – 25%  of cumulative population
  mod_upper    = 0.50,    # 25%  – 50%
  high_upper   = 0.75     # 50%  – 75%  (above = "Very High")
)

FLAG_HIGH_PRED_MIN_SURVEYS    <- 20   # More than 20 surveys conducted but species was not detected
FLAG_ABSENCE_MIN_SURVEYS      <- 2    # Species must have been observed twice in a hexagon categorized as absent
FLAG_LINEWIDTH                <- 0.2
