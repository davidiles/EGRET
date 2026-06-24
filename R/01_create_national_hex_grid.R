# ==============================================================================
# Create national 10-km hexagon grid
# ==============================================================================

library(tidyverse)
library(sf)
library(lwgeom)
library(rnaturalearth)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_CACHE, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Output paths
# ------------------------------------------------------------------------------

hex_grid_path <- file.path(
  OUT_CACHE,
  paste0("hex_grid_national_", HEX_WIDTH_KM, "_km.rds")
)

hex_grid_gpkg_path <- file.path(
  OUT_CACHE,
  paste0("hex_grid_national_", HEX_WIDTH_KM, "_km.gpkg")
)

region_national_path <- file.path(
  OUT_CACHE,
  "region_national.rds"
)

region_national_gpkg_path <- file.path(
  OUT_CACHE,
  "region_national.gpkg"
)

# ------------------------------------------------------------------------------
# Study area
# ------------------------------------------------------------------------------

study_area <- rnaturalearth::ne_states(country = "Canada") |>
  sf::st_transform(3978) |>
  sf::st_make_valid()

bcr <- sf::st_read(BCR_GPKG_PATH, quiet = TRUE) |>
  sf::st_transform(sf::st_crs(study_area)) |>
  sf::st_make_valid() |>
  sf::st_intersection(study_area) |>
  dplyr::filter(!(bcr_label %in% BCR_EXCLUDE))

region_national <- sf::st_union(bcr) |>
  sf::st_as_sf() |>
  remove_multipolygon_holes()

# ------------------------------------------------------------------------------
# Create hex grid
# ------------------------------------------------------------------------------

hex_grid_national <- make_hex_grid(
  study_area    = region_national,
  hex_width_km  = HEX_WIDTH_KM
)

# Add metadata used later for thresholding / summaries
hex_grid_national <- precompute_hex_metadata(
  hex_grid   = hex_grid_national,
  study_area = region_national
)

# Basic checks
stopifnot(inherits(hex_grid_national, "sf"))
stopifnot("hex_id" %in% names(hex_grid_national))
stopifnot("hex_area_km2" %in% names(hex_grid_national))

message("Created hex grid with ", nrow(hex_grid_national), " hexagons.")

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

saveRDS(hex_grid_national, hex_grid_path)

sf::st_write(
  hex_grid_national,
  hex_grid_gpkg_path,
  layer = "hex_grid_national",
  delete_dsn = TRUE,
  quiet = TRUE
)

saveRDS(region_national, region_national_path)

sf::st_write(
  region_national,
  region_national_gpkg_path,
  layer = "region_national",
  delete_dsn = TRUE,
  quiet = TRUE
)

message("Saved hex grid RDS:      ", hex_grid_path)
message("Saved hex grid GPKG:     ", hex_grid_gpkg_path)
message("Saved region RDS:        ", region_national_path)
message("Saved region GPKG:       ", region_national_gpkg_path)
