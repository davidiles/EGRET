# ==============================================================================
# Extract BirdLife range-boundary overlap by national hexagon
#
# Purpose:
#   For each species in the national species lookup table, identify whether each
#   national hexagon intersects the species' BirdLife resident or breeding range.
#
# Output:
#   One RDS file per species in OUT_HEX_SUMMARIES_BIRDLIFE, analogous in
#   structure to the eBird/BAM species-specific outputs, but with 0/1 range
#   indicators rather than relative-abundance predictions or categories.
#   BirdLife seasonal codes used:
#     1 = resident
#     2 = breeding
# ==============================================================================

library(tidyverse)
library(sf)
library(readxl)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_HEX_SUMMARIES_BIRDLIFE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_CACHE,                  recursive = TRUE, showWarnings = FALSE)

sf::sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# Species list
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
  "english_name"    %in% names(species_table),
  "scientific_name" %in% names(species_table)
)

optional_name_cols <- c("AOS_common_name", "ebird_common_name", "bam_common_name")
for (nm in optional_name_cols) {
  if (!nm %in% names(species_table)) species_table[[nm]] <- NA_character_
}

species_table <- species_table |>
  dplyr::distinct(scientific_name, .keep_all = TRUE) |>
  dplyr::filter(!is.na(scientific_name))

species_vec <- species_table$scientific_name
message("Species to process: ", length(species_vec))

# ------------------------------------------------------------------------------
# National study area and hex grid
# ------------------------------------------------------------------------------

region_national_path <- file.path(OUT_CACHE, "region_national.rds")
if (!file.exists(region_national_path)) {
  stop("National region file does not exist: ", region_national_path)
}
region_national <- readRDS(region_national_path)
stopifnot(inherits(region_national, "sf"))

hex_grid_path <- file.path(
  OUT_CACHE,
  paste0("hex_grid_national_", HEX_WIDTH_KM, "_km.rds")
)
if (!file.exists(hex_grid_path)) {
  stop("Hex grid does not exist: ", hex_grid_path)
}

hex_grid_national <- readRDS(hex_grid_path) |>
  sf::st_transform(sf::st_crs(region_national))

stopifnot(
  inherits(hex_grid_national, "sf"),
  "hex_id"       %in% names(hex_grid_national),
  "hex_area_km2" %in% names(hex_grid_national)
)

message("Loaded hex grid with ", nrow(hex_grid_national), " hexagons.")

hex_template <- hex_grid_national |>
  dplyr::select(hex_id, hex_area_km2, geometry)

# ------------------------------------------------------------------------------
# Load and prepare BirdLife layer
# ------------------------------------------------------------------------------

message("Reading BirdLife range polygons...")

birdlife <- sf::st_read(BIRDLIFE_RANGE_PATH, quiet = TRUE) |>
  dplyr::filter(seasonal %in% c(1L, 2L))

# Normalise the species ID column name.
if (!"sisid" %in% names(birdlife)) {
  candidates <- names(birdlife)[tolower(names(birdlife)) %in% c("sisid", "sisrecid")]
  if (length(candidates) == 0) {
    stop("Cannot find 'sisid' or 'SISRecID' column in BirdLife range file.")
  }
  names(birdlife)[names(birdlife) == candidates[1]] <- "sisid"
}
if (!"seasonal" %in% names(birdlife)) {
  stop("Column 'seasonal' not found in BirdLife range file.")
}

birdlife <- birdlife |>
  dplyr::rename(geometry = geom) |>
  dplyr::select(sisid, seasonal, geometry) |>
  dplyr::mutate(sisid = as.character(sisid))

message("BirdLife polygons loaded (resident + breeding, global): ", nrow(birdlife))

# ------------------------------------------------------------------------------
# Bbox crop and validity repair — done ONCE on the whole layer
# ------------------------------------------------------------------------------

message("Cropping BirdLife layer to Canadian bounding box...")

canada_bbox <- sf::st_bbox(
  sf::st_transform(hex_template, sf::st_crs(birdlife))
) |>
  sf::st_as_sfc()
sf::st_crs(canada_bbox) <- sf::st_crs(birdlife)

birdlife <- suppressWarnings(sf::st_crop(birdlife, canada_bbox))
message("After bbox crop: ", nrow(birdlife), " polygons remain.")

birdlife <- suppressWarnings(sf::st_make_valid(birdlife))
birdlife <- sf::st_transform(birdlife, sf::st_crs(hex_template))

# Attach a row index that will survive the join below.
birdlife$bl_row <- seq_len(nrow(birdlife))

# ------------------------------------------------------------------------------
# SINGLE GLOBAL st_intersects() CALL
#
# This is the only spatial predicate in the entire script.
# GEOS builds one STR-tree over all ~n_polygons BirdLife features; all
# 120,000 hex queries run against that index in one call.
# The result is a sparse Matrix (hex × polygon). We convert it immediately
# to a compact long data frame and drop the geometry objects.
# ------------------------------------------------------------------------------

message(
  "Running st_intersects() across all ", nrow(birdlife),
  " BirdLife polygons and ", nrow(hex_template), " hexagons..."
)

hex_centroids <- sf::st_centroid(hex_template)
hits_sparse   <- sf::st_intersects(hex_centroids, birdlife)

message("Spatial index query complete. Building lookup table...")

# Convert sparse list to a two-column data frame.
# hits_sparse[[i]] is an integer vector of birdlife row indices that hex i hits.
hex_ids     <- hex_template$hex_id
hits_long <- tibble::tibble(
  hex_id = rep(hex_ids, lengths(hits_sparse)),
  bl_row = unlist(hits_sparse, use.names = FALSE)
)

# Join sisid and seasonal from the (now geometry-free) birdlife attribute table.
birdlife_attr <- sf::st_drop_geometry(birdlife)  # drop heavy geometry column

hits_long <- hits_long |>
  dplyr::left_join(
    birdlife_attr[, c("bl_row", "sisid", "seasonal")],
    by = "bl_row"
  )

# Pre-split into resident and breeding subsets for fast per-species lookup.
hits_any      <- hits_long                           # seasons 1 and 2 combined
hits_resident <- hits_long[hits_long$seasonal == 1L, ]
hits_breeding <- hits_long[hits_long$seasonal == 2L, ]

message(
  "Lookup table built: ", nrow(hits_long), " hex-polygon intersections recorded."
)

# Free the geometry-heavy BirdLife sf object; it is no longer needed.
rm(birdlife, hits_sparse)

# ------------------------------------------------------------------------------
# BirdLife species attribute table
# ------------------------------------------------------------------------------

birdlife_species <- readxl::read_xlsx(BIRDLIFE_SPECIES_PATH) |>
  dplyr::select(`Common name`, `Scientific name`, SISRecID) |>
  dplyr::distinct() |>
  dplyr::mutate(
    SISRecID          = as.character(SISRecID),
    `Common name`     = as.character(`Common name`),
    `Scientific name` = as.character(`Scientific name`)
  )

# ------------------------------------------------------------------------------
# Build BirdLife matching table
# ------------------------------------------------------------------------------

birdlife_lookup <- species_table |>
  dplyr::mutate(
    birdlife_match_names = purrr::pmap_chr(
      list(english_name, AOS_common_name, ebird_common_name, bam_common_name),
      function(english_name, AOS_common_name, ebird_common_name, bam_common_name) {
        paste(unique(na.omit(c(
          english_name, AOS_common_name, ebird_common_name, bam_common_name
        ))), collapse = " | ")
      }
    )
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    birdlife_sisids = list({
      common_names <- unique(na.omit(c(
        english_name, AOS_common_name, ebird_common_name, bam_common_name
      )))
      matched <- birdlife_species |>
        dplyr::filter(
          `Scientific name` %in% scientific_name |
            `Common name`     %in% common_names
        )
      unique(as.character(matched$SISRecID))
    }),
    n_birdlife_matches = length(birdlife_sisids),
    in_birdlife        = as.integer(n_birdlife_matches > 0)
  ) |>
  dplyr::ungroup()

message(
  "Species in BirdLife: ", sum(birdlife_lookup$in_birdlife),
  " / ", nrow(birdlife_lookup)
)

readr::write_csv(
  birdlife_lookup |>
    dplyr::select(
      english_name, scientific_name, birdlife_match_names,
      birdlife_sisids, n_birdlife_matches, in_birdlife
    ),
  file.path(OUT_HEX_SUMMARIES_BIRDLIFE, "birdlife_species_matching_status.csv")
)

# All hex_ids as a reference vector for building the output template.
all_hex_ids <- hex_template$hex_id

# ------------------------------------------------------------------------------
# Per-species loop — NO spatial operations
#
# For each species:
#   1. Look up its BirdLife SISRecID(s).
#   2. Filter the pre-built hits_long table to those IDs.
#   3. Assign 0/1 flags by simple set membership.
#   4. Join onto hex_template and save.
# ------------------------------------------------------------------------------

status_path <- file.path(
  OUT_HEX_SUMMARIES_BIRDLIFE,
  "birdlife_hex_extraction_status.csv"
)

status_list <- vector("list", length(species_vec))

for (i in seq_along(species_vec)) {
  
  sp_scientific <- species_vec[i]
  
  sp_row <- dplyr::filter(birdlife_lookup, scientific_name == sp_scientific) |>
    dplyr::slice(1)
  
  sp_english    <- sp_row$english_name
  sp_file_stub  <- make.names(sp_scientific)
  hex_summary_path <- file.path(
    OUT_HEX_SUMMARIES_BIRDLIFE,
    paste0(sp_file_stub, "_BirdLife_hex_summary.rds")
  )
  
  message(
    "[", i, " / ", length(species_vec), "] ",
    sp_english, " (", sp_scientific, ")"
  )
  
  if (file.exists(hex_summary_path)) {
    status_list[[i]] <- tibble::tibble(
      scientific_name    = sp_scientific,
      sp_english         = sp_english,
      status             = "skipped_existing",
      n_birdlife_matches = sp_row$n_birdlife_matches,
      n_hex_any          = NA_integer_,
      n_hex_resident     = NA_integer_,
      n_hex_breeding     = NA_integer_,
      elapsed_min        = 0,
      error              = NA_character_
    )
    next
  }
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    
    sp_sisids <- unlist(sp_row$birdlife_sisids)
    
    # Hexes hit by ANY seasonal polygon for this species.
    hex_any      <- unique(hits_any     [hits_any$sisid      %in% sp_sisids, "hex_id", drop = TRUE])
    hex_resident <- unique(hits_resident[hits_resident$sisid %in% sp_sisids, "hex_id", drop = TRUE])
    hex_breeding <- unique(hits_breeding[hits_breeding$sisid %in% sp_sisids, "hex_id", drop = TRUE])
    
    hex_summary <- hex_template |>
      dplyr::mutate(
        in_range_birdlife_any      = as.integer(hex_id %in% hex_any),
        in_range_birdlife_resident = as.integer(hex_id %in% hex_resident),
        in_range_birdlife_breeding = as.integer(hex_id %in% hex_breeding)
      )
    
    species_output <- list(
      sp_english         = sp_english,
      scientific_name    = sp_scientific,
      birdlife_sisids    = sp_sisids,
      n_birdlife_matches = sp_row$n_birdlife_matches,
      range_source       = BIRDLIFE_RANGE_PATH,
      species_source     = BIRDLIFE_SPECIES_PATH,
      seasons_used       = c("resident" = 1L, "breeding" = 2L),
      hex_summary        = hex_summary
    )
    
    saveRDS(species_output, hex_summary_path)
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    sp_status <- if (length(sp_sisids) == 0) "no_birdlife_match" else "success"
    
    tibble::tibble(
      scientific_name    = sp_scientific,
      sp_english         = sp_english,
      status             = sp_status,
      n_birdlife_matches = sp_row$n_birdlife_matches,
      n_hex_any          = sum(hex_summary$in_range_birdlife_any),
      n_hex_resident     = sum(hex_summary$in_range_birdlife_resident),
      n_hex_breeding     = sum(hex_summary$in_range_birdlife_breeding),
      elapsed_min        = elapsed,
      error              = NA_character_
    )
    
  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    tibble::tibble(
      scientific_name    = sp_scientific,
      sp_english         = sp_english,
      status             = "failed",
      n_birdlife_matches = sp_row$n_birdlife_matches,
      n_hex_any          = NA_integer_,
      n_hex_resident     = NA_integer_,
      n_hex_breeding     = NA_integer_,
      elapsed_min        = elapsed,
      error              = conditionMessage(e)
    )
  })
  
  status_list[[i]] <- result
  
  # Write progress after every species so the run is resumable.
  readr::write_csv(dplyr::bind_rows(status_list), status_path)
}

extraction_status <- dplyr::bind_rows(status_list)
readr::write_csv(extraction_status, status_path)

message("Finished BirdLife range-boundary extraction.")
message("Successful:   ", sum(extraction_status$status == "success"))
message("Skipped:      ", sum(extraction_status$status == "skipped_existing"))
message("No match:     ", sum(extraction_status$status == "no_birdlife_match"))
message("Failed:       ", sum(extraction_status$status == "failed"))
message("Status file:  ", status_path)

# ------------------------------------------------------------------------------
# Example plot for a single species
# ------------------------------------------------------------------------------

successful_species <- extraction_status |>
  dplyr::filter(status %in% c("success", "skipped_existing")) |>
  dplyr::slice(10)

if (nrow(successful_species) > 0) {
  
  example_file <- file.path(
    OUT_HEX_SUMMARIES_BIRDLIFE,
    paste0(make.names(successful_species$scientific_name[1]), "_BirdLife_hex_summary.rds")
  )
  
  if (file.exists(example_file)) {
    example_result <- readRDS(example_file)
    print(
      ggplot(example_result$hex_summary) +
        geom_sf(aes(fill = factor(in_range_birdlife_any)), col = NA) +
        scale_fill_manual(
          values = c("0" = "white", "1" = "dodgerblue"),
          name   = "BirdLife range",
          labels = c("0" = "Outside range", "1" = "Intersects range")
        ) +
        geom_sf(data = region_national, col = "black", fill = "transparent") +
        ggtitle(paste0(
          example_result$sp_english,
          " — BirdLife resident/breeding range"
        )) +
        theme_bw()
    )
  }
}
