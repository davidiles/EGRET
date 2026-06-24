# ==============================================================================
# Extract predictions from BAM abundance rasters
# ==============================================================================

library(tidyverse)
library(terra)
library(sf)
library(exactextractr)
library(future)
library(furrr)
library(progressr)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_HEX_SUMMARIES_BAM, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_CACHE, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Species list
# ------------------------------------------------------------------------------

species_lookup_path <- file.path(
  OUT_CACHE,
  "species_lookup_ebird_bam.rds"
)

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
  "bam_species_code" %in% names(species_table),
  "bam_common_name" %in% names(species_table)
)

species_table <- species_table |>
  filter(!is.na(bam_species_code)) |>
  distinct(bam_species_code, .keep_all = TRUE)

species_vec <- species_table$bam_species_code

message("Species with BAM products: ", length(species_vec))

# ------------------------------------------------------------------------------
# National study area
# ------------------------------------------------------------------------------

region_national <- readRDS(
  file.path(OUT_CACHE, "region_national.rds")
)

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
# Parallel setup
# ------------------------------------------------------------------------------

n_workers <- 16

future::plan(
  future::multisession,
  workers = n_workers
)

message("Using ", n_workers, " parallel workers.")

progressr::handlers(global = TRUE)

# ------------------------------------------------------------------------------
# Species processing function
# ------------------------------------------------------------------------------

process_bam_species <- function(spcd_bam,
                                species_table,
                                hex_grid_national,
                                bam_raster_dir,
                                out_dir,
                                cumpop_breaks) {
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(terra)
    library(sf)
    library(exactextractr)
  })
  
  sp_row <- species_table |>
    filter(bam_species_code == spcd_bam) |>
    slice(1)
  
  sp_english <- sp_row$english_name
  sp_scientific <- sp_row$scientific_name
  bam_common_name <- sp_row$bam_common_name
  
  hex_summary_path <- file.path(
    out_dir,
    paste0(spcd_bam, "_BAM_hex_summary.rds")
  )
  
  if (file.exists(hex_summary_path)) {
    return(tibble(
      spcd_bam = spcd_bam,
      sp_english = sp_english,
      status = "skipped_existing",
      elapsed_min = 0,
      error = NA_character_
    ))
  }
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    
    terra_tmp <- file.path(tempdir(), paste0("terra_", Sys.getpid()))
    dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
    terra::terraOptions(tempdir = terra_tmp)
    
    bam_raster_path_sp <- file.path(
      bam_raster_dir,
      paste0(spcd_bam, "_Canada_2020.tif")
    )
    
    if (!file.exists(bam_raster_path_sp)) {
      stop("BAM raster not found: ", bam_raster_path_sp)
    }
    
    bam_rast_all <- terra::rast(bam_raster_path_sp)
    
    if (!"mean" %in% names(bam_rast_all)) {
      stop(
        "Layer 'mean' not found in BAM raster. Available layers: ",
        paste(names(bam_rast_all), collapse = ", ")
      )
    }
    
    bam_rast <- bam_rast_all[["mean"]]
    
    hex_summary <- extract_hex_predictions(
      hex_grid     = hex_grid_national,
      rast         = bam_rast,
      min_coverage = 0.25
    )
    
    thresholds <- get_thresholds_cumpop(
      hex_summary   = hex_summary,
      col           = "pred_mean",
      coverage_col  = "pred_coverage",
      area_col      = "hex_area_km2",
      cumpop_breaks = cumpop_breaks
    )
    
    hex_summary <- apply_thresholds(
      hex_summary = hex_summary,
      thresholds  = thresholds,
      col         = "pred_mean",
      out_col     = "category"
    )
    
    species_output <- list(
      sp_english      = sp_english,
      scientific_name = sp_scientific,
      spcd_bam        = spcd_bam,
      bam_common_name = bam_common_name,
      raster_path     = bam_raster_path_sp,
      raster_layer    = "mean",
      thresholds      = thresholds,
      hex_summary     = hex_summary
    )
    
    saveRDS(species_output, hex_summary_path)
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    tibble(
      spcd_bam = spcd_bam,
      sp_english = sp_english,
      status = "success",
      elapsed_min = elapsed,
      error = NA_character_
    )
    
  }, error = function(e) {
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    tibble(
      spcd_bam = spcd_bam,
      sp_english = sp_english,
      status = "failed",
      elapsed_min = elapsed,
      error = conditionMessage(e)
    )
  })
  
  gc()
  result
}

# ------------------------------------------------------------------------------
# Run in parallel
# ------------------------------------------------------------------------------

status_path <- file.path(
  OUT_HEX_SUMMARIES_BAM,
  "bam_hex_extraction_status.csv"
)

with_progress({
  p <- progressor(along = species_vec)
  
  extraction_status <- furrr::future_map_dfr(
    species_vec,
    function(spcd) {
      out <- process_bam_species(
        spcd_bam          = spcd,
        species_table     = species_table,
        hex_grid_national = hex_grid_national,
        bam_raster_dir    = BAM_RASTER_PATH,
        out_dir           = OUT_HEX_SUMMARIES_BAM,
        cumpop_breaks     = CUMPOP_BREAKS
      )
      
      p(message = spcd)
      out
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c(
        "tidyverse",
        "terra",
        "sf",
        "exactextractr"
      )
    )
  )
})

readr::write_csv(extraction_status, status_path)

message("Finished BAM extraction.")
message("Successful: ", sum(extraction_status$status == "success"))
message("Skipped:    ", sum(extraction_status$status == "skipped_existing"))
message("Failed:     ", sum(extraction_status$status == "failed"))
message("Status file: ", status_path)

future::plan(future::sequential)

# ------------------------------------------------------------------------------
# Example plot for a single species
# ------------------------------------------------------------------------------

colpal <- c("white", "#FBF7E2", "#CEF2B0", "#18A065", "#006344")

example_result <- readRDS(
  file.path(OUT_HEX_SUMMARIES_BAM, "CANGOO_BAM_hex_summary.rds")
)

ggplot(example_result$hex_summary) +
  geom_sf(aes(fill = category), col = NA) +
  scale_fill_manual(values = colpal, na.value = "gray85", name = "Category") +
  geom_sf(data = region_national, col = "black", fill = "transparent") +
  ggtitle(paste0(example_result$sp_english, " - BAM mean")) +
  theme_bw()