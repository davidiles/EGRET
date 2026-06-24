# ==============================================================================
# Extract predictions from eBird seasonal abundance rasters
# ==============================================================================

library(tidyverse)
library(terra)
library(ebirdst)
library(rnaturalearth)
library(sf)
library(lwgeom)
library(lubridate)
library(exactextractr)
library(readxl)
library(future)
library(furrr)
library(progressr)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_HEX_SUMMARIES_EBIRD, recursive = TRUE, showWarnings = FALSE)
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
  "ebird_species_code" %in% names(species_table)
)

species_table <- species_table |>
  dplyr::filter(!is.na(ebird_species_code)) |>
  dplyr::distinct(ebird_species_code, .keep_all = TRUE)

species_vec <- species_table$ebird_species_code

message("Species with eBird products: ", length(species_vec))

# ------------------------------------------------------------------------------
# National study area
# ------------------------------------------------------------------------------

region_national_path <- file.path(
  OUT_CACHE,
  "region_national.rds"
)

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

process_ebird_species <- function(spcd_ebird,
                                  species_table,
                                  hex_grid_national,
                                  out_dir,
                                  cumpop_breaks) {
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(terra)
    library(ebirdst)
    library(sf)
    library(exactextractr)
  })
  
  sp_row <- species_table |>
    dplyr::filter(ebird_species_code == spcd_ebird) |>
    dplyr::slice(1)
  
  sp_english <- sp_row$english_name
  
  hex_summary_path <- file.path(
    out_dir,
    paste0(spcd_ebird, "_eBird_hex_summary.rds")
  )
  
  
  if (file.exists(hex_summary_path)) {
    return(tibble(
      spcd_ebird = spcd_ebird,
      sp_english = sp_english,
      status = "skipped_existing",
      season = NA_character_,
      elapsed_min = 0,
      error = NA_character_
    ))
  }
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    
    # Optional: separate terra temp files by worker/process
    terra_tmp <- file.path(tempdir(), paste0("terra_", Sys.getpid()))
    dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
    terra::terraOptions(tempdir = terra_tmp)
    
    rast_ebird <- ebirdst::load_raster(
      species    = spcd_ebird,
      product    = "abundance",
      period     = "seasonal",
      metric     = "mean",
      resolution = "3km"
    )
    
    if ("breeding" %in% names(rast_ebird)) {
      rast_ebird <- rast_ebird[["breeding"]]
      selected_layer <- "breeding"
    } else if ("resident" %in% names(rast_ebird)) {
      rast_ebird <- rast_ebird[["resident"]]
      selected_layer <- "resident"
    } else {
      stop("No breeding or resident seasonal layer found.")
    }
    
    hex_summary <- extract_hex_predictions(
      hex_grid     = hex_grid_national,
      rast         = rast_ebird,
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
      sp_english  = sp_english,
      spcd_ebird  = spcd_ebird,
      season      = selected_layer,
      thresholds  = thresholds,
      hex_summary = hex_summary
    )
    
    saveRDS(species_output, hex_summary_path)
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    tibble(
      spcd_ebird  = spcd_ebird,
      sp_english  = sp_english,
      status      = "success",
      season      = selected_layer,
      elapsed_min = elapsed,
      error       = NA_character_
    )
    
  }, error = function(e) {
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    tibble(
      spcd_ebird = spcd_ebird,
      sp_english = sp_english,
      status = "failed",
      season = NA_character_,
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
  OUT_HEX_SUMMARIES_EBIRD,
  "ebird_hex_extraction_status.csv"
)

with_progress({
  p <- progressor(along = species_vec)
  
  extraction_status <- furrr::future_map_dfr(
    species_vec,
    function(spcd) {
      out <- process_ebird_species(
        spcd_ebird        = spcd,
        species_table     = species_table,
        hex_grid_national = hex_grid_national,
        out_dir           = OUT_HEX_SUMMARIES_EBIRD,
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
        "ebirdst",
        "sf",
        "exactextractr"
      )
    )
  )
})

readr::write_csv(extraction_status, status_path)

message("Finished eBird extraction.")
message("Successful: ", sum(extraction_status$status == "success"))
message("Skipped:    ", sum(extraction_status$status == "skipped_existing"))
message("Failed:     ", sum(extraction_status$status == "failed"))
message("Status file: ", status_path)

future::plan(future::sequential)



# ------------------------------------------------------------------------------
# Example plot for a single species
# ------------------------------------------------------------------------------

colpal <-  c("white","#FBF7E2", "#CEF2B0", "#18A065", "#006344")
example_result <- readRDS(
  file.path(OUT_HEX_SUMMARIES_EBIRD, "cacgoo1_eBird_hex_summary.rds")
)

ggplot(example_result$hex_summary) +
  geom_sf(aes(fill = category), col = NA) +
  scale_fill_manual(values = colpal, na.value = "gray85", name = "Category") +
  geom_sf(data = region_national, col = "black", fill = "transparent") +
  ggtitle(paste0(example_result$sp_english, " - ", example_result$season)) +
  theme_bw()
