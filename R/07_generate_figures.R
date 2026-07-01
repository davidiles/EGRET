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
library(ggspatial)
library(patchwork)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

# ------------------------------------------------------------------------------
# Path to save figures
# ------------------------------------------------------------------------------

out_path_figs <- file.path(
  OUT_ROOT,
  "maps"
)

dir.create(out_path_figs, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Species list
# ------------------------------------------------------------------------------

species_lookup_path <- file.path(
  OUT_CACHE,
  "species_lookup_ebird_bam.rds"
)

species_table <- readRDS(species_lookup_path)


species_table <- species_table |>
  dplyr::filter(!is.na(ebird_species_code)) |>
  dplyr::distinct(ebird_species_code, .keep_all = TRUE)

species_vec <- species_table$ebird_species_code

message("Species with eBird products: ", length(species_vec))

bam_hex_path <- file.path(OUT_CACHE, "bam_sf_with_hex_id.rds")


# ------------------------------------------------------------------------------
# National study area
# ------------------------------------------------------------------------------

region_national_path <- file.path(
  OUT_CACHE,
  "region_national.gpkg"
)

region_national <- st_read(region_national_path)

# ------------------------------------------------------------------------------
# Hex grid
# ------------------------------------------------------------------------------

hex_grid_path <- file.path(
  OUT_CACHE,
  paste0("hex_grid_national_", HEX_WIDTH_KM, "_km.gpkg")
)

hex_grid_national <- st_read(hex_grid_path)


head(hex_grid_national)

# ------------------------------------------------------------------------------
# Load processed data
# ------------------------------------------------------------------------------

hex_estimates_path <- file.path(
  OUT_CACHE,
  "hex_species_abundance_categories_ebird_bam.csv")

hex_estimates <- read.csv(hex_estimates_path)

# ------------------------------------------------------------------------------
# Generate figures for every species
# ------------------------------------------------------------------------------

colpal <- c("gray85", "white", "#FBF7E2", "#CEF2B0", "#18A065", "#006344")

# CVD-safe binary palette for BirdLife in/out-of-range panel (deuteranomaly-safe:
# grey vs. blue rather than red/green).
colpal_birdlife <- c("Outside range" = "gray90", "In range" = "#0072B2")

levels_order <- c(
  "Not Modeled",
  "Effectively Absent",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

species_list <- unique(hex_estimates$english_name)

# Optional: simplify boundary only
region_national_plot <- st_simplify(region_national, dTolerance = 1000)

make_category_polygons <- function(x, category_col) {
  x |>
    select(category = {{ category_col }}, geometry) |>
    group_by(category) |>
    summarise(geometry = st_union(geometry), .groups = "drop")
}

for (species_name in species_list) {
  
  print(species_name)
  
  species_df <- hex_estimates %>%
    filter(english_name == species_name)
  
  season_name <- species_df$season[1]
  
  out_file <- file.path(
    out_path_figs,
    paste0(
      make.names(species_name),
      "_",
      season_name,
      "_eBird_BAM.pdf"
    )
  )
  
  if (file.exists(out_file)) next
  
  species_sf <- hex_grid_national %>%
    left_join(species_df, by = "hex_id") |>
    mutate(
      ebird_category = factor(
        coalesce(ebird_category, "Not Modeled"),
        levels = levels_order,
        ordered = TRUE
      ),
      bam_category = factor(
        coalesce(bam_category, "Not Modeled"),
        levels = levels_order,
        ordered = TRUE
      ),
      range_BirdLife = factor(
        coalesce(range_BirdLife, "Outside range"),
        levels = c("Outside range", "In range")
      )
    ) %>%
    dplyr::rename(geometry = geom)
  
  # Dissolve hexes by category before plotting
  ebird_sf_plot    <- make_category_polygons(species_sf, ebird_category)
  bam_sf_plot      <- make_category_polygons(species_sf, bam_category)
  birdlife_sf_plot <- make_category_polygons(species_sf, range_BirdLife)
  
  base_map_theme <- theme_bw() +
    theme(
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.8),
        colour = "grey50"
      ),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  ebird_plot <- ggplot() +
    geom_sf(data = ebird_sf_plot, aes(fill = category), col = NA) +
    geom_sf(
      data = region_national,
      colour = scales::alpha("black", 0.4),
      fill = NA,
      linewidth = 0.1
    )+
    annotation_scale(location = "bl", width_hint = 0.25) +
    scale_fill_manual(values = colpal, drop = FALSE, name = "Category") +
    labs(title = "eBird") +
    base_map_theme
  
  bam_plot <- ggplot() +
    geom_sf(data = bam_sf_plot, aes(fill = category), col = NA) +
    geom_sf(
      data = region_national,
      colour = scales::alpha("black", 0.4),
      fill = NA,
      linewidth = 0.1
    )+
    annotation_scale(location = "bl", width_hint = 0.25) +
    scale_fill_manual(values = colpal, drop = FALSE, name = "Category") +
    labs(title = "BAM") +
    base_map_theme
  
  birdlife_plot <- ggplot() +
    geom_sf(data = birdlife_sf_plot, aes(fill = category), col = NA) +
    geom_sf(
      data = region_national,
      colour = scales::alpha("black", 0.4),
      fill = NA,
      linewidth = 0.1
    )+
    annotation_scale(location = "bl", width_hint = 0.25) +
    scale_fill_manual(values = colpal_birdlife, drop = FALSE, name = "BirdLife range") +
    labs(title = "BirdLife") +
    base_map_theme
  
  combined_plot <- ebird_plot + bam_plot + birdlife_plot +
    plot_layout(nrow = 1) +
    plot_annotation(
      title = paste0(species_name, " - ", season_name)
    )
  
  ggsave(
    out_file,
    combined_plot,
    device = cairo_pdf,
    width = 28,
    height = 8,
    units = "in",
    bg = "white"
  )
}
