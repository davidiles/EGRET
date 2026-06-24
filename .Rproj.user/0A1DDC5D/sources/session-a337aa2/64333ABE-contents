# ==============================================================================
# Create species lookup table for avian core, eBird, and BAM
# ==============================================================================

library(tidyverse)
library(readxl)
library(ebirdst)
library(BAMexploreR)

rm(list = ls())

source("R/config.R")
source("R/functions_new.R")

dir.create(OUT_CACHE, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Output paths
# ------------------------------------------------------------------------------

species_lookup_path <- file.path(
  OUT_CACHE,
  "species_lookup_ebird_bam.rds"
)

species_lookup_csv_path <- file.path(
  OUT_CACHE,
  "species_lookup_ebird_bam.csv"
)

species_lookup_flagged_path <- file.path(
  OUT_CACHE,
  "species_lookup_ebird_bam_flagged.csv"
)

# ------------------------------------------------------------------------------
# 4-letter bird codes from American Ornithological Society
# Downloaded from: https://www.birdpop.org/pages/birdSpeciesCodes.php
# ------------------------------------------------------------------------------
AOS <- read.csv(AOS_PATH) %>%
  dplyr::select(COMMONNAME,SCINAME,SPEC) %>%
  dplyr::rename(AOS_CODE_4L = SPEC) %>%
  distinct()

AOS_sci <- AOS |>
  transmute(
    scientific_name = SCINAME,
    AOS_CODE_4L = AOS_CODE_4L,
    AOS_common_name = COMMONNAME
  ) |>
  distinct(scientific_name, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Avian core species list
# ------------------------------------------------------------------------------

avian_core <- read_xlsx(AVIAN_CORE_PATH) |>
  filter(
    Taxon == "Aves",
    Full_Species__Espèce_complète == "Yes - Oui",
    CDN_Status__Statut_CDN %in% c("BRE", "BRE_OCC", "RNB")
  ) |>
  transmute(
    english_name    = English_Name__Nom_Anglais,
    scientific_name = Scientific_Name__Nom_Scientifique,
    bbl_code        = Alpha_Code_BBL__Code_Alpha_BBL
  ) |>
  distinct(scientific_name, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# eBird lookup
# ------------------------------------------------------------------------------

ebird_lookup <- ebirdst::ebirdst_runs |>
  filter(species_code != "yebsap-example") |>
  transmute(
    scientific_name,
    ebird_species_code    = species_code,
    ebird_common_name     = common_name,
    ebird_status_version  = status_version_year
  ) |>
  arrange(scientific_name, desc(ebird_status_version)) |>
  distinct(scientific_name, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# BAM lookup
# ------------------------------------------------------------------------------

bam_lookup <- BAMexploreR::spp_tbl |>
  transmute(
    scientific_name  = scientificName,
    bam_species_code = speciesCode,
    bam_common_name  = commonName
  ) |>
  
  # Correct a spelling mistake in the BAM table
  mutate(
    scientific_name = recode(
      scientific_name,
      "Leiothylpis ruficapilla" = "Leiothlypis ruficapilla"
    )
  ) |>
  distinct(scientific_name, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Combined lookup
# ------------------------------------------------------------------------------

species_lookup <- avian_core |>
  left_join(AOS_sci, by = "scientific_name") |>
  left_join(ebird_lookup, by = "scientific_name") |>
  left_join(bam_lookup, by = "scientific_name") |>
  mutate(
    in_aos   = !is.na(AOS_CODE_4L),
    in_ebird = !is.na(ebird_species_code),
    in_bam   = !is.na(bam_species_code),
    in_final_table = in_ebird & in_bam,
    missing_from = case_when(
      in_ebird & in_bam   ~ NA_character_,
      !in_ebird & !in_bam ~ "eBird and BAM",
      !in_ebird           ~ "eBird",
      !in_bam             ~ "BAM"
    ),
    name_mismatch_flag = case_when(
      in_ebird & !is.na(ebird_common_name) &
        str_to_lower(english_name) != str_to_lower(ebird_common_name) ~ TRUE,
      in_bam & !is.na(bam_common_name) &
        str_to_lower(english_name) != str_to_lower(bam_common_name) ~ TRUE,
      !is.na(AOS_common_name) &
        str_to_lower(english_name) != str_to_lower(AOS_common_name) ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  arrange(english_name)

# ------------------------------------------------------------------------------
# Species included in both products
# ------------------------------------------------------------------------------

species_lookup_final <- species_lookup |>
  filter(in_final_table)

# ------------------------------------------------------------------------------
# Species needing review
# ------------------------------------------------------------------------------

species_lookup_flagged <- species_lookup |>
  filter(!in_final_table | name_mismatch_flag)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

saveRDS(species_lookup, species_lookup_path)

readr::write_csv(species_lookup, species_lookup_csv_path)
readr::write_csv(species_lookup_flagged, species_lookup_flagged_path)

message("Total avian-core species: ", nrow(species_lookup))
message("Included in both eBird and BAM: ", sum(species_lookup$in_final_table))
message("Missing from eBird: ", sum(!species_lookup$in_ebird))
message("Missing from BAM: ", sum(!species_lookup$in_bam))
message("Name mismatch flags: ", sum(species_lookup$name_mismatch_flag))
message("Missing from AOS: ", sum(!species_lookup$in_aos))

message("Saved full lookup: ", species_lookup_csv_path)
message("Saved flagged lookup: ", species_lookup_flagged_path)

a = species_lookup %>%
  filter(!is.na(bam_species_code))

subset(a,AOS_CODE_4L != bam_species_code)
