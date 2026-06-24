library(googledrive)
library(dplyr)
library(glue)
library(fs)
library(purrr)
library(BAMexploreR)

# Note: may have to clear cached auths using: drive_deauth()
# Then use:
# drive_auth(
#   email = "your_email_here@gmail.com",
#   scopes = "drive.readonly"
# )

drive_auth(scopes = "drive.readonly")

root <- drive_get(as_id("https://drive.google.com/drive/folders/0AMsqxXlPq2e-Uk9PVA"))

packaged_root <- googledrive::drive_ls(root) |>
  dplyr::filter(.data$name == "output") |>
  googledrive::drive_ls() |>
  dplyr::filter(.data$name == "10_packaged")

download_species_tif <- function(sp,
                                 year = 2020,
                                 stratum = "Canada",
                                 packaged_root,
                                 local_dir = "bam_downloads",
                                 overwrite = FALSE,
                                 verbose = TRUE) {
  
  fs::dir_create(local_dir)
  
  tif_name <- glue::glue("{sp}_{stratum}_{year}.tif")
  out_file <- fs::path(local_dir, tif_name)
  
  # ---- NEW: skip if already exists ----
  if (fs::file_exists(out_file) && !overwrite) {
    if (verbose) message("Skipping (already exists): ", tif_name)
    return(out_file)
  }
  
  # ---- proceed to Drive only if needed ----
  species_folder <- googledrive::drive_ls(packaged_root) |>
    dplyr::filter(.data$name == sp)
  
  if (nrow(species_folder) == 0) {
    warning("Could not find species folder: ", sp)
    return(NA_character_)
  }
  
  stratum_folder <- googledrive::drive_ls(species_folder) |>
    dplyr::filter(.data$name == stratum)
  
  if (nrow(stratum_folder) == 0) {
    warning("Could not find stratum folder '", stratum, "' for species ", sp)
    return(NA_character_)
  }
  
  tif_file <- googledrive::drive_ls(stratum_folder) |>
    dplyr::filter(.data$name == tif_name)
  
  if (nrow(tif_file) == 0) {
    warning("Could not find file: ", tif_name)
    return(NA_character_)
  }
  
  if (verbose) message("Downloading: ", tif_name)
  
  googledrive::drive_download(
    file = tif_file,
    path = out_file,
    overwrite = overwrite
  )
  
  return(out_file)
}

sp_list <- BAMexploreR::spp_tbl$speciesCode

strata_to_download <- "Canada"
year_to_download <- 2020
local_dir <- paste0("C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/BAMv5_model_outputs/",year_to_download)
dir_create(local_dir)

downloaded <- purrr::map_chr(
  sp_list,
  ~ download_species_tif(
    sp = .x,
    year = year_to_download,
    stratum = strata_to_download,
    packaged_root = packaged_root,
    local_dir = local_dir
  )
)
