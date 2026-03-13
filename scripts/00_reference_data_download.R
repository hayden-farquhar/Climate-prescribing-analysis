# ==============================================================================
# Script 00: Download ABS Reference Data and Shapefiles
# ==============================================================================
# Downloads SEIFA 2021, population estimates, remoteness areas, and
# LGA/SA4 boundary shapefiles from the ABS. All URLs are direct downloads
# requiring no authentication.
#
# Outputs:
#   data/reference/SEIFA_2021_LGA.xlsx
#   data/reference/ERP_LGA_2001-2024.xlsx
#   data/spatial/LGA_2021_AUST_GDA2020_SHP.zip  (+ extracted shapefiles)
#   data/spatial/SA4_2021_AUST_SHP_GDA2020.zip   (+ extracted shapefiles)
#   data/spatial/RA_2021_AUST_GDA2020.zip         (+ extracted shapefiles)
#
# Runtime: ~2-5 minutes depending on connection speed
# ==============================================================================

library(here)

# Set project root (adjust if not using {here})
project_dir <- here::here()
ref_dir     <- file.path(project_dir, "data", "reference")
spatial_dir <- file.path(project_dir, "data", "spatial")

dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

# Use a longer timeout for large shapefile downloads
options(timeout = 600)

# ------------------------------------------------------------------------------
# 1. SEIFA 2021 — Index of Relative Socio-economic Disadvantage (IRSD) at LGA
# ------------------------------------------------------------------------------
cat("Downloading SEIFA 2021 (LGA)...\n")
download.file(
  url = "https://www.abs.gov.au/statistics/people/people-and-communities/socio-economic-indexes-areas-seifa-australia/2021/Local%20Government%20Area,%20Indexes,%20SEIFA%202021.xlsx",
  destfile = file.path(ref_dir, "SEIFA_2021_LGA.xlsx"),
  mode = "wb"
)
cat("  -> SEIFA_2021_LGA.xlsx downloaded.\n\n")

# ------------------------------------------------------------------------------
# 2. Estimated Resident Population by LGA, 2001-2024
# ------------------------------------------------------------------------------
cat("Downloading ERP by LGA (2001-2024)...\n")
download.file(
  url = "https://www.abs.gov.au/statistics/people/population/regional-population/2023-24/32180DS0004_2001-24.xlsx",
  destfile = file.path(ref_dir, "ERP_LGA_2001-2024.xlsx"),
  mode = "wb"
)
cat("  -> ERP_LGA_2001-2024.xlsx downloaded.\n\n")

# ------------------------------------------------------------------------------
# 3. LGA Boundary Shapefiles (2021 ASGS Edition 3, GDA2020)
# ------------------------------------------------------------------------------
lga_zip <- file.path(spatial_dir, "LGA_2021_AUST_GDA2020_SHP.zip")
cat("Downloading LGA boundary shapefiles...\n")
download.file(
  url = "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/LGA_2021_AUST_GDA2020_SHP.zip",
  destfile = lga_zip,
  mode = "wb"
)
unzip(lga_zip, exdir = file.path(spatial_dir, "LGA_2021"))
cat("  -> LGA shapefiles extracted to data/spatial/LGA_2021/\n\n")

# ------------------------------------------------------------------------------
# 4. SA4 Boundary Shapefiles (2021 ASGS Edition 3, GDA2020)
# ------------------------------------------------------------------------------
sa4_zip <- file.path(spatial_dir, "SA4_2021_AUST_SHP_GDA2020.zip")
cat("Downloading SA4 boundary shapefiles...\n")
download.file(
  url = "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/SA4_2021_AUST_SHP_GDA2020.zip",
  destfile = sa4_zip,
  mode = "wb"
)
unzip(sa4_zip, exdir = file.path(spatial_dir, "SA4_2021"))
cat("  -> SA4 shapefiles extracted to data/spatial/SA4_2021/\n\n")

# ------------------------------------------------------------------------------
# 5. Remoteness Areas (2021 ASGS Edition 3, GDA2020)
# ------------------------------------------------------------------------------
ra_zip <- file.path(spatial_dir, "RA_2021_AUST_GDA2020.zip")
cat("Downloading Remoteness Areas shapefiles...\n")
download.file(
  url = "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/RA_2021_AUST_GDA2020.zip",
  destfile = ra_zip,
  mode = "wb"
)
unzip(ra_zip, exdir = file.path(spatial_dir, "RA_2021"))
cat("  -> Remoteness Areas extracted to data/spatial/RA_2021/\n\n")

# ------------------------------------------------------------------------------
# Verify downloads
# ------------------------------------------------------------------------------
cat("=== Download Summary ===\n")
ref_files <- list.files(ref_dir, full.names = FALSE)
cat("Reference files:\n")
for (f in ref_files) cat("  ", f, "\n")

spatial_dirs <- list.dirs(spatial_dir, recursive = FALSE, full.names = FALSE)
cat("Spatial directories:\n")
for (d in spatial_dirs) {
  n <- length(list.files(file.path(spatial_dir, d)))
  cat("  ", d, sprintf("(%d files)\n", n))
}

cat("\nDone. All ABS reference data downloaded.\n")
