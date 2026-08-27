# =============================================================================
# Assign each lion sighting a protection status by overlaying it on the
# national park boundary
# =============================================================================
# Reads  data/lion_sighting_points.csv   (written by R/01_capture_histories.R)
# Writes data/sighting_protection.csv    (read back by R/01_capture_histories.R
#                                         on its next run to build area[i,t])
#
# Run order: 01 -> 02 -> 01 again -> 03.  The first pass of 01 has no
# protection file and falls back to a placeholder stratum; the second pass
# picks up this script's output.

library(terra)
library(tmap)
library(readr)
library(dplyr)

source("R/00_paths.R")

zmb <- vect(WDPA_GDB, layer = WDPA_LAYER)
summary(zmb)
luangwa <- zmb[zmb$NAME_ENG %in% c("North Luangwa", "South Luangwa", "Luambe")]

sightings <- read_csv("data/lion_sighting_points.csv")
# State the CRS explicitly rather than relying on the points and the polygons
# happening to share one: an unset CRS makes the intersect silently wrong if
# the WDPA layer is ever supplied projected rather than in lon/lat.
sightings <- vect(sightings, geom = c("lon", "lat"), crs = "EPSG:4326")
luangwa <- project(luangwa, crs(sightings))
plot(sightings)

# Find intersects
idx <- is.related(sightings, luangwa, "intersects")
sightings$protected_area <- idx

tm_shape(luangwa) +
  tm_polygons(col = "lightgray") +
  tm_shape(sightings) +
  tm_dots(col = "protected_area", palette = "Viridis")

sighting_protection <- sightings |>
  as.data.frame() |>
  mutate(
    protected = case_when(
      protected_area == TRUE ~ 1,
      protected_area == FALSE ~ 0
    )
  ) |>
  select(SightingID, protected) |>
  readr::write_csv("data/sighting_protection.csv")
