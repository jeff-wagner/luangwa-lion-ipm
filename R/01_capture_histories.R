# KNP WD CMR analysis
# processing tables from Access to create input files specifically
# for annual Bayesian closed capture models with adata augmentation

rm(list = ls())

# Install any required packages that are missing, then load them
required_packages <- c(
  "dplyr",
  "lubridate",
  "reshape",
  "tidyr",
  "RODBC",
  "stringr",
  "IPMbook",
  "purrr",
  "openxlsx",
  "terra"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

lapply(required_packages, library, character.only = TRUE)

source("R/00_paths.R")

#read the relevant tables from Access
#alternatively could use the dbi and dbplyr packages to connect to the Access
#database directly, but if the files fit in memory that would be slower
#and have no advantages, see: https://db.rstudio.com/dplyr/

# Open DB connection
dbCon <- odbcConnectAccess2007(ACCESS_DB)
#sqlTables(dbCon)
#lionSight <- sqlQuery(dbCon, "
# SELECT
#   AnimSight.*,
#   Sight.*,
#   Anim.*
# FROM (tblAnimalsatSighting AS AnimSight
#   LEFT JOIN tblSightings AS Sight
#     ON AnimSight.SightingID = Sight.SightingID)
#   LEFT JOIN tblAnimals AS Anim
#     ON AnimSight.AnimalID = Anim.AnimalID
# WHERE Sight.LionSighting = TRUE
# ")

lionSight <- sqlQuery(
  dbCon,
  "
                      SELECT
                        AnimSight.*,
                        Sight.*,
                        Anim.*,
                        [Group].*
                      FROM ((tblAnimalsatSighting AS AnimSight
                        LEFT JOIN tblSightings AS Sight 
                          ON AnimSight.SightingID = Sight.SightingID)
                        LEFT JOIN tblAnimals AS Anim 
                          ON AnimSight.AnimalID = Anim.AnimalID)
                        LEFT JOIN tblAnimalSightings AS [Group] 
                          ON AnimSight.SightingID = [Group].SightingID
                      WHERE Sight.LionSighting = TRUE
                      "
)

# Close DB connection
odbcClose(dbCon)

#note that this does not eliminate some rows for other species that were present at dog sightings
#but will eliminate those cases at final join of capture histories with covariates

#convert date from character to date format (m/d/yyyy)
#NOTE THAT THIS STEP HAS NOT BEEN CONSISTENT BETWEEN IMPORTS FROM ACCESS DEPENDING ON THE COMPUTER USED
#tweak the format argument as needed
lionSight$SightingDate <- as.Date.character(
  x = lionSight$SightingDate,
  format = "%Y-%m-%d"
)
lionSight$DOB <- as.Date.character(x = lionSight$DOB, format = "%Y-%m-%d")
# QC for incorrect dates, AnimalIDs and species
lionSight %>%
  filter(lubridate::year(SightingDate) < 1980, lubridate::year(DOB) < 1980)
# Looking in database, seems like the one odd date is a typo and (maybe should be 2017-06-17?) but lion was seen in the same month again so just delete this sighting (should have not effect on inference)
lionSight %>%
  filter(AnimalID == "ELI-126") %>%
  pull(SightingDate, SightingID)
# delete the problematic sighting
lionSight <- lionSight[!lionSight$SightingID %in% "E-18435", ]


# Remove NA AnimalIDs
lionSight <- lionSight %>%
  filter(!is.na(AnimalID))

# Fix case in AnimalIDs
has_lowercase <- str_detect(lionSight$AnimalID, "[[:lower:]]")
lionSight$AnimalID[has_lowercase] <- toupper(lionSight$AnimalID[has_lowercase])
length(unique(lionSight$AnimalID))

# Remove other species (some wild dog obs in there)
lionSight <- lionSight %>%
  filter(Species == "Lion")
#table(unique(lionSight$AnimalID) %in% animal.ids)

# filter out GPS collar data
table(lionSight$`How found`)

hist(
  lionSight %>%
    filter(`How found` == "Satellite") %>%
    pull(SightingDate),
  breaks = 20
)


lionSight <- lionSight %>%
  filter(!`How found` %in% c("Collar Data"))

#some lion ID's are errors. Not all may be findable.
#some lions were 'seen' before their date of birth. these are obvious errors
sum(lionSight$SightingDate < lionSight$DOB, na.rm = T)
problemdates <- lionSight %>%
  # Filter for sightings outside the valid lifespan
  filter(SightingDate > `Date of Death` | SightingDate <= DOB) %>%
  select(
    AnimalID,
    SightingID,
    SightingDate,
    DOB,
    `Date of Death`,
    `Date of Death error`
  ) %>%
  arrange(AnimalID, SightingDate) %>%
  mutate(
    beforebirth = SightingDate <= DOB,
    afterdeath = SightingDate > `Date of Death`
  )
problemdates

# example of hashing out issues with lion ELI-647
problemlions <- lionSight %>%
  # Filter for sightings outside the valid lifespan
  filter(AnimalID == "ELI-143") %>%
  select(
    AnimalID,
    SightingID,
    SightingDate,
    DOB,
    `Date of Death`,
    `Date of Death error`,
    AnimalGroup
  ) %>%
  arrange(AnimalID, SightingDate) %>%
  mutate(
    beforebirth = SightingDate <= DOB,
    afterdeath = SightingDate > `Date of Death`
  )
problemlions

#View(lionSight[lionSight$AnimalGroup %in% c("Mwamba II")&lionSight$DOB %in% c("2013-01-31"),])
#View(lionSight[lionSight$AnimalGroup %in% c("Mwamba II")&year(lionSight$SightingDate) %in% c(2016),])
#View(lionSight[lionSight$SightingDate %in% c("2016-08-21"),])
#View(lionSight[lionSight$SightingID %in% c("E-15549"),])

#lionSight %>%
#  filter(AnimalID=="ELI-173") %>%
#  pull(SightingDate,SightingID)

### Corrections - see Luangwa DB CHange Log for details
### Once permanent changes to DB are made, these corrections should not be necessary

lionSight$`Date of Death`[lionSight$AnimalID %in% c("ELI-647")] <- NULL # wrong lion cub recorded as dead in infanticide, but still being sighted

lionSight$`Date of Death`[lionSight$AnimalID %in% c("ELI-653")] <- "2013-08-06" # lion cub recorded as surviving an infanticide, but really died

# add sighting of lion at is collaring
newsightingELI138 <- lionSight[lionSight$AnimalID %in% c("ELI-138"), ][1, ] # lion was collared in 2013 but not recorded as at immobilization(sighting E-15886)
newsightingELI138[,] <- NA
# newsightingELI138[, c(
#   "AnimalID",
#   "SightingID",
#   "AnimalGroup",
#   "ObservationStrategy",
#   "LionSighting",
#   "Immobilisation",
#   "SightingDate",
#   "Time",
#   "GPS"
# )] <- c(
#   "ELI-138",
#   "E-15886",
#   "Mwamba I",
#   "Collar: Opportunistic",
#   "1",
#   "1",
#   "2013-12-05",
#   "13:45",
#   "MSB"
# )
# lionSight <- rbind(newsightingELI138, lionSight)

# The above code is breaking the date column, use add_row() instead
lionSight <- lionSight |>
  add_row(
    "AnimalID" = "ELI-138",
    "SightingID" = "E-15886",
    "AnimalGroup" = "Mwamba I",
    "ObservationStrategy" = "Collar: Opportunistic",
    "LionSighting" = 1,
    "Immobilisation" = 1,
    "SightingDate" = as.Date("2013-12-05"),
    "Time" = as.POSIXct("1899-12-30 13:45"),
    "GPS" = "MSB"
  )

lionSight$`Date of Death`[lionSight$AnimalID %in% c("ELI-122")] <- NA #recorded as dead but seen again in same group, not observed dead? wrong cub?

lionSight$`Date of Death`[lionSight$AnimalID %in% c("ELI-119")] <- NA #recorded as dead but seen again in same group, not observed dead?

lionSight$`Date of Death`[lionSight$AnimalID %in% c("ELI-169")] <- NA #recorded as dead but seen again in same group, not observed dead?

# lionSight <- lionSight[
#   -which(
#     lionSight$AnimalID %in% c("ELI-173") & lionSight$SightingDate > "2018-08-28"
#   ),
# ] # unfixable misidentification of a lion, this lion sighted after its observed death hunting an elephant calf.  These must be sightings of some other lion.
# Fixing the above code as it was getting rid of all rows in lionSight, though this animal ID doesn't seem to be present?
lionSight <- lionSight |>
  filter_out(AnimalID %in% c("ELI-173") & SightingDate > "2018-08-28") # Removed 0 rows?

lionSight$SightingDate[lionSight$SightingID %in% c("E-17254")] <- "2017-07-24" # year typo in sighting date in DB


##############################################################################

# add a detected variable, all ones at this point for the sightings themselves
lionSight$detect <- rep(1, length(lionSight$SightingID))


# Get recoded, canonical identifers for sequential time, base on annual, monthly, or two-month intervals
# note that if some intervals can be discardede.g. jan/feb, then discard those intervals and record using as.numeric(as.factor()) to restore a canonical id.
lionSight <- lionSight %>%
  mutate(
    year = year(SightingDate),
    annualid = year - min(year) + 1,
    month = month(SightingDate),
    monthlyid = as.numeric(as.factor(year * 10 + month)),
    twomonth = (month - 1) %/% 2 + 1, # records time of year into 6, 2-month periods, jan-feb, mar-apr,...
    twomonthlyid = as.numeric(as.factor(year * 10 + twomonth))
  )

#####
#####
#####
#made with help from chatgpt
# df must have:
#   - an individual id column (e.g., id)
#   - a date column (e.g., date)
# note that each month of the year is represented, even if no data were collected

monthly_tally <- lionSight %>%
  mutate(
    month = floor_date(as.Date(SightingDate), unit = "month")
  ) %>%
  count(AnimalID, month, name = "n_sightings") %>%
  complete(
    AnimalID,
    month = seq.Date(
      floor_date(min(month), "year"),
      floor_date(max(month), "year") + months(11),
      by = "month"
    ),
    fill = list(n_sightings = 0)
  ) %>%
  arrange(AnimalID, month)

head(monthly_tally)

# Wide table with one column per month (chronological)
monthly_tally_wide <- monthly_tally %>%
  mutate(month_label = format(month, "%Y-%m")) %>%
  select(AnimalID, month_label, n_sightings) %>%
  pivot_wider(
    names_from = month_label,
    values_from = n_sightings,
    values_fill = 0
  ) %>%
  arrange(AnimalID)

# Bimonthly (Jan-Feb, Mar-Apr, ..., Nov-Dec) tallies for each individual
all_periods <- expand.grid(
  year = seq(
    min(year(as.Date(lionSight$SightingDate))),
    max(year(as.Date(lionSight$SightingDate)))
  ),
  bimonth = 1:6
)

bimonthly_tally <- lionSight %>%
  mutate(
    year = year(as.Date(SightingDate)),
    bimonth = ceiling(month(as.Date(SightingDate)) / 2)
  ) %>%
  count(AnimalID, year, bimonth, name = "n_sightings") %>%
  complete(
    AnimalID,
    nesting(year = all_periods$year, bimonth = all_periods$bimonth),
    fill = list(n_sightings = 0)
  ) %>%
  mutate(
    period_label = paste0(year, "-P", bimonth)
  ) %>%
  arrange(AnimalID, year, bimonth)

bimonthly_tally_wide <- bimonthly_tally %>%
  select(AnimalID, period_label, n_sightings) %>%
  pivot_wider(
    names_from = period_label,
    values_from = n_sightings,
    values_fill = 0
  ) %>%
  arrange(AnimalID)

head(bimonthly_tally_wide)

#####
# Data from first Thandiwe Mweetwa PLOS paper takes precedent and should overwrite all captures histories for Luangwa 2008-2015
# additionally, sex and age data for lions 2008-2015 should come from the PLOS one data set given the time put into editting it by TM
# original files
# script created with Claude sonnet 4.6
# folder <- "C:/Users/christianson/OneDrive - University of Wyoming/Jeff Wagner's files - Luangwa Lion/Code/first luangwa paper/2008-2015excel"
folder <- file.path(FIRST_PAPER_DIR, "2008-2015excel")
month_cols <- c(
  "jan",
  "feb",
  "mar",
  "apr",
  "may",
  "jun",
  "jul",
  "aug",
  "sep",
  "oct",
  "nov",
  "dec"
)

read_lion_TMPLOS <- function(path) {
  yr <- as.integer(str_extract(basename(path), "\\d{4}"))
  raw <- openxlsx::read.xlsx(path)

  names(raw) <- tolower(names(raw))
  names(raw)[names(raw) == "sept"] <- "sep"

  if (yr %in% 2013:2014) {
    raw <- raw[-1, ]
    names(raw)[names(raw) == "x2"] <- "AnimalID"
  } else if (yr == 2015) {
    names(raw)[names(raw) == "lionid"] <- "AnimalID"
  } else {
    names(raw)[names(raw) == "id"] <- "AnimalID"
  }

  raw |>
    select(AnimalID, all_of(month_cols)) |>
    filter(!is.na(AnimalID), AnimalID != "ID") |>
    mutate(
      across(all_of(month_cols), ~ as.integer(replace_na(as.numeric(.x), 0))),
      year = yr
    )
}

all_long_TMPLOS <- list.files(
  folder,
  pattern = "\\.xlsx$",
  full.names = TRUE
) |>
  map(read_lion_TMPLOS) |>
  list_rbind() |>
  pivot_longer(
    all_of(month_cols),
    names_to = "month_name",
    values_to = "detected"
  ) |>
  mutate(month_num = match(month_name, tolower(month.abb)))

all_animals_TMPLOS <- sort(unique(all_long_TMPLOS$AnimalID))

monthly_tally_TMPLOS <- all_long_TMPLOS |>
  tidyr::complete(
    AnimalID = all_animals_TMPLOS,
    year = 2008:2015,
    month_num = 1:12,
    fill = list(detected = 0L)
  ) |>
  mutate(period_label = paste0(year, "-", sprintf("%02d", month_num))) |>
  arrange(AnimalID, year, month_num)

monthly_tally_wide_TMPLOS <- monthly_tally_TMPLOS |>
  select(AnimalID, period_label, detected) |>
  tidyr::pivot_wider(names_from = period_label, values_from = detected) |>
  mutate(across(where(is.integer), ~ replace_na(.x, 0L))) |>
  arrange(AnimalID)

bimonthly_tally_TMPLOS <- monthly_tally_TMPLOS |>
  mutate(
    bimonth = ceiling(month_num / 2),
    bp_label = paste0(year, "-P", bimonth)
  ) |>
  group_by(AnimalID, year, bimonth, bp_label) |>
  summarise(detected = sum(detected, na.rm = TRUE), .groups = "drop") |>
  arrange(AnimalID, year, bimonth)

bimonthly_tally_wide_TMPLOS <- bimonthly_tally_TMPLOS |>
  select(AnimalID, bp_label, detected) |>
  tidyr::pivot_wider(names_from = bp_label, values_from = detected) |>
  mutate(across(where(is.integer), ~ replace_na(.x, 0L))) |>
  arrange(AnimalID)
#####
# Find ID's in TM's PLOS data that are not in the current Access Database
monthly_tally_wide_TMPLOS$AnimalID[
  !monthly_tally_wide_TMPLOS$AnimalID %in% monthly_tally_wide$AnimalID
]


(c("ELI-131") %in% lionSight$AnimalID)


##############################################################################
##############################################################################
##
## FINAL FORMATTING: build the data bundles consumed by the IPM
## (structure follows Creel et al. 2024, Conserv Sci Pract, and the JAGS code
##  in R_JAGS_IPM_CODE.docx)
##
## The IPM's survival sub-model needs:
##     y[i,t]     detection (0/1) of individual i on occasion t
##     f[i]       occasion of first detection
##     sex[i]     1 = female, 2 = male   (beta.sex is the "male effect")
##     age[i,t]   age class 1-4 of individual i on occasion t
##                (cub 0-1.99, subadult 2-3.99, prime adult 4-5.99, old 6+)
##     area[i,t]  protection/region stratum of individual i on occasion t
##     nind, n.occasions, nageclass
##
## The closed-capture sub-model needs one capture history matrix per year
## (chc), fit OUTSIDE the IPM; its posterior medians become Pop[t].
##
##############################################################################
##############################################################################

# ---------------------------------------------------------------------------
# 0. Analysis settings -- change these, not the code below
# ---------------------------------------------------------------------------

STUDY_YEARS <- 2008:2023

# Creel et al. used three 2-month bins per year covering May-Oct.  These data
# run a little later in the year -- individual-months per bimonth are Jan-Feb
# 184, Mar-Apr 526, May-Jun 1381, Jul-Aug 1946, Sep-Oct 2100, Nov-Dec 1131 --
# so a fourth Nov-Dec bin is included, raising retention from ~75% to ~90%.
# Everything downstream, including the age matrix, adapts to whatever is set
# here; the JAGS model annualises survival as phi^n_bins, so n_bins must be
# changed there too.
BIN_MONTHS <- list(
  "1" = 5:6, # May-Jun
  "2" = 7:8, # Jul-Aug
  "3" = 9:10, # Sep-Oct
  "4" = 11:12 # Nov-Dec
)

# Years for which Thandiwe Mweetwa's PLOS ONE capture histories take precedence
PLOS_YEARS <- 2008:2015

# Age-class boundaries (years); 4 classes as in the paper
get_age_class <- function(age_years) {
  cut(
    age_years,
    breaks = c(-Inf, 2, 4, 6, Inf),
    labels = FALSE,
    right = FALSE
  )
}

set.seed(20260820) # unknown sexes are drawn at random -- keep this reproducible
# ---------------------------------------------------------------------------
# 1. Combine Access and PLOS detection records; PLOS wins for 2008-2015
# ---------------------------------------------------------------------------
# Both sources are reduced to long, binary, individual-by-month records first.
# Merging in long form (rather than aligning wide matrices) means an individual
# present in only one source cannot be silently mis-aligned.
#
# Precedence is applied PER LION-YEAR, not per lion.  This matters: the
# 2008-2012 spreadsheets are the curated PLOS histories and cover 78-116 lions
# each, closely matching Access.  The 2013-2015 files are smaller working
# files with a different column layout (45-53 lions against 84-108 in Access).
# Overriding a lion's whole 2008-2015 block whenever it appears in ANY PLOS
# file therefore blanks out real 2013-2015 detections for lions those later
# files simply do not list -- 767 individual-months deleted, against 219 when
# precedence is scoped to the lion-years actually present in a source file.
# The year-by-year table printed below shows the effect; check it after any
# change to the source spreadsheets.

access_long <- monthly_tally |>
  transmute(
    AnimalID,
    year = year(month),
    month_num = month(month),
    detected = as.integer(n_sightings > 0)
  )

# "ELI-NC" is not a real mark ID -- drop it
plos_long <- monthly_tally_TMPLOS |>
  filter(AnimalID != "ELI-NC") |>
  transmute(AnimalID, year, month_num, detected = as.integer(detected > 0))

# The lion-years actually listed in a PLOS spreadsheet.  all_long_TMPLOS is
# read straight off the files, before tidyr::complete() pads every lion out to
# all eight years -- that padding is what makes per-lion precedence destructive.
plos_lion_years <- all_long_TMPLOS |>
  filter(AnimalID != "ELI-NC", year %in% PLOS_YEARS) |>
  distinct(AnimalID, year)

plos_long <- plos_long |>
  semi_join(plos_lion_years, by = c("AnimalID", "year"))

combined_long <- bind_rows(
  access_long |>
    anti_join(plos_lion_years, by = c("AnimalID", "year")) |>
    mutate(source = "access"),
  plos_long |>
    mutate(source = "plos")
)

# --- what did the override actually change? --------------------------------
overwrite_check <- full_join(
  access_long |>
    semi_join(plos_lion_years, by = c("AnimalID", "year")) |>
    dplyr::rename(access = detected),
  plos_long |> dplyr::rename(plos = detected),
  by = c("AnimalID", "year", "month_num")
) |>
  mutate(across(c(access, plos), ~ replace_na(.x, 0L)))

cat("\n--- PLOS override ---\n")
cat("lion-years overridden:          ", nrow(plos_lion_years), "\n")
cat(
  "lions in PLOS but not Access:   ",
  sum(!unique(plos_long$AnimalID) %in% access_long$AnimalID),
  "\n"
)
cat(
  "individual-months gained (0->1):",
  sum(overwrite_check$access == 0 & overwrite_check$plos == 1),
  "\n"
)
cat(
  "individual-months lost   (1->0):",
  sum(overwrite_check$access == 1 & overwrite_check$plos == 0),
  "\n"
)
cat(
  "individual-months agreeing:     ",
  sum(overwrite_check$access == overwrite_check$plos),
  "\n"
)
cat("by year (watch for years where 'lost' dwarfs 'gained'):\n")
print(
  overwrite_check |>
    group_by(year) |>
    summarise(
      lions = n_distinct(AnimalID),
      lost = sum(access == 1 & plos == 0),
      gained = sum(access == 0 & plos == 1),
      n_access = sum(access),
      n_plos = sum(plos),
      .groups = "drop"
    ) |>
    as.data.frame(),
  row.names = FALSE
)


# ---------------------------------------------------------------------------
# 2. Collapse to the CJS occasion structure
# ---------------------------------------------------------------------------

n_bins <- length(BIN_MONTHS)
n_occasions <- length(STUDY_YEARS) * n_bins

bin_lookup <- tibble(
  month_num = unlist(BIN_MONTHS, use.names = FALSE),
  bin = rep(seq_len(n_bins), lengths(BIN_MONTHS))
)

# occasion key: which year/bin each of the 1..n_occasions columns refers to,
# plus the midpoint date of the bin, which is used below to age individuals
occasion_key <- tidyr::expand_grid(
  year = STUDY_YEARS,
  bin = seq_len(n_bins)
) |>
  mutate(
    occasion = row_number(),
    label = paste0(year, "-B", bin),
    bin_start = as.Date(paste0(
      year,
      "-",
      vapply(bin, function(b) BIN_MONTHS[[b]][1], 0),
      "-01"
    )),
    mid_date = bin_start +
      floor(30.44 * vapply(bin, function(b) length(BIN_MONTHS[[b]]), 0) / 2)
  )

occ_long <- combined_long |>
  filter(year %in% STUDY_YEARS) |>
  inner_join(bin_lookup, by = "month_num") |>
  group_by(AnimalID, year, bin) |>
  summarise(detected = as.integer(any(detected > 0)), .groups = "drop") |>
  left_join(
    occasion_key |> select(year, bin, occasion),
    by = c("year", "bin")
  )

ch_wide <- occ_long |>
  select(AnimalID, occasion, detected) |>
  complete(
    AnimalID,
    occasion = seq_len(n_occasions),
    fill = list(detected = 0L)
  ) |>
  arrange(AnimalID, occasion) |>
  pivot_wider(
    names_from = occasion,
    values_from = detected,
    names_sort = TRUE
  ) |>
  arrange(AnimalID)

# Lions never seen inside the retained bins carry no information
never_seen <- ch_wide$AnimalID[rowSums(as.matrix(ch_wide[, -1])) == 0]

# A lion first seen on the very last occasion also carries no information, and
# breaks the JAGS loop `for (t in (f[i]+1):n.occasions)` -- when f equals
# n.occasions that sequence counts backwards -- so it must go too.
first_occ_all <- apply(as.matrix(ch_wide[, -1]), 1, function(x) {
  if (all(x == 0)) NA_integer_ else min(which(x != 0))
})
last_occ_only <- ch_wide$AnimalID[
  !is.na(first_occ_all) & first_occ_all == n_occasions
]

bin.histories <- ch_wide |>
  filter(!AnimalID %in% c(never_seen, last_occ_only)) |>
  as.data.frame()

cat("\n--- capture history ---\n")
cat(
  "occasions:                    ",
  n_occasions,
  sprintf(" (%d years x %d bins)\n", length(STUDY_YEARS), n_bins)
)
cat("lions with any record:         ", nrow(ch_wide), "\n")
cat("dropped, never seen in bins:   ", length(never_seen), "\n")
cat("dropped, first seen at t=last: ", length(last_occ_only), "\n")
cat("lions retained:                ", nrow(bin.histories), "\n")
cat(
  "total detections:              ",
  sum(as.matrix(bin.histories[, -1])),
  "\n"
)


# ---------------------------------------------------------------------------
# 3. y and f
# ---------------------------------------------------------------------------

animal_ids <- bin.histories$AnimalID # row order for EVERY covariate below
y <- unname(as.matrix(bin.histories[, -1]))
nind <- nrow(y)

get_first <- function(x) min(which(x != 0))
f <- apply(y, 1, get_first)


# ---------------------------------------------------------------------------
# 4. Individual covariates: sex and date of birth
# ---------------------------------------------------------------------------
# Sex and DOB for the 2008-2015 lions come from TM's curated PLOS ID file,
# which takes precedence over the Access database for the same reason the
# capture histories do.

plos_ids_file <- file.path(dirname(folder), "Lion ID 2008-2016 with DOBs.xlsx")

plos_covars <- openxlsx::read.xlsx(plos_ids_file, sheet = "All Lions") |>
  transmute(
    AnimalID = toupper(trimws(MARKID)),
    sex_chr = Sex,
    DOB = openxlsx::convertToDate(DOB),
    err_DOB = suppressWarnings(as.numeric(error.DOB))
  ) |>
  filter(!is.na(AnimalID)) |>
  distinct(AnimalID, .keep_all = TRUE)

access_covars <- lionSight |>
  group_by(AnimalID) |>
  summarise(
    sex_chr = first(na.omit(Sex)),
    DOB = first(na.omit(as.Date(DOB))),
    err_DOB = suppressWarnings(first(na.omit(as.numeric(`error DOB`)))),
    .groups = "drop"
  )

covars <- tibble(AnimalID = animal_ids) |>
  left_join(plos_covars, by = "AnimalID") |>
  left_join(access_covars, by = "AnimalID", suffix = c(".plos", ".access")) |>
  mutate(
    from_plos = AnimalID %in% plos_covars$AnimalID,
    sex_chr = coalesce(sex_chr.plos, sex_chr.access),
    DOB = coalesce(DOB.plos, DOB.access),
    err_DOB = coalesce(err_DOB.plos, err_DOB.access)
  ) |>
  select(AnimalID, from_plos, sex_chr, DOB, err_DOB)

stopifnot(identical(covars$AnimalID, animal_ids)) # row order must not drift

# --- sex: 1 = female, 2 = male ---------------------------------------------
covars <- covars |>
  mutate(
    sex_known = tolower(sex_chr) %in% c("male", "female"),
    sex = case_when(
      tolower(sex_chr) == "female" ~ 1L,
      tolower(sex_chr) == "male" ~ 2L,
      # unknown sex assigned at random (see set.seed above).  The notes at the
      # end of R_JAGS_IPM_CODE.docx say sex should eventually be a random
      # variable estimated inside the model instead.
      TRUE ~ sample(1:2, size = n(), replace = TRUE)
    )
  )

# --- DOB, corrected for the recorded error window --------------------------
# DOB is entered as the earliest plausible date, so the correction shifts it
# toward the middle of the stated error window.  A missing error is treated as
# 0.  (The commented-out version of this had two problems: an NA `error DOB`
# produced an NA corrected DOB, discarding a usable birth date, and the
# `> 270` branch was unreachable because `<= 365` caught it first.)
covars <- covars |>
  mutate(
    err_DOB = replace_na(err_DOB, 0),
    DOB_corr = case_when(
      is.na(DOB) ~ as.Date(NA),
      err_DOB <= 30 ~ DOB,
      err_DOB <= 90 ~ DOB + 45,
      err_DOB <= 180 ~ DOB + 135,
      err_DOB <= 365 ~ DOB + 270,
      TRUE ~ DOB # error > 1 year: no defensible shift
    )
  )

cat("\n--- individual covariates ---\n")
cat(
  "sex known: ",
  sum(covars$sex_known),
  "of",
  nind,
  sprintf("(%d assigned at random)\n", sum(!covars$sex_known))
)
cat("DOB known: ", sum(!is.na(covars$DOB_corr)), "of", nind, "\n")
cat("covariates taken from PLOS ID file:", sum(covars$from_plos), "\n")


# ---------------------------------------------------------------------------
# 5. Age at first detection, and the age-class matrix
# ---------------------------------------------------------------------------

covars$age_at_first <- as.numeric(
  occasion_key$mid_date[f] - covars$DOB_corr
) /
  365.25
covars$age_at_first[covars$age_at_first < 0] <- 0 # bin-midpoint rounding

# Unknown age:
#   known sex   -> median age at first detection among non-cubs
#   unknown sex -> treat as a first-year cub
# This is the rule in the commented-out code.  The notes at the end of
# R_JAGS_IPM_CODE.docx suggest an alternative -- draw a random adult age --
# on the grounds that unaged lions are probably not cubs.  Switch here if you
# prefer that.
median_age_noncub <- median(
  covars$age_at_first[covars$age_at_first >= 2],
  na.rm = TRUE
)

covars <- covars |>
  mutate(
    age_known = !is.na(age_at_first),
    age_at_first = case_when(
      !is.na(age_at_first) ~ age_at_first,
      sex_known ~ median_age_noncub,
      TRUE ~ 1
    )
  )

# Project age forward from first detection across every occasion.  Working
# from the bin midpoints keeps ages correct whatever BIN_MONTHS is set to.
occ_years <- as.numeric(occasion_key$mid_date) / 365.25
years_since_first <- matrix(
  occ_years,
  nrow = nind,
  ncol = n_occasions,
  byrow = TRUE
) -
  matrix(occ_years[f], nrow = nind, ncol = n_occasions)

age_years <- years_since_first +
  matrix(covars$age_at_first, nrow = nind, ncol = n_occasions)
age_years[age_years < 0] <- 0
age <- matrix(get_age_class(age_years), nrow = nind, ncol = n_occasions)

cat("\n--- age ---\n")
cat("age at first detection known:", sum(covars$age_known), "of", nind, "\n")
cat(
  "median age at first detection among non-cubs:",
  round(median_age_noncub, 2),
  "yr\n"
)
cat("age class at first detection:\n")
print(table(
  `age class` = age[cbind(seq_len(nind), f)],
  sex = c("F", "M")[covars$sex]
))
# ---------------------------------------------------------------------------
# 6. Protection stratum, from a spatial overlay of sightings on the park
# ---------------------------------------------------------------------------
# area[i,t] indexes protection, following Creel et al.: 1 = high (inside the
# national park), 2 = low (outside).  It is what identifies beta.prot.
#
# This runs in two passes:
#   pass 1  no assignment file present -> writes SIGHTING_POINTS_FILE, one row
#           per sighting with clean lat/lon, and falls back to all-1s.
#   pass 2  you overlay those points on the park polygon, save the result as
#           PROTECTION_FILE with columns SightingID and protected (1 = inside
#           the park, 0 = outside), and re-run.  area[] is then built from it.

# How a lion's protection stratum is decided for each occasion.  See the
# comment at the assignment below for what the two rules do and why the annual
# one is usually preferable.
PROTECTION_RULE <- "annual"   # "annual" or "occasion"

SIGHTING_POINTS_FILE <- "data/lion_sighting_points.csv"
PROTECTION_FILE <- "data/sighting_protection.csv"

# Plausible bounding box for South Luangwa; anything outside is a data-entry
# error (zeros, truncated values, misplaced decimal points) and is dropped
# rather than guessed at.
LAT_RANGE <- c(-14, -11.5)
LON_RANGE <- c(30.9, 33.5)

sighting_points <- lionSight |>
  transmute(
    SightingID,
    AnimalID,
    SightingDate = as.Date(SightingDate),
    year = year(SightingDate),
    month_num = month(SightingDate),
    lat = suppressWarnings(as.numeric(Lat)),
    lon = suppressWarnings(as.numeric(Lon))
  ) |>
  left_join(bin_lookup, by = "month_num") |>
  left_join(
    occasion_key |> select(year, bin, occasion),
    by = c("year", "bin")
  ) |>
  mutate(
    # A whole-degree coordinate means the minutes were lost on data entry, not
    # that the lion was at that spot.  These pass the bounding box (e.g. the
    # 12 sightings recorded at Lon 31, 11 of them at Lat -13.00000) but would
    # be silently assigned a protection stratum by the overlay, so exclude
    # them explicitly rather than relying on where the box edge happens to sit.
    placeholder = (lat == round(lat)) | (lon == round(lon)),
    coord_ok = !is.na(lat) &
      !is.na(lon) &
      lat > LAT_RANGE[1] &
      lat < LAT_RANGE[2] &
      lon > LON_RANGE[1] &
      lon < LON_RANGE[2] &
      !placeholder
  )

cat("\n--- sighting coordinates ---\n")
cat("animal-sighting rows:        ", nrow(sighting_points), "\n")
cat("with usable coordinates:     ", sum(sighting_points$coord_ok), "\n")
cat(
  "dropped, whole-degree placeholder:",
  n_distinct(sighting_points$SightingID[
    sighting_points$placeholder %in% TRUE
  ]),
  "sightings\n"
)
cat(
  "distinct sightings:          ",
  n_distinct(sighting_points$SightingID),
  "\n"
)
cat(
  "distinct sightings w/ coords:",
  n_distinct(sighting_points$SightingID[sighting_points$coord_ok]),
  "\n"
)

# one row per sighting -- this is what gets overlaid on the park polygon
points_for_overlay <- sighting_points |>
  filter(coord_ok) |>
  distinct(SightingID, .keep_all = TRUE) |>
  select(SightingID, SightingDate, year, month_num, bin, occasion, lat, lon) |>
  arrange(SightingDate)

write.csv(points_for_overlay, SIGHTING_POINTS_FILE, row.names = FALSE)
cat("wrote", nrow(points_for_overlay), "points to", SIGHTING_POINTS_FILE, "\n")

if (file.exists(PROTECTION_FILE)) {
  protection <- read.csv(PROTECTION_FILE, stringsAsFactors = FALSE)
  stopifnot(all(c("SightingID", "protected") %in% names(protection)))

  # 1 = high protection (inside park), 2 = low protection (outside)
  protection <- protection |>
    transmute(
      SightingID,
      area_obs = ifelse(as.integer(protected) == 1L, 1L, 2L)
    ) |>
    filter(!is.na(area_obs))

  # Both rules take the stratum a lion was seen in most often.  Ties are broken
  # by a coin flip rather than by always choosing one stratum: a deterministic
  # tie-break would push every evenly-split lion the same way and bias the
  # composition.  This makes area[] stochastic, so it depends on the seed set at
  # the top of this script -- re-running with a different seed will move a
  # handful of classifications.
  #
  # The two rules differ in the window the majority is taken over.
  #
  # "occasion": majority within each 2-month bin.  Keeps within-year movement,
  #             but a bin with a single sighting is classified on that one
  #             sighting, and bins with none are filled by carry-forward.
  #
  # "annual":   majority over ALL of a lion's sightings in the year, applied to
  #             every occasion in that year.  Each classification rests on more
  #             sightings (median 4 per lion-year against 2 per lion-occasion),
  #             so it is far less sensitive to a single opportunistic sighting
  #             outside a lion's usual range.  The cost is that within-year
  #             movement is averaged away.
  #
  # The annual rule is a better match for how the protection effect is
  # interpreted -- as a property of where a lion lived that year, not where it
  # happened to be on one afternoon.
  protection_by_sighting <- sighting_points |>
    inner_join(protection, by = "SightingID") |>
    filter(!is.na(occasion), AnimalID %in% animal_ids)

  area_obs_mat <- matrix(NA_integer_, nrow = nind, ncol = n_occasions)

  if (PROTECTION_RULE == "occasion") {
    animal_occ_area <- protection_by_sighting |>
      count(AnimalID, occasion, area_obs) |>
      group_by(AnimalID, occasion) |>
      slice_max(order_by = n, n = 1, with_ties = TRUE) |>
      summarise(
        # area_obs[sample.int(...)], never sample(area_obs, 1): with a single
        # tied row holding the value 2, sample(2, 1) draws from 1:2 instead of
        # returning 2
        area_obs = area_obs[sample.int(length(area_obs), 1)],
        .groups = "drop"
      )

    area_obs_mat[cbind(
      match(animal_occ_area$AnimalID, animal_ids),
      animal_occ_area$occasion
    )] <- animal_occ_area$area_obs
  } else if (PROTECTION_RULE == "annual") {
    animal_year_area <- protection_by_sighting |>
      count(AnimalID, year, area_obs) |>
      group_by(AnimalID, year) |>
      slice_max(order_by = n, n = 1, with_ties = TRUE) |>
      summarise(
        area_obs = area_obs[sample.int(length(area_obs), 1)],
        .groups = "drop"
      )

    # one classification per lion-year, written to every occasion in that year
    year_occ <- occasion_key |> select(year, occasion)
    expanded <- animal_year_area |>
      inner_join(year_occ, by = "year", relationship = "many-to-many")

    area_obs_mat[cbind(
      match(expanded$AnimalID, animal_ids),
      expanded$occasion
    )] <- expanded$area_obs
  } else {
    stop("PROTECTION_RULE must be 'occasion' or 'annual'")
  }

  # Carry the last known stratum forward across occasions when the lion was
  # not seen, then back-fill before its first assignment.
  area <- t(apply(area_obs_mat, 1, function(r) {
    if (all(is.na(r))) {
      return(rep(1L, length(r))) # never located: default to high protection
    }
    known <- which(!is.na(r))
    r[seq_len(min(known))] <- r[min(known)] # back-fill
    for (k in seq_along(r)[-1]) {
      if (is.na(r[k])) r[k] <- r[k - 1] # carry forward
    }
    as.integer(r)
  }))

  cat("\n--- protection stratum (from", PROTECTION_FILE, ") ---\n")
  cat("lion-occasions assigned directly:", sum(!is.na(area_obs_mat)), "\n")
  cat(
    "lions never located:             ",
    sum(apply(area_obs_mat, 1, function(r) all(is.na(r)))),
    "\n"
  )
  cat("lion-occasions high protection:  ", sum(area == 1), "\n")
  cat("lion-occasions low protection:   ", sum(area == 2), "\n")
  cat(
    "lions that switched stratum:     ",
    sum(apply(area, 1, function(r) length(unique(r)) > 1)),
    "\n"
  )
} else {
  area <- matrix(1L, nrow = nind, ncol = n_occasions)
  cat(
    "\n*** ",
    PROTECTION_FILE,
    " not found: area[] is all 1s (placeholder).\n",
    "*** beta.prot is NOT identifiable until you supply it. Overlay\n",
    "*** ",
    SIGHTING_POINTS_FILE,
    " on the park polygon and save the result\n",
    "*** with columns SightingID, protected (1 = inside, 0 = outside).\n",
    sep = ""
  )
}


# ---------------------------------------------------------------------------
# 7. Bundle for the CJS / IPM survival component
# ---------------------------------------------------------------------------

# Exponents that annualise the per-interval survival estimated by the CJS.
# With n_bins occasions per year, annual survival is phi^n_bins.
#
# Cubs need a larger exponent.  The cub age class is 24 months long, but a cub
# is not detectable for roughly the first CUB_UNSEEN_MONTHS of life (denning),
# so phi for that class is estimated over a shorter window and has to be
# stretched to cover the full 24 months.  Creel et al. used 4.41 unseen months
# with 3 bins, giving 3 * 24/(24-4.41) = 3.672; the same rule with 4 bins
# gives 4.900.  These are passed as data so they follow BIN_MONTHS rather than
# sitting as magic numbers in the JAGS file.
CUB_UNSEEN_MONTHS <- 4.41

surv_exp <- n_bins
surv_exp_cub <- n_bins * 24 / (24 - CUB_UNSEEN_MONTHS)

cat("\n--- survival annualisation ---\n")
cat("bins per year:            ", n_bins, "\n")
cat("exponent, ages 2-4:       ", surv_exp, "\n")
cat("exponent, cubs:           ", round(surv_exp_cub, 3), "\n")

jags.data <- list(
  y = y,
  f = f,
  nind = nind,
  n.occasions = n_occasions,
  nageclass = 4,
  age = age,
  sex = covars$sex,
  area = area,
  n.bins = n_bins,
  surv.exp = surv_exp,
  surv.exp.cub = surv_exp_cub
)

# Row order and dimensions must line up with y, or the model silently fits one
# lion's covariates to another lion's history.
stopifnot(
  nrow(jags.data$age) == nind,
  ncol(jags.data$age) == n_occasions,
  length(jags.data$sex) == nind,
  all(dim(jags.data$area) == c(nind, n_occasions)),
  all(jags.data$f < n_occasions),
  all(!is.na(jags.data$age)),
  all(jags.data$sex %in% 1:2),
  all(jags.data$area %in% 1:2)
)
str(jags.data)


# ---------------------------------------------------------------------------
# 8. Closed-capture histories, one matrix per year (chc)
# ---------------------------------------------------------------------------
# Fit these outside the IPM, one model per year, using the same bins and the
# same logit-normal individual random effect on p as the CJS model.  The
# posterior medians of N become Pop[t] in the state-space component.

chc <- lapply(STUDY_YEARS, function(yr) {
  cols <- occasion_key$occasion[occasion_key$year == yr]
  m <- y[, cols, drop = FALSE]
  m[rowSums(m) > 0, , drop = FALSE] # closed capture: seen at least once
})
names(chc) <- STUDY_YEARS

# Monthly version of the same within-year histories.  The CJS needs 2-month
# bins so that survival intervals are equal, but the closed-capture model has
# no such constraint -- it only assumes closure within the year.  Four
# occasions leaves the Mh model badly under-identified (the individual
# heterogeneity SD and N trade off along a ridge, so N is barely pinned down);
# splitting the same months into single-month occasions doubles the occasions
# at no cost in assumptions.  Same months as BIN_MONTHS, one column each.
chc_months <- sort(unlist(BIN_MONTHS, use.names = FALSE))

chc_monthly <- lapply(STUDY_YEARS, function(yr) {
  m <- combined_long |>
    filter(year == yr, month_num %in% chc_months) |>
    group_by(AnimalID, month_num) |>
    summarise(detected = as.integer(any(detected > 0)), .groups = "drop") |>
    filter(AnimalID %in% animal_ids) |>
    complete(
      AnimalID = animal_ids,
      month_num = chc_months,
      fill = list(detected = 0L)
    ) |>
    arrange(AnimalID, month_num) |>
    pivot_wider(
      names_from = month_num,
      values_from = detected,
      names_sort = TRUE
    )
  out <- as.matrix(m[, -1])
  rownames(out) <- m$AnimalID
  out[rowSums(out) > 0, , drop = FALSE]
})
names(chc_monthly) <- STUDY_YEARS

cat("\n--- monthly closed-capture histories ---\n")
cat("occasions per year:", length(chc_months), "\n")
print(
  data.frame(
    year = STUDY_YEARS,
    n_known = vapply(chc_monthly, nrow, 0L),
    n_detections = vapply(chc_monthly, sum, 0L)
  ),
  row.names = FALSE
)

cat("\n--- closed-capture histories by year ---\n")
print(
  data.frame(
    year = STUDY_YEARS,
    n_known = vapply(chc, nrow, 0L),
    n_detections = vapply(chc, sum, 0L)
  ),
  row.names = FALSE
)


# ---------------------------------------------------------------------------
# 9. Fecundity: first-year cubs per adult female per year
# ---------------------------------------------------------------------------
# Feeds the zero-inflated Poisson component:  C[i] ~ dpois(mu.f[i]),
# log(lambda[i]) <- alpha + gamma[region[i]].  One row per adult-female-year.
#
# Litter records come from 'first luangwa paper/litters.xlsx' (2006-2015).
# What is observed reliably is the NUMBER and SIZE of litters in each
# pride-year; what is mostly missing is which female produced each litter.
# Mother is named for only 33 of 88 litters in 2008-2015, covering 76 of 223
# cubs, so attributing cubs to named mothers alone and calling every other
# adult female a zero puts the mean at 0.31 cubs/female/yr against 1.18 from
# the pride-level totals -- a four-fold downward bias, with 94% zeros driving
# the zero-inflation parameter psi.
#
# Instead the litters are allocated within their pride: females with a named
# litter get it, and the remaining litters in that pride-year are assigned at
# random among the remaining adult females, one litter each (lions raise one
# litter at a time).  Observed litter counts and sizes are preserved exactly;
# only the identity of the mother is imputed.  Set FEC_IMPUTE_LITTERS to FALSE
# to use named mothers only.

LITTERS_FILE <- file.path(dirname(folder), "litters.xlsx")
FEC_YEARS <- 2008:2015 # extent of usable litter records
FEC_MIN_AGE <- 3 # breeding females are Leslie age classes 4-7
FEC_IMPUTE_LITTERS <- TRUE

litters_raw <- openxlsx::read.xlsx(
  LITTERS_FILE,
  sheet = "litters(<2yrs at first detect)"
)

clean_mark_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[!str_detect(x, "^ELI-[[:space:]]*[[:digit:]]+$")] <- NA_character_
  str_replace(x, "ELI-[[:space:]]*", "ELI-")
}

# Pride names differ between sources ("bigpride"/"Big", "mwamba2"/"Mwamba II").
# Strip punctuation and case, and fold trailing Roman numerals to digits.
normalise_pride <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- str_replace_all(x, "[^a-z0-9]", "")
  x <- str_replace(x, "iii$", "3")
  x <- str_replace(x, "ii$", "2")
  x <- str_replace(x, "i$", "1")
  str_replace(x, "pride$", "")
}

litters <- litters_raw |>
  mutate(
    litter = row_number(),
    n_cubs = coalesce(size, males + females + unknown),
    pride_key = normalise_pride(pride),
    mother_direct = clean_mark_id(mother)
  )

# --- pool the cub -> mother links from all three sources -------------------
link_litters <- litters |>
  select(litter, mother_direct, ID1:ID5) |>
  pivot_longer(ID1:ID5, values_to = "cub", names_to = NULL) |>
  filter(!is.na(cub)) |>
  transmute(cub = toupper(trimws(cub)), mom = mother_direct, source = "litters")

link_plos <- openxlsx::read.xlsx(plos_ids_file, sheet = "All Lions") |>
  transmute(
    cub = toupper(trimws(MARKID)),
    mom = clean_mark_id(Suspected.Mother),
    source = "plos_id"
  )

link_access <- lionSight |>
  group_by(AnimalID) |>
  summarise(m = first(na.omit(`Suspected Mother`)), .groups = "drop") |>
  transmute(cub = AnimalID, mom = clean_mark_id(m), source = "access")

cub_mother <- bind_rows(link_litters, link_plos, link_access) |>
  filter(!is.na(mom), !is.na(cub))

cat("\n--- cub-to-mother links ---\n")
print(table(source = cub_mother$source))
cat(
  "cubs with conflicting mothers:",
  nrow(cub_mother |> distinct(cub, mom) |> count(cub) |> filter(n > 1)),
  "\n"
)

cub_mother <- cub_mother |> distinct(cub, .keep_all = TRUE)

litter_mother <- litters |>
  select(litter, ID1:ID5) |>
  pivot_longer(ID1:ID5, values_to = "cub", names_to = NULL) |>
  filter(!is.na(cub)) |>
  mutate(cub = toupper(trimws(cub))) |>
  left_join(cub_mother |> select(cub, mom), by = "cub") |>
  group_by(litter) |>
  summarise(
    mom_from_cubs = {
      v <- unique(na.omit(mom))
      if (length(v) == 1) v else NA_character_
    },
    .groups = "drop"
  )

litters <- litters |>
  left_join(litter_mother, by = "litter") |>
  mutate(mother_final = coalesce(mother_direct, mom_from_cubs))

litters_used <- litters |>
  filter(yearofbirth %in% FEC_YEARS, !is.na(pride_key), !is.na(n_cubs))

cat("\n--- litters ---\n")
cat("litters in file:           ", nrow(litters), "\n")
cat(
  "usable,",
  min(FEC_YEARS),
  "-",
  max(FEC_YEARS),
  ":       ",
  nrow(litters_used),
  "\n"
)
cat("with an identified mother: ", sum(!is.na(litters_used$mother_final)), "\n")
cat(
  "cubs total / attributed:   ",
  sum(litters_used$n_cubs),
  "/",
  sum(litters_used$n_cubs[!is.na(litters_used$mother_final)]),
  "\n"
)

# --- pride membership per lion-year ----------------------------------------
pride_by_year <- lionSight |>
  filter(!is.na(AnimalGroup), !is.na(SightingDate)) |>
  mutate(year = year(SightingDate)) |>
  group_by(AnimalID, year) |>
  summarise(
    pride_key = names(sort(table(AnimalGroup), decreasing = TRUE))[1],
    .groups = "drop"
  ) |>
  mutate(pride_key = normalise_pride(pride_key))

litter_prides <- sort(unique(litters_used$pride_key))
cat(
  "\nprides in litters.xlsx matched to AnimalGroup:",
  sum(litter_prides %in% pride_by_year$pride_key),
  "of",
  length(litter_prides),
  "\n"
)

# --- the adult-female-year frame -------------------------------------------
# A female enters for year t if she was detected in t, was at least
# FEC_MIN_AGE at the midpoint of t, and belonged to a pride whose litters were
# recorded.  Females in prides absent from litters.xlsx are excluded: we have
# no reproductive information for them, so counting them as zeros would be the
# same false-zero problem in another form.
occ_first_of_year <- occasion_key |>
  group_by(year) |>
  summarise(first_occ = min(occasion), .groups = "drop")

seen_in_year <- vapply(
  FEC_YEARS,
  function(yr) {
    cols <- occasion_key$occasion[occasion_key$year == yr]
    rowSums(y[, cols, drop = FALSE]) > 0
  },
  logical(nind)
)
colnames(seen_in_year) <- FEC_YEARS

female_years <- expand_grid(
  row = which(covars$sex == 1L),
  year = FEC_YEARS
) |>
  left_join(occ_first_of_year, by = "year") |>
  mutate(
    AnimalID = animal_ids[row],
    age_yr = age_years[cbind(row, first_occ)],
    seen = seen_in_year[cbind(row, match(year, FEC_YEARS))]
  ) |>
  filter(seen, age_yr >= FEC_MIN_AGE) |>
  left_join(pride_by_year, by = c("AnimalID", "year")) |>
  filter(!is.na(pride_key), pride_key %in% litter_prides)

cat("\n--- adult-female-year frame ---\n")
cat("female-years:      ", nrow(female_years), "\n")
cat("females:           ", n_distinct(female_years$AnimalID), "\n")
cat("prides:            ", n_distinct(female_years$pride_key), "\n")

# --- allocate litters to females within pride-years ------------------------
assign_litters <- function(fy, lt, impute) {
  fy$C <- 0L
  # named mothers first
  named <- lt |> filter(!is.na(mother_final))
  hit <- match(named$mother_final, fy$AnimalID)
  for (k in which(!is.na(hit))) {
    fy$C[hit[k]] <- fy$C[hit[k]] + named$n_cubs[k]
  }
  claimed <- unique(na.omit(hit))

  if (!impute) {
    return(list(
      fy = fy,
      unplaced = sum(named$n_cubs[is.na(hit)]) +
        sum(lt$n_cubs[is.na(lt$mother_final)])
    ))
  }

  # litters with no usable mother, plus named mothers not in this frame
  spare <- c(
    lt$n_cubs[is.na(lt$mother_final)],
    named$n_cubs[is.na(hit)]
  )
  if (length(spare) == 0) {
    return(list(fy = fy, unplaced = 0))
  }

  free <- setdiff(seq_len(nrow(fy)), claimed)
  n_take <- min(length(spare), length(free))
  if (n_take > 0) {
    take <- if (length(free) == 1) free else sample(free, n_take)
    fy$C[take] <- fy$C[take] + spare[seq_len(n_take)]
  }
  # more litters than adult females in the frame: pile the rest onto random
  # females rather than discarding observed cubs
  leftover <- spare[-seq_len(n_take)]
  if (length(leftover) > 0) {
    extra <- if (nrow(fy) == 1) {
      rep(1L, length(leftover))
    } else {
      sample(seq_len(nrow(fy)), length(leftover), replace = TRUE)
    }
    for (k in seq_along(leftover)) {
      fy$C[extra[k]] <- fy$C[extra[k]] + leftover[k]
    }
  }
  list(fy = fy, unplaced = 0)
}

fec_parts <- female_years |>
  group_by(pride_key, year) |>
  group_split()

unplaced_total <- 0
fec_data <- lapply(fec_parts, function(fy) {
  lt <- litters_used |>
    filter(pride_key == fy$pride_key[1], yearofbirth == fy$year[1])
  out <- assign_litters(as.data.frame(fy), lt, FEC_IMPUTE_LITTERS)
  unplaced_total <<- unplaced_total + out$unplaced
  out$fy
}) |>
  bind_rows()

# cubs observed vs cubs represented in C[] -- these must agree when imputing
cubs_in_frame <- litters_used |>
  semi_join(female_years, by = c("pride_key", "yearofbirth" = "year")) |>
  summarise(s = sum(n_cubs)) |>
  pull(s)

if (!FEC_IMPUTE_LITTERS) {
  cat(
    "cubs dropped for want of a named mother:",
    unplaced_total,
    "(FEC_IMPUTE_LITTERS is FALSE)\n"
  )
}
cat("cubs in matched pride-years:", cubs_in_frame, "\n")
cat("cubs represented in C[]:    ", sum(fec_data$C), "\n")
cat(
  "cubs in litters with no matching pride-year in the frame:",
  sum(litters_used$n_cubs) - cubs_in_frame,
  "\n"
)

# --- region[i] and ordering ------------------------------------------------
fec_data <- fec_data |>
  mutate(region = area[cbind(row, first_occ)]) |>
  select(row, AnimalID, year, pride_key, age_yr, C, region)

# The JAGS code partitions mu.f by region with hard-coded index ranges
# (fec.prot1 <- mean(mu.f[1:74])), so rows must be sorted by region and the
# boundaries read off rather than assumed.
fec_data <- fec_data |> arrange(region, year, pride_key, AnimalID)

fec_bounds <- fec_data |>
  mutate(i = row_number()) |>
  group_by(region) |>
  summarise(from = min(i), to = max(i), n = n(), .groups = "drop")

cat("\n--- fecundity data ---\n")
cat("adult-female-years:", nrow(fec_data), "\n")
cat("distribution of C (first-year cubs per female-year):\n")
print(table(C = fec_data$C))
cat(
  "mean C:",
  round(mean(fec_data$C), 3),
  "  proportion zero:",
  round(mean(fec_data$C == 0), 3),
  "\n"
)
cat("by year:\n")
print(
  fec_data |>
    group_by(year) |>
    summarise(
      female_years = n(),
      cubs = sum(C),
      mean_C = round(mean(C), 2),
      .groups = "drop"
    ) |>
    as.data.frame(),
  row.names = FALSE
)
cat(
  "index ranges for the JAGS model (fec.prot<region> <- mean(mu.f[from:to])):\n"
)
print(as.data.frame(fec_bounds), row.names = FALSE)

if (FEC_IMPUTE_LITTERS) {
  cat(
    "\nNOTE: mother identity was imputed within pride-years for ",
    sum(litters_used$n_cubs) -
      sum(litters_used$n_cubs[!is.na(litters_used$mother_final)]),
    " of ",
    sum(litters_used$n_cubs),
    " cubs.\n",
    "Litter counts and sizes are observed; only which female produced each\n",
    "litter is imputed.  Re-run with a different set.seed() to check that the\n",
    "fecundity posterior is not sensitive to the draw.\n",
    sep = ""
  )
}

# Pride-level totals, as a cross-check and as the basis for the Poisson GLM of
# cubs per pride per year that the paper reports alongside the ZIP model.
fec_pride <- litters_used |>
  group_by(pride_key, year = yearofbirth) |>
  summarise(n_litters = n(), cubs = sum(n_cubs), .groups = "drop")

jags.data.fec <- list(
  C = fec_data$C,
  region = fec_data$region,
  n = nrow(fec_data)
)


# ---------------------------------------------------------------------------
# 10. Save
# ---------------------------------------------------------------------------

lion_ipm_data <- list(
  jags.data = jags.data,
  jags.data.fec = jags.data.fec,
  chc = chc,
  chc_monthly = chc_monthly,
  animal_ids = animal_ids,
  covars = covars,
  occasion_key = occasion_key,
  fec_data = fec_data,
  fec_bounds = fec_bounds,
  fec_pride = fec_pride,
  litters = litters,
  sighting_points = points_for_overlay,
  settings = list(
    STUDY_YEARS = STUDY_YEARS,
    BIN_MONTHS = BIN_MONTHS,
    PLOS_YEARS = PLOS_YEARS,
    FEC_YEARS = FEC_YEARS,
    FEC_MIN_AGE = FEC_MIN_AGE
  )
)

saveRDS(jags.data, "data/jags_data_cjs_ipm.RDS")
saveRDS(lion_ipm_data, "data/lion_ipm_data.RDS")


# ---------------------------------------------------------------------------
# 11. Remaining gaps
# ---------------------------------------------------------------------------
cat(
  "
=============================================================================
REMAINING GAPS
=============================================================================
1. Pop[t] -- run the closed-capture model on chc (one model per year), then
   pass the posterior medians of N into the state-space component, along with
   Area for the density calculation.

2. Fecundity years.  litters.xlsx stops at 2015, so C[] covers 2008-2015
   while the capture histories run to 2023.  Either restrict the IPM to
   2008-2015, or extend the litter records.  What the Access database offers:
     - Group composition (Ad Fem, Male Juv, Fem Juv, Juv unk in
       tblAnimalSightings) is ~100% complete for all 16 years and gives cubs
       and adult females per pride per sighting.  This supports a pride-level
       fecundity series to 2023, but not individual-level C[i].
     - 'Suspected Mother' in tblAnimals adds mother links, but covers 171 of
       604 lions and thins out badly after 2017.
     - tblPregnancy holds lactation/pregnancy observations, which could
       identify breeding females but not litter sizes.
     - tblDen and tblDenVisit are wild dog only.

3. Cub attribution.  See the warning above if it printed.  Options: accept
   the downward bias, restrict the female frame to prides and years where
   every litter has a named mother, or move the fecundity component to the
   pride level (fec_pride) with an offset for the number of adult females.
=============================================================================
"
)
