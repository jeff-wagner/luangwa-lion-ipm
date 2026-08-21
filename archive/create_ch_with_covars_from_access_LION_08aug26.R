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
  "openxlsx"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

lapply(required_packages, library, character.only = TRUE)

#read the relevant tables from Access
#alternatively could use the dbi and dbplyr packages to connect to the Access
#database directly, but if the files fit in memory that would be slower
#and have no advantages, see: https://db.rstudio.com/dplyr/

# Open DB connection
dbCon <- odbcConnectAccess2007("Luangwa Database - 2024Nov22 - MDF.accdb")
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
folder <- "H:/Biometric/Github/luangwa-lion-ipm/first luangwa paper/2008-2015excel"
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


# Bundle data for JAGS
y <- as.matrix(bin.histories[, -1]) # Remove AnimalID column
y <- unname(y)
n_ind <- nrow(y)
n_occasions <- ncol(y)

# Get first capture occasion for each individual
get_first <- function(x) min(which(x != 0))
first <- apply(y, 1, get_first)

# Bundle data for JAGS
jags.data <- list(
  y = y,
  nind = n_ind,
  nyears = n_occasions,
  f = first
)
str(jags.data)

# # Age and sex data --------------------------------------------------------
# # Get sex for each individual
# sex <- lionSight %>%
#   select(AnimalID, Sex) %>%
#   arrange(AnimalID) %>%
#   distinct()
# #0 is female, 1 is male, NA is unkown
# sexcode <- 1 * (tolower(sex$Sex) == 'male')
# sexcode[!(tolower(sex$Sex) %in% c('male', 'female'))] <- NA
# # Get age at first detection for CMR models by comparison of DOB and date of first detection

# #age (in months) at first detection
# age_at_first_months <- lionSight %>%
#   mutate(SightingDate = as.Date(SightingDate), DOB = as.Date(DOB)) %>%
#   group_by(AnimalID, DOB) %>%
#   summarise(first_detection = min(SightingDate), .groups = "drop") %>%
#   mutate(
#     age_months_at_first = round(
#       as.numeric(first_detection - DOB) / (365.25 / 12),
#       1
#     )
#   )

# age_at_first_months$AnimalID[which(
#   age_at_first_months$age_months_at_first <= 0
# )]
# age_at_first_months$age_months_at_first[which(
#   age_at_first_months$age_months_at_first <= 0
# )]

# #year of first detection
# age_at_first <- data.frame(AnimalID = bin.histories$AnimalID, sex = sex$Sex)
# age_at_first$year_first_cap <- apply(bin.histories[, -1], 1, get_first) + 2007

# # Add DOB and DOB error
# age_at_first <- age_at_first %>%
#   left_join(
#     lionSight %>% select(AnimalID, DOB, `error DOB`),
#     by = "AnimalID"
#   ) %>%
#   mutate(
#     DOB = as_date(DOB),
#     midp_first_cap = as_date(paste0(year_first_cap, "-07-01")),
#     year = year(DOB),
#     age_at_first = midp_first_cap - DOB
#   ) %>%
#   distinct()

# # Decide how to handle DOB error - assign corrected DOBs
# age_at_first <- age_at_first %>%
#   mutate(
#     DOBcorrected = case_when(
#       `error DOB` <= 30 ~ DOB,
#       `error DOB` > 30 & `error DOB` <= 90 ~ DOB + 45,
#       `error DOB` > 90 & `error DOB` <= 180 ~ DOB + 135,
#       `error DOB` > 180 & `error DOB` <= 365 ~ DOB + 270,
#       `error DOB` > 270 ~ DOB
#     )
#   )

# # Add ages in days for all years
# # Function to assign age class based on age in years
# get_age_class <- function(age_years) {
#   case_when(
#     age_years < 2 ~ 1, # cub (0-1.99)
#     between(age_years, 2, 3.99) ~ 2, # subadult (2-3.99)
#     between(age_years, 4, 5.99) ~ 3, # young adult (4-5.99)
#     age_years >= 6 ~ 4
#   ) # old adult (6+)
# }

# # Convert age at first capture from days to years
# age_at_first$age_at_first_years <- as.numeric(
#   age_at_first$age_at_first,
#   units = "days"
# ) /
#   365.25

# # For unknown age, known sex: assign median age of non-cubs
# # For unknown age, unknown sex: make 1st year cubs
# # FOr unknown sex, assign random
# median.age.adult <- age_at_first %>%
#   filter(age_at_first_years >= 2) %>%
#   pull(age_at_first_years) %>%
#   median(na.rm = TRUE)
# age_at_first <- age_at_first %>%
#   mutate(
#     age_at_first_years = case_when(
#       is.na(age_at_first_years) & sex != "Unknown" ~ median.age.adult,
#       is.na(age_at_first_years) & sex == "Unknown" ~ 1,
#       .default = age_at_first_years
#     ),
#     sex = case_when(
#       sex == "Unknown" ~ sample(
#         c("Male", "Female"),
#         1,
#         replace = TRUE,
#         prob = c(0.5, 0.5)
#       ),
#       .default = sex
#     )
#   )

# # Get the range of years in your study
# all_years <- sort(unique(lionSight$year)) # Adjust to match your nyears=16
# nyears <- length(all_years)
# nind <- nrow(age_at_first)

# # Initialize age matrix
# age <- matrix(NA, nrow = nind, ncol = nyears)

# # Fill in age class for each individual at each time point
# for (i in 1:nind) {
#   first_year <- age_at_first$year_first_cap[i]
#   first_year_index <- which(all_years == first_year)
#   first_year_age <- age_at_first$age_at_first_years[i]

#   for (t in first_year_index:nyears) {
#     # Calculate age at time t
#     years_since_first <- t - first_year_index
#     current_age <- first_year_age + years_since_first

#     # Assign age class
#     age[i, t] <- get_age_class(current_age)
#   }
# }

# age.mat <- createAge(
#   f = first,
#   age = get_age_class(age_at_first$age_at_first_years),
#   nyears = nyears,
#   mAge = 4
# )

# #Month of first detection, assigned to the middle of each bin/window in the ch's (e.g. May would be middle of April, May, June bin)
# #NOTE need a number of if()s equal to number of bins/yr, and
# #number of logicals within each if equal to number of years, and
# #will need to tweak if ch is not may/jun jul/aug sept/oct within each year
# # for (i in 1:length(capt.hist$ch))

# ## 3 3-mo bins 8 years for a total of 24 bins. --BG

# for (i in 1:length(capt.hist$ch)) {
#   if (
#     regexpr("1", capt.hist$ch[i]) == 1 |
#       regexpr("1", capt.hist$ch[i]) == 4 |
#       regexpr("1", capt.hist$ch[i]) == 7 |
#       regexpr("1", capt.hist$ch[i]) == 10 |
#       regexpr("1", capt.hist$ch[i]) == 13 |
#       regexpr("1", capt.hist$ch[i]) == 16 |
#       regexpr("1", capt.hist$ch[i]) == 19 |
#       regexpr("1", capt.hist$ch[i]) == 22
#   ) {
#     capt.hist$monthatfirst[i] <- 5
#   }

#   if (
#     regexpr("1", capt.hist$ch[i]) == 2 |
#       regexpr("1", capt.hist$ch[i]) == 5 |
#       regexpr("1", capt.hist$ch[i]) == 8 |
#       regexpr("1", capt.hist$ch[i]) == 11 |
#       regexpr("1", capt.hist$ch[i]) == 14 |
#       regexpr("1", capt.hist$ch[i]) == 17 |
#       regexpr("1", capt.hist$ch[i]) == 20 |
#       regexpr("1", capt.hist$ch[i]) == 23
#   ) {
#     capt.hist$monthatfirst[i] <- 8
#   }

#   if (
#     regexpr("1", capt.hist$ch[i]) == 3 |
#       regexpr("1", capt.hist$ch[i]) == 6 |
#       regexpr("1", capt.hist$ch[i]) == 9 |
#       regexpr("1", capt.hist$ch[i]) == 12 |
#       regexpr("1", capt.hist$ch[i]) == 15 |
#       regexpr("1", capt.hist$ch[i]) == 18 |
#       regexpr("1", capt.hist$ch[i]) == 21 |
#       regexpr("1", capt.hist$ch[i]) == 24
#   ) {
#     capt.hist$monthatfirst[i] <- 11
#   }
# }

# #^ 5, 8, and 11 correspond to the middle of each of the 3 month bins (e.g. in this example it is may, august, and november)

# ## 2 6-mo bins --BG
# # for (i in 1:length(capt.hist$ch)) {
# #
# #   if (regexpr("1",capt.hist$ch[i]) == 1 |regexpr("1",capt.hist$ch[i]) == 3| regexpr("1",capt.hist$ch[i]) == 5|regexpr("1",capt.hist$ch[i]) == 7 |regexpr("1",capt.hist$ch[i]) == 9|
# #       regexpr("1",capt.hist$ch[i]) == 11 |regexpr("1",capt.hist$ch[i]) == 13) capt.hist$monthatfirst[i] <- 3
# #
# #   if (regexpr("1",capt.hist$ch[i]) == 2 | regexpr("1",capt.hist$ch[i]) == 4| regexpr("1",capt.hist$ch[i]) == 6 | regexpr("1",capt.hist$ch[i]) == 8 | regexpr("1",capt.hist$ch[i]) == 10|
# #       regexpr("1",capt.hist$ch[i]) == 4| regexpr("1",capt.hist$ch[i]) == 12| regexpr("1",capt.hist$ch[i]) == 1) capt.hist$monthatfirst[i] <- 9
# # }

# #RENAME DOB COLUMN
# capt.hist$DOB <- capt.hist$'The date the animal is born in (DD/MM/YYYY), if only a year is known note: yyyy/01/01'

# capt.hist$monthatfirst
# #convert DOB in master ID file to an R date (NOTE POSIXlt, not DATE)
# capt.hist$DOB <- strptime(capt.hist$DOB, "%m/%d/%Y")

# # estimate age at first detection using the month and year of first detection and DOB.

# # Assume day of month at first detection was the 15th to avoid negative ages.
# capt.hist$ageatfirst <- as.numeric(
#   difftime(
#     strptime(
#       paste(capt.hist$yearatfirst, capt.hist$monthatfirst, "28", sep = "-"),
#       "%Y-%m-%d"
#     ),
#     capt.hist$DOB,
#     units = "days"
#   ) /
#     365
# )
# #a few animals have 1990 DOB presumably to indicate they were in oldest age class when detected, e.g. 10 years old
# capt.hist$ageatfirst[which(capt.hist$ageatfirst > 10)] <- 10
# capt.hist
