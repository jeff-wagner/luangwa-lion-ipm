# KNP WD CMR analysis
# processing tables from Access to create input files specifically
# for annual Bayesian closed capture models with adata augmentation


rm(list = ls()) 

library(dplyr)
library(lubridate)
#library(purrr)
library(reshape)
#library(readr)
#library(tidyr)
#library(tibble)
library(RODBC)
library(stringr)
library(IPMbook)

#read the relevant tables from Access
#alternatively could use the dbi and dbplyr packages to connect to the Access
#database directly, but if the files fit in memory that would be slower
#and have no advantages, see: https://db.rstudio.com/dplyr/

# Open DB connection
dbCon <- odbcConnectAccess2007("./Luangwa Database - 2024Nov22 - MDF.accdb")
sqlTables(dbCon)

lionSight <- sqlQuery(dbCon, "
                      SELECT
                        AnimSight.*,
                        Sight.*,
                        Anim.*
                      FROM (tblAnimalsatSighting AS AnimSight
                        LEFT JOIN tblSightings AS Sight ON AnimSight.SightingID = Sight.SightingID)
                        LEFT JOIN tblAnimals AS Anim ON AnimSight.AnimalID = Anim.AnimalID
                      WHERE Sight.LionSighting = TRUE
                      ")

# Close DB connection
odbcClose(dbCon)

#note that this does not eliminate some rows for other species that were present at dog sightings
#but will eliminate those cases at final join of capture histories with covariates

#convert date from character to date format (m/d/yyyy)
#NOTE THAT THIS STEP HAS NOT BEEN CONSISTENT BETWEEN IMPORTS FROM ACCESS DEPENDING ON THE COMPUTER USED
#tweak the format argument as needed

lionSight$SightingDate <- as.Date.character(x=lionSight$SightingDate, format = "%Y-%m-%d")

# QC for incorrect dates, AnimalIDs and species
lionSight %>% 
  filter(lubridate::year(SightingDate) < 1900) # Looking in database, seems like this entry was a typo. Should be 2017-06-17
  
lionSight[lionSight$SightingID %in% "E-18435", "SightingDate"] <- "2017-06-17"

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
table(unique(lionSight$AnimalID) %in% animal.ids)

# add a detected variable, all ones at this point for the sightings themselves
lionSight$detect <- rep(1, length(lionSight$SightingID))

#reshape to pivot the data
# temp<-melt(lionsight,id.var=c("AnimalID","SightingDate"),measure.var="detect")
# histories=cast(temp,AnimalID ~SightingDate)

#structurally, looks good (animals as rows, sightings as columns, 0/1 for not/detected in the body
#but one column for every date with a lion sighting
#so first need to collapse columns into regular set of time windows
#in correct part of each year

#create date bins first, and then melt by AnimalID and bin
#first need to identify the best set of time bins
hist(lionSight$SightingDate, breaks = 20)

#explore distribution of dates to decide on bins
hist(lionSight$SightingDate, breaks = 100, freq = TRUE,
     xlab = 'Month', ylab = 'Sightings', main= 'KNP AWD')

#how many sightings per month?
month<- month(lionSight$SightingDate)
hist(month, breaks = 11)  
#months 5/6, 7/8, 9/10 as 3 bins/year include almost all of the data, pretty even among bins,
#for 2013,2014,2014,2016

# 3 bins of 3-4 months per year was selected by looking at alternative binnings (code deleted)
# in each year, bin 1 = months 3,4,5, bin 2 = months 6,7,8 and bin 3 = months 9,10,11,12

# lionSight$bin[lionSight$SightingDate >= as.Date("2012-03-01") & 
#                     lionSight$SightingDate < as.Date("2012-06-01")] <- 1
# lionSight$bin[lionSight$SightingDate >= as.Date("2012-06-01") & 
#                     lionSight$SightingDate < as.Date("2012-09-01")] <- 2
# lionSight$bin[lionSight$SightingDate >= as.Date("2012-09-01") & 
#                     lionSight$SightingDate <= as.Date("2012-12-31")] <- 3
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2013-03-01") & 
#                     lionSight$SightingDate < as.Date("2013-6-01")] <- 4
# lionSight$bin[lionSight$SightingDate >= as.Date("2013-06-01") & 
#                     lionSight$SightingDate < as.Date("2013-09-01")] <- 5
# lionSight$bin[lionSight$SightingDate >= as.Date("2013-09-01") & 
#                     lionSight$SightingDate <= as.Date("2013-12-31")] <- 6
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2014-03-01") & 
#                     lionSight$SightingDate < as.Date("2014-06-01")] <- 7
# lionSight$bin[lionSight$SightingDate >= as.Date("2014-06-01") & 
#                     lionSight$SightingDate < as.Date("2014-09-01")] <- 8
# lionSight$bin[lionSight$SightingDate >= as.Date("2014-09-01") & 
#                     lionSight$SightingDate <= as.Date("2014-12-31")] <- 9
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2015-03-01") & 
#                     lionSight$SightingDate < as.Date("2015-06-01")] <- 10
# lionSight$bin[lionSight$SightingDate >= as.Date("2015-06-01") & 
#                     lionSight$SightingDate < as.Date("2015-09-01")] <- 11
# lionSight$bin[lionSight$SightingDate >= as.Date("2015-09-01") & 
#                     lionSight$SightingDate <= as.Date("2015-12-31")] <- 12
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2016-03-01") & 
#                     lionSight$SightingDate < as.Date("2016-06-01")] <- 13
# lionSight$bin[lionSight$SightingDate >= as.Date("2016-06-01") & 
#                     lionSight$SightingDate < as.Date("2016-09-01")] <- 14
# lionSight$bin[lionSight$SightingDate >= as.Date("2016-09-01") & 
#                     lionSight$SightingDate <= as.Date("2016-12-31")] <- 15
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2017-03-01") & 
#                     lionSight$SightingDate < as.Date("2017-06-01")] <- 16
# lionSight$bin[lionSight$SightingDate >= as.Date("2017-06-01") & 
#                     lionSight$SightingDate < as.Date("2017-09-01")] <- 17
# lionSight$bin[lionSight$SightingDate >= as.Date("2017-09-01") & 
#                     lionSight$SightingDate <= as.Date("2017-12-31")] <- 18
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2018-03-01") & 
#                     lionSight$SightingDate < as.Date("2018-06-01")] <- 19
# lionSight$bin[lionSight$SightingDate >= as.Date("2018-06-01") & 
#                     lionSight$SightingDate < as.Date("2018-09-01")] <- 20
# lionSight$bin[lionSight$SightingDate >= as.Date("2018-09-01") & 
#                     lionSight$SightingDate < as.Date("2018-12-31")] <- 21
# 
# lionSight$bin[lionSight$SightingDate >= as.Date("2019-03-01") & 
#                     lionSight$SightingDate < as.Date("2019-06-01")] <- 22
# lionSight$bin[lionSight$SightingDate >= as.Date("2019-06-01") & 
#                     lionSight$SightingDate < as.Date("2019-09-01")] <- 23
# lionSight$bin[lionSight$SightingDate >= as.Date("2019-09-01") & 
#                     lionSight$SightingDate < as.Date("2019-12-31")] <- 24

# More efficient binning
  # Annual bins
lionSight <- lionSight %>%
  mutate(
    year = year(SightingDate),
    bin = year - min(year) + 1
  )

  # Quarterly
# lionSight <- lionSight %>%
#   mutate(
#     year = year(SightingDate),
#     quarter = quarter(SightingDate),
#     bin = (year - 2008) * 4 + quarter
#   )

  # Monthly
lionSight_monthly <- lionSight %>%
  mutate(
    year = year(SightingDate),
    month = month(SightingDate),
    # Create monthly occasion: months since study start
    month_occasion = (year - min(year)) * 12 + month
  )

#check efficiency of this binning scheme
table(!is.na(lionSight$bin))
table(lionSight$bin)
hist(lionSight$bin)

# 3 3-4 month bins 
# drops 79 sightings out of 4270
# retains 98% of the data

#again, reshape the data, this time using the bins just created
temp2 <- melt(lionSight,id.var=c("AnimalID","bin"),measure.var="detect")
bin.histories = cast(temp2,AnimalID ~bin)

temp.monthly <- melt(lionSight_monthly, id.var = c("AnimalID", "month_occasion"), measure.var = "detect")
bin.histories.monthly <- cast(temp.monthly, AnimalID ~month_occasion)

# Remove NA row
bin.histories <- bin.histories %>% 
  filter(!is.na(AnimalID))
bin.histories.monthly <- bin.histories.monthly %>% 
  filter(!is.na(AnimalID))

# # remove "NA" column
# bin.histories <- select(bin.histories, -c(25))
# 
# ## Adding missing column for bin 4, which is absent because there were no observations, in the correct position
# bin.histories <- add_column(bin.histories, '4' = rep(0, length(bin.histories$AnimalID)), .after = '3')


#replace values of >=1 in bin.histories with 1 (i.e., detected is detected no matter how often)
#first need to save the IDs so they don't all get recoded, then put them back after recoding

  # better way
bin.histories <- bin.histories %>% 
  mutate(across(where(is.numeric), ~ifelse(.x > 1, 1, .x))) %>% 
  as.data.frame()

bin.histories.monthly <- bin.histories.monthly %>% 
  mutate(across(where(is.numeric), ~ifelse(.x > 1, 1, .x))) %>% 
  as.data.frame()
# ID <- bin.histories$AnimalID
# bin.histories[bin.histories >= 1]<-1
# bin.histories$AnimalID <- ID

# # Remove other species
# bin.histories <- bin.histories[-c(1:17), ]

# add pack information to bin.histories to allow subsetting the packs with sufficient monitoring for
# year-specific density estimates
# 
# add.packs <- subset.data.frame(x=lionSight, select = c(AnimalID, AnimalGroup))
# add.packs <- distinct(add.packs)
# bin.histories <- left_join(x= bin.histories, y=add.packs)
# 
# # check the exact spelling and capitalization for pack (=AnimalGroup) names to allow subsetting
# # note - must combine some of them that appear more than one way
# levels(as.factor(bin.histories$AnimalGroup))
# 
# bin.histories.2016 <- subset(bin.histories, AnimalGroup %in% c('Tateyoyo','Kasabushi','Chunga', 'Chunga Pack'))
# #confirm subset is correct packs
# levels(as.factor(bin.histories.2016$AnimalGroup))
# #any dog that changed packs in its entire life has multiple rows with identical capture histories, once for each pack
# #pack ID does not matter for this analysis so just remove duplicated rows, keeping in mind that they are NOT identical until
# #pack is deleted
# bin.histories.2016 <- bin.histories.2016[-26]
# bin.histories.2016 <- distinct(bin.histories.2016)
# 
# bin.histories.2017 <- subset(bin.histories, AnimalGroup %in% c('Lushimba','Tateyoyo','Kasabushi','Chunga', 'Chunga Pack','Mapunga Pack','twin palm'))
# levels(as.factor(bin.histories.2017$AnimalGroup))
# bin.histories.2017 <- bin.histories.2017[-26]
# bin.histories.2017 <- distinct(bin.histories.2017)
# 
# 
# bin.histories.2018 <- subset(bin.histories, AnimalGroup %in% c('Lushimba','Chunga', 'Chunga Pack','twin palm','Shishamba'))
# levels(as.factor(bin.histories.2018$AnimalGroup))
# bin.histories.2018 <- bin.histories.2018[-26]
# bin.histories.2018 <- distinct(bin.histories.2018)
# 
# bin.histories.2019 <- subset(bin.histories, AnimalGroup %in% c('Chunga', 'Chunga Pack','twin palm','Shishamba',"Musanza","Chinese Bridge Pack","Lufupa","Musekwa", "Katinti Pack","Panthera Dispersers"))
# levels(as.factor(bin.histories.2019$AnimalGroup))
# bin.histories.2019 <- bin.histories.2019[-26]
# bin.histories.2019 <- distinct(bin.histories.2019)

# concatenate the bins into a capture history named ch, (THIS WAS DESIGNED to format for RMark code and MARK, but code just below
# provides format for JAGS code) 
# so comment out the one not needed

#function to read the matrix and create the capture history strings
past.ch<-function(x)
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    y<-(x[i,]>0)*1

    out[i]<-paste(y[1],y[2],y[3],y[4],y[5],
                  y[6],y[7],y[8],y[9],
                  y[10],y[11],y[12],y[13],y[14],
                  y[15],y[16],y[17],y[18],
                  y[19],y[20],y[21],y[22],y[23],y[24], sep="")
  }
  return(out)
}

#Use paste.ch() just created to create an RMARK ch file (will need to tweak function
# and call to match number of time bins - note offset in column indices to account
# for AnimalID in column 1 of bin.histories)
# Also note the last column named NA in bin histories (from sightings in periods not included
# in the selected date bins) is discarded here, so it is not a problem in subsequent analysis.

## Keep bin histories separated by columns for bayesian models in JAGS following Kery and Schuab examples

capt.hist<-data.frame(ch=(bin.histories[,2:17]))

# capt.hist.2016 <- data.frame(ch=(bin.histories.2016[,2:25]))
# capt.hist.2017 <- data.frame(ch=(bin.histories.2017[,2:25]))
# capt.hist.2018 <- data.frame(ch=(bin.histories.2018[,2:25]))
# capt.hist.2019 <- data.frame(ch=(bin.histories.2019[,2:25]))
# 
# write.csv(capt.hist, "capt_hist.csv")
# write.csv(capt.hist.2016, "capt_hist_2016.csv")
# write.csv(capt.hist.2017, "capt_hist_2017.csv")
# write.csv(capt.hist.2018, "capt_hist_2018.csv")
# write.csv(capt.hist.2019, "capt_hist_2019.csv")

# Bundle data for JAGS
y <- as.matrix(bin.histories[,-1])  # Remove AnimalID column
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



# Age and sex data --------------------------------------------------------
# Get sex for each individual 
sex <- lionSight %>% 
  filter(AnimalID %in% bin.histories$AnimalID) %>% 
  select(AnimalID, Sex) %>%
  arrange(AnimalID) %>%
  distinct()

# Get age at first detection for CMR models by comparison of DOB and date of first detection

#year of first detection 
age_at_first <- data.frame(AnimalID = bin.histories$AnimalID,
                           sex = sex$Sex)
age_at_first$year_first_cap <- apply(bin.histories[,-1], 1, get_first) + 2007

# Add DOB and DOB error
age_at_first <- age_at_first %>% 
  left_join(lionSight %>% select(AnimalID, DOB, `error DOB`), by = "AnimalID") %>% 
  mutate(DOB = as_date(DOB),
         midp_first_cap = as_date(paste0(year_first_cap, "-07-01")),
         year = year(DOB),
         age_at_first = midp_first_cap-DOB) %>% 
  distinct()
  

# Decide how to handle DOB error - assign corrected DOBs
age_at_first <- age_at_first %>%
  mutate(DOBcorrected = case_when(
    `error DOB` <= 30 ~ DOB,
    `error DOB` > 30 & `error DOB` <= 90 ~ DOB + 45,
    `error DOB` > 90 & `error DOB` <=180 ~ DOB + 135,
    `error DOB` > 180 & `error DOB` <= 365 ~ DOB + 270,
    `error DOB` > 270 ~ DOB
    )
  )

# Add ages in days for all years
# Function to assign age class based on age in years
get_age_class <- function(age_years) {
  case_when(age_years < 2 ~ 1,                # cub (0-1.99)
            between(age_years, 2, 3.99) ~ 2,  # subadult (2-3.99)
            between(age_years, 4, 5.99) ~ 3,  # young adult (4-5.99)
            age_years >= 6 ~ 4)               # old adult (6+)
}

# Convert age at first capture from days to years
age_at_first$age_at_first_years <- as.numeric(age_at_first$age_at_first, units = "days") / 365.25

# For unknown age, known sex: assign median age of non-cubs
# For unknown age, unknown sex: make 1st year cubs
# FOr unknown sex, assign random
median.age.adult <- age_at_first %>% 
  filter(age_at_first_years >= 2) %>% 
  pull(age_at_first_years) %>% 
  median(na.rm = TRUE)
age_at_first <- age_at_first %>%
  mutate(age_at_first_years = case_when(
    is.na(age_at_first_years) & sex != "Unknown" ~ median.age.adult,
    is.na(age_at_first_years) & sex == "Unknown" ~ 1,
    .default = age_at_first_years
    ),
    sex = case_when(
      sex == "Unknown" ~ sample(c("Male", "Female"), 1, replace = TRUE, prob = c(0.5, 0.5)),
      .default = sex
    ))

# Get the range of years in your study
all_years <- sort(unique(lionSight$year))  # Adjust to match your nyears=16
nyears <- length(all_years)
nind <- nrow(age_at_first)

# Initialize age matrix
age <- matrix(NA, nrow = nind, ncol = nyears)

# Fill in age class for each individual at each time point
for (i in 1:nind) {
  first_year <- age_at_first$year_first_cap[i]
  first_year_index <- which(all_years == first_year)
  first_year_age <- age_at_first$age_at_first_years[i]
  
  for (t in first_year_index:nyears) {
    # Calculate age at time t
    years_since_first <- t - first_year_index
    current_age <- first_year_age + years_since_first
    
    # Assign age class
    age[i, t] <- get_age_class(current_age)
  }
}

age.mat <- createAge(f = first, age = get_age_class(age_at_first$age_at_first_years), nyears = nyears, mAge = 4)

#Month of first detection, assigned to the middle of each bin/window in the ch's (e.g. May would be middle of April, May, June bin)
#NOTE need a number of if()s equal to number of bins/yr, and
#number of logicals within each if equal to number of years, and
#will need to tweak if ch is not may/jun jul/aug sept/oct within each year
# for (i in 1:length(capt.hist$ch)) 

## 3 3-mo bins 8 years for a total of 24 bins. --BG

for (i in 1:length(capt.hist$ch)){

  if (regexpr("1",capt.hist$ch[i]) == 1 |regexpr("1",capt.hist$ch[i]) == 4| regexpr("1",capt.hist$ch[i]) == 7|
      regexpr("1",capt.hist$ch[i]) == 10 |regexpr("1",capt.hist$ch[i])==13|regexpr("1",capt.hist$ch[i])==16|regexpr("1",capt.hist$ch[i])==19 | regexpr("1",capt.hist$ch[i])==22) capt.hist$monthatfirst[i] <- 5

  if (regexpr("1",capt.hist$ch[i]) == 2 |regexpr("1",capt.hist$ch[i]) == 5| regexpr("1",capt.hist$ch[i]) == 8|
      regexpr("1",capt.hist$ch[i]) == 11 | regexpr("1",capt.hist$ch[i])==14|regexpr("1",capt.hist$ch[i])==17|regexpr("1",capt.hist$ch[i])==20 | regexpr("1",capt.hist$ch[i])==23) capt.hist$monthatfirst[i] <- 8

  if (regexpr("1",capt.hist$ch[i]) == 3 |regexpr("1",capt.hist$ch[i]) == 6| regexpr("1",capt.hist$ch[i]) == 9|
      regexpr("1",capt.hist$ch[i]) == 12|regexpr("1",capt.hist$ch[i])==15|regexpr("1",capt.hist$ch[i])==18|regexpr("1",capt.hist$ch[i])==21 | regexpr("1",capt.hist$ch[i])==24) capt.hist$monthatfirst[i] <- 11
}

#^ 5, 8, and 11 correspond to the middle of each of the 3 month bins (e.g. in this example it is may, august, and november)


## 2 6-mo bins --BG
# for (i in 1:length(capt.hist$ch)) {
#
#   if (regexpr("1",capt.hist$ch[i]) == 1 |regexpr("1",capt.hist$ch[i]) == 3| regexpr("1",capt.hist$ch[i]) == 5|regexpr("1",capt.hist$ch[i]) == 7 |regexpr("1",capt.hist$ch[i]) == 9|
#       regexpr("1",capt.hist$ch[i]) == 11 |regexpr("1",capt.hist$ch[i]) == 13) capt.hist$monthatfirst[i] <- 3
#
#   if (regexpr("1",capt.hist$ch[i]) == 2 | regexpr("1",capt.hist$ch[i]) == 4| regexpr("1",capt.hist$ch[i]) == 6 | regexpr("1",capt.hist$ch[i]) == 8 | regexpr("1",capt.hist$ch[i]) == 10|
#       regexpr("1",capt.hist$ch[i]) == 4| regexpr("1",capt.hist$ch[i]) == 12| regexpr("1",capt.hist$ch[i]) == 1) capt.hist$monthatfirst[i] <- 9
# }


#RENAME DOB COLUMN
capt.hist$DOB<-capt.hist$'The date the animal is born in (DD/MM/YYYY), if only a year is known note: yyyy/01/01'

capt.hist$monthatfirst
#convert DOB in master ID file to an R date (NOTE POSIXlt, not DATE)
capt.hist$DOB<-strptime(capt.hist$DOB,"%m/%d/%Y")


# estimate age at first detection using the month and year of first detection and DOB.

# Assume day of month at first detection was the 15th to avoid negative ages.
capt.hist$ageatfirst<-as.numeric(difftime(strptime(paste(capt.hist$yearatfirst,
                                                         capt.hist$monthatfirst,"28",sep="-"),"%Y-%m-%d"),capt.hist$DOB,
                                          units="days")/365)
#a few animals have 1990 DOB presumably to indicate they were in oldest age class when detected, e.g. 10 years old
capt.hist$ageatfirst[which(capt.hist$ageatfirst > 10)] <- 10
capt.hist
