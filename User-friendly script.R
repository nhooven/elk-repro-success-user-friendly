# Title: Using low-fix rate GPS telemetry to expand estimates of ungulate reproductive success
# Subtitle: User-friendly R code
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 13 Sep 2021
# Date completed: 13 Sep 2021
# Date modified: 
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(sp)             # work with spatial objects
library(amt)            # work with tracks and generate steps
library(lubridate)      # work with time 
library(adehabitatHR)   # fit MCPs
library(mefa4)          # notin function
library(randomForest)   # random forest classification

#_____________________________________________________________________________________________________________
# 0. Inputs ----
#_____________________________________________________________________________________________________________

# these are meant to be modified by the analyst to reflect differences in data collection, species biology, etc.

# define projection
projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")

# expected collar sampling rate (in h)
fix.rate <- 13

# temporal tolerance (in h; how many hours +/- the sampling rate should be considered consecutive?)
temp.tol <- 1

# window size (in days) for mean step length (how many days shold be included for calculations?)
sl.period <- 3

# window size (in days) for minimum convex polygon 
mcp.period <- 7

# time zone
time.zone <- "America/New_York"

# beginning of study period
start.date <- as.POSIXct("2021-05-15 00:00:00", tz = time.zone)

# end of study period
end.date <- as.POSIXct("2021-07-15 23:59:59", tz = time.zone)

#_____________________________________________________________________________________________________________
# 1. Raw data cleaning ----
#_____________________________________________________________________________________________________________

# this section requires timestamped GPS relocation data that has already been screened for inaccurate relocations

# read in data
raw.data <- read.csv("all_data.csv")

# only keep columns we need
raw.data <- raw.data %>% dplyr::select(x, y, t, AnimalID)

# make sure "t" is a POSIXct variable
raw.data$t <- as.POSIXct(raw.data$t)

# vector of AnimalIDs
animal.ids <- unique(raw.data$AnimalID)

# create bursts
all.data <- data.frame()

for (y in animal.ids) {
  
  indivID <- y
  
  # subset data
  indiv.data <- raw.data %>% filter(AnimalID == indivID)
  
  # make a track
  indiv.track <- indiv.data %>% make_track(x, y, t, all_cols = TRUE)
  
  # create "burst_" column (resample to fix rate +/- tolerance)
  indiv.track.res <- indiv.track %>% track_resample(rate = hours(fix.rate), tolerance = hours(temp.tol))
  
  # bind to "all.data" data.frame
  all.data <- rbind(all.data, indiv.track.res)
  
}

#_____________________________________________________________________________________________________________
# 2. Moving window mean step length ----
#_____________________________________________________________________________________________________________

all.steps <- data.frame()

for (z in unique(all.data$AnimalID)) {
  
  indivID.2 <- z
  
  # subset all.data
  indiv.data.2 <- all.data %>% filter(AnimalID == indivID.2)
  
  # create a track
  indiv.track.2 <- indiv.data.2 %>% make_track(x_, y_, t_, all_cols = TRUE)
  
  # generate steps
  indiv.steps <- indiv.track.2 %>% steps_by_burst()
  
  # determine how many days before and after period we need to keep based upon the window size
  buffer.days <- (sl.period - 1) / 2
  
  # filter steps across study period
  indiv.steps.1 <- indiv.steps %>% filter(t1_ < (end.date + buffer.days*24*60*60) & 
                                          t1_ >= (start.date - buffer.days*24*60*60))
  
  # create a day of the year variable
  day.seq <- seq(start.date, end.date, by = 24*60*60)
  
  # note that we add an extra day because of Daylight Saving Time in the U.S.
  DOY.seq <- as.integer(difftime(day.seq, 
                                 as.POSIXct("2021-01-01 00:00:00", 
                                            tz = time.zone), 
                                 units = "days") + 2)

  # add a "DOY" variable to the steps tibble
  indiv.steps.2 <- indiv.steps.1 %>% mutate(DOY = as.integer(difftime(t1_, 
                                                                      as.POSIXct("2021-01-01 00:00:00", 
                                                                                 tz = time.zone), 
                                                                      units = "days") + 2))
  
  # compute mean step length within a 3-day moving window
  indiv.steps.2.summary <- data.frame()
  
  # for loop which calculates 3-day averages of average daily sl
  for (w in DOY.seq) {
    
    # subset data
    focal.steps <- indiv.steps.2 %>% filter(DOY %in% c(w - buffer.days, w, w + buffer.days))
    
    # calculate mean sl for focal period
    focal.mean <- mean(focal.steps$sl_, na.rm = TRUE)
    
    # bind into a df with the DOY
    focal.summary <- data.frame(Animal = indivID.2,
                                sl.3day = focal.mean,
                                DOY = w)
    
    # bind to master df
    indiv.steps.2.summary <- rbind(indiv.steps.2.summary, focal.summary)
    
  }
  
  # bind to master df
  all.steps <- rbind(all.steps, indiv.steps.2.summary)
  
}

#_____________________________________________________________________________________________________________
# 3. Moving window MCPs ----
#_____________________________________________________________________________________________________________

all.mcps <- data.frame()

for (v in unique(all.data$AnimalID)) {
  
  indivID.3 <- v
  
  # subset all.data to individual
  indiv.data.3 <- all.data %>% filter(AnimalID == indivID.3)
  
  # subset for study period +/- days needed for calculations
  mcp.days <- (mcp.period - 1) / 2
  
  indiv.data.4 <- indiv.data.3 %>% filter(t_ < (end.date + mcp.days*24*60*60) & 
                                          t_ >= (start.date - mcp.days*24*60*60))
  
  # create a DOY column (as in part 2)
  indiv.data.4 <- indiv.data.4 %>% mutate(DOY = as.integer(difftime(t_, 
                                                                    as.POSIXct("2021-01-01 00:00:00", 
                                                                               tz = time.zone), 
                                                                      units = "days") + 2))
  
  # calculate 100% MCP areas and add to data frame
  Win.MCP <- data.frame(MCP = NA,
                        DOY = DOY.seq)
  
  for (p in DOY.seq) {
    
    # define relocations for moving window p
    indiv.data.5 <- indiv.data.4 %>% dplyr::filter(DOY >= (p - mcp.days) & DOY <= (p + mcp.days))
    
    # fit MCP to relocations
    focal.sp <- SpatialPoints(coords = indiv.data.5[ ,c("x_", "y_")], 
                              proj4string = projection)
    
    # fit an MCP (use an NA if there are not enough relocations)
    focal.mcp <- ifelse(nrow(focal.sp@coords) > 4,
                        mcp.area(focal.sp, percent = 100, unin = "m", unout = "km2", plotit = FALSE),
                        NA)
    
    Win.MCP$MCP[Win.MCP$DOY == p] <- ifelse(is.na(focal.mcp) == FALSE,
                                            focal.mcp[[1]],
                                            focal.mcp)
    
  }
  
  # add CollarID and bind to master data frame
  Win.MCP <- Win.MCP %>% mutate(Animal = indivID.3)
  
  # bind to master df
  all.mcps <- rbind(all.mcps, Win.MCP)
  
}

# merge steps and MCPs dfs
all.metrics <- merge(all.steps, all.mcps)

#_____________________________________________________________________________________________________________
# 4. Add in "Case" variable from confirmed parturition dates ----
#_____________________________________________________________________________________________________________

# read in parturition dates file
part.dates <- read.csv("part_dates.csv")

all.metrics.2 <- data.frame()

for (q in unique(all.metrics$Animal)) {
  
  indivID.4 <- q
  
  # subset individual's parturition date
  indiv.date <- part.dates$DOY[part.dates$Animal == indivID.4]
  
  # subset individual's data only and add a "Case" column (0 or 1)
  indiv.metrics <- all.metrics %>% filter(Animal == indivID.4) %>%
                                   mutate(Case = ifelse(DOY == indiv.date, 1, 0))
  
  # bind to master df
  all.metrics.2 <- rbind(all.metrics.2, indiv.metrics)
  
}

#_____________________________________________________________________________________________________________
# 5. Examine distributions of predictors based upon Case ----
#_____________________________________________________________________________________________________________

# sl.3day
ggplot(data = all.metrics.2, aes(x = sl.3day)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case)),
                    alpha = 0.5)

# mcp
ggplot(data = all.metrics.2, aes(x = MCP)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case)),
                    alpha = 0.5)

#_____________________________________________________________________________________________________________
# 6. Random forest classification ----
#_____________________________________________________________________________________________________________

# train RF model
rf.model <- randomForest(as.factor(Case) ~ sl.3day + MCP,
                                     na.action = na.omit,
                                     sampsize = c(50, 50),
                                     ntree = 1000,
                                     data = all.metrics.2)

# confusion matrix
rf.model

# assess variable importance
importance(rf.model, type = 2)

# generate predictions for each focal day in the time series
predictions <- as.data.frame(predict(rf.model, type = "prob"))

pred.prob <- predictions[ ,2]

# subset only the complete cases (i.e. rows without NA values for predictors)
all.metrics.complete <- all.metrics.2[complete.cases(all.metrics.2), ]

# bind probabilities to complete cases
all.metrics.complete <- cbind(all.metrics.complete, pred.prob)

# plot all time series
ggplot(data = all.metrics.complete, aes(DOY, pred.prob)) +
       theme_bw() +
       facet_wrap(~Animal) +
       geom_point(color = "darkgreen", alpha = 0.5) +
       ylab("") +
       xlab("") +
       ggtitle("Probability time series")

# save model for testing
save(rf.model, file = "rf_model.Rdata")
