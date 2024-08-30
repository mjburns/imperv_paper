## Preparation of LSC data for paper XXX
## Script written by co-author MJB
## August, 2024

## Load packages and set timezone
library(lubridate) 
library(tidyverse)

Sys.setenv(TZ = "UTC")

# The rainfall runoff data used in Li et al. (XXX) is loaded here
load(paste(here::here(""), "/Input/Flow_Data_Sampled.rda", sep = ""))
head(flow_data_to_sample)

## Only select the LSC catchments with SCMs, so exclude the controls and reference
## "Br" "D4" "Fe" "L1" "L4" "Ln" "Ls" "Ly" "Ol" "Sa"
flow_data_to_sample_screened <-  flow_data_to_sample[(!flow_data_to_sample$new_sites_name %in% (c("Br", "Fe", "Ly", "Ol", "Sa"))),]

## Specify max and min precip depths to include in analysis (match data units), rainfall in mm
precip_min <- 2.54
precip_max <- 999999999999

## For this paper, screen the lsc data to only get rainfall runoff data from 2014 onwards
## Roughly the time point by which most SCMs were installed
#date_start <- ymd_hms("2000-01-01 00:00:00", tz = "UTC")
date_start <- ymd_hms("2014-01-01 00:00:00", tz = "UTC")
date_end <- max(flow_data_to_sample_screened$RainEnd_UTC)

## Format the lsc rainfall-runoff data similar to Clarkesburg
data.df <- data.frame(Precip_EndTime = flow_data_to_sample_screened$RainEnd_UTC, Site = flow_data_to_sample_screened$new_sites_name,
                      Precip_Total_mm = flow_data_to_sample_screened$event_rain_mm, Q_Runoff_mm = flow_data_to_sample_screened$new_QF_mm)

## Make into a tibble?
data.df.tibble <- tibble(data.df)

## Filter based on date and rainfall

data <- data.df.tibble %>%  
  filter(Precip_EndTime >= date_start & Precip_EndTime <= date_end) %>%
  filter(Precip_Total_mm >= precip_min & Precip_Total_mm <= precip_max) %>%
  select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)

## Save the lsc events

save(data, file = paste(here::here(""), "/Input/lsc_events_no_precip_max.Rdata", sep = ""))
##save(data, file = paste(here::here(""), "/Input/lsc_events_no_precip_max_all_data.Rdata", sep = ""))

## Get the LSC land-use data
landuse <- get(load("~/Documents/imperv_paper/Input/ei_ts_11_sites_correct.rda"))
landuse$lubri_date <- ymd(landuse$date, tz = "UTC")

#levels(factor(landuse$sitecode))

## Screen landuse to only get the sites of interest for the impervious paper
## Get lower and upper bound based on date_start and date_end
landuse_screened <- landuse[(!landuse$sitecode %in% (c("BRS0015", "DBS0008", "FER0006", "LYR0007", "OLN0009", "SAS0002"))),]

## Floor date to day
date_start_day <- floor_date(date_start, "day")
date_end_day <- floor_date(date_end, "day")

## Declare null variable
landuse_postscm <- NULL

for(i in 1:NROW(levels(factor(landuse_screened$sitecode)))) {
  
  get_site_i_data <- landuse[which(landuse_screened$sitecode == levels(factor(landuse_screened$sitecode))[i]),]
  
  ## Screen by dates
  get_site_i_post_scm_early <- get_site_i_data[which(get_site_i_data$lubri_date == date_start_day),]
  get_site_i_post_scm_later <- get_site_i_data[which(get_site_i_data$lubri_date == date_end_day),]

  ## Build a data.frame for the sites
  landuse_i_df <- data.frame(sitecode = get_site_i_post_scm_early$sitecode[1], ti_early = get_site_i_post_scm_early$ti[1], ti_later = get_site_i_post_scm_later$ti[1],
                             ei_early = get_site_i_post_scm_early$ei[1], ei_later = get_site_i_post_scm_later$ei[1], s_early = get_site_i_post_scm_early$ro[1],
                             s_later = get_site_i_post_scm_later$ro[1])
  
  #Combine into the above
  landuse_postscm <- rbind(landuse_postscm, landuse_i_df) } #End the loop now


## Get the same site codes as the paper
## "Br" "D4" "Fe" "L1" "L4" "Ln" "Ls" "Ly" "Ol" "Sa"

sites_match <- data.frame(sitecode = c("LIS0001", "LSN0001", "LSS0001", "LIS0004", "DBS0004"), shorter_name = c("L1", "Ln", "Ls", "L4", "D4"))


## Use a match
landuse_postscm$imperv_sitecode <-  sites_match$shorter_name[match(landuse_postscm$sitecode, sites_match$sitecode)]


plot(c(1,1), c(landuse_postscm$ti_early[1], landuse_postscm$ti_later[1]), type = 'o', ylim = c(0,0.3), axes = FALSE, xlim = c(-1,6), xlab = "Site", ylab = "Metric")
lines(c(1.1,1.1), c(landuse_postscm$ei_early[1], landuse_postscm$ei_later[1]), col = 'red', type = 'o')
lines(c(1.2,1.2), c(landuse_postscm$s_early[1], landuse_postscm$s_later[1]), col = 'blue', type = 'o')

lines(c(2.0,2.0), c(landuse_postscm$ti_early[2], landuse_postscm$ti_later[2]), col = 'black', type = 'o')
lines(c(2.1,2.1), c(landuse_postscm$ei_early[2], landuse_postscm$ei_later[2]), col = 'red', type = 'o')
lines(c(2.2,2.2), c(landuse_postscm$s_early[2], landuse_postscm$s_later[2]), col = 'blue', type = 'o')


lines(c(3.0,3.0), c(landuse_postscm$ti_early[3], landuse_postscm$ti_later[3]), col = 'black', type = 'o')
lines(c(3.1,3.1), c(landuse_postscm$ei_early[3], landuse_postscm$ei_later[3]), col = 'red', type = 'o')
lines(c(3.2,3.2), c(landuse_postscm$s_early[3], landuse_postscm$s_later[3]), col = 'blue', type = 'o')


lines(c(4.0,4.0), c(landuse_postscm$ti_early[4], landuse_postscm$ti_later[4]), col = 'black', type = 'o')
lines(c(4.1,4.1), c(landuse_postscm$ei_early[4], landuse_postscm$ei_later[4]), col = 'red', type = 'o')
lines(c(4.2,4.2), c(landuse_postscm$s_early[4], landuse_postscm$s_later[4]), col = 'blue', type = 'o')


lines(c(5.0,5.0), c(landuse_postscm$ti_early[5], landuse_postscm$ti_later[5]), col = 'black', type = 'o')
lines(c(5.1,5.1), c(landuse_postscm$ei_early[5], landuse_postscm$ei_later[5]), col = 'red', type = 'o')
lines(c(5.2,5.2), c(landuse_postscm$s_early[5], landuse_postscm$s_later[5]), col = 'blue', type = 'o')


axis(1, at = c(-1, 1, 2, 3, 4, 5, 6), lab = c("", landuse_postscm$imperv_sitecode, ""))
axis(2)
legend("topleft", c("ti", "ei", "s"), col = c('black', 'red', 'blue'), lty = c(1,1,1))

## Avg
landuse_postscm$ti_mean <- (landuse_postscm$ti_early + landuse_postscm$ti_later)/2
landuse_postscm$ei_mean <- (landuse_postscm$ei_early + landuse_postscm$ei_later)/2
landuse_postscm$s_mean <- (landuse_postscm$s_early + landuse_postscm$s_later)/2

## Add new site names and round
## Table 3 Melb data
order_paper <- c("D4", "L4", "L1", "Ln", "Ls")
landuse_postscm_ordered <- landuse_postscm[match(order_paper, landuse_postscm$imperv_sitecode),]

#Generate table 3
table_three_melb <- data.frame(Treatment = landuse_postscm_ordered$imperv_sitecode, method_one_es = round((landuse_postscm_ordered$s_mean * 100),1), method_two_ei = c(0,0,0,0,0), ti = round((landuse_postscm_ordered$ti_mean * 100),1))



## Checking paired catchment work
## Based on discussions with Chris and Peter Poelsma
## This text can be removed at some point
#control catchment area total = 2.21 ha
#Impact is 2.62 ha
#pair <- read.csv("~/Documents/imperv_paper/Input/Paired.csv", stringsAsFactors = FALSE, header = TRUE, na.strings = -9999)
#pair$control_mm <- ((pair$control_L/1000)/(2.21*10000))*1000
#pair$impact_mm <- ((pair$impact_L/1000)/(2.62*10000))*1000

#plot(pair$rainfall_mm, pair$control_mm, col = 'red')
#points(pair$rainfall_mm, pair$impact_mm, col = 'blue')
#legend("bottomright", c("Control (Bradman)", "Impact (Heath)"), col = c('red', 'blue'), pch = c(1,1))

#lm_control <- lm(control_mm ~ rainfall_mm, data = pair)
#lm_impact <- lm(impact_mm ~ rainfall_mm, data = pair)


