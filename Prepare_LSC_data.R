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
       
     



