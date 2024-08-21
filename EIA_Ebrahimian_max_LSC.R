## Estimate Effective Impervious Area Using Rainfall-Runoff Data
## Based on Ebrahimian and others (2016) Weighted Linear Regression

# Charlie Stillwell
# July 11, 2024
# Updated 26th July, 2024 by MJB to include LSC catchments
 


### Step 00: Set Up ------------------------------------------------------------

# Load packages
library(tidyverse)
library(broom)
library(lubridate) # MJB edit
library(gridExtra) # MJB edit

#Set the timezone
Sys.setenv(TZ = "UTC")

#Here also load the flow data to be sampled, name is flow_data_to_sample
load(paste(here::here(""), "/Flow_Data_Sampled.rda", sep = ""))
head(flow_data_to_sample)

#As discussed, probably should remove Ly and Ol, basically very little impervious area in these catchments
#levels(factor(flow_data_to_sample$new_sites_name))
#"Br" "D4" "Fe" "L1" "L4" "Ln" "Ls" "Ly" "Ol" "Sa"

flow_data_to_sample_screened <-  flow_data_to_sample[(!flow_data_to_sample$new_sites_name %in% (c("Ly", "Ol"))),]
#Think above looks OK


# Specify max and min precip depths to include in analysis (match data units), rainfall in mm
# Updated to essentially have no maximum precip as per discussions
precip_min <- 2.54
precip_max <- 999999999999

# Choose types of events (only runoff-generating events or all precip events)
event_type <- "all"     # one of "all" or "nonzero"

# Specify date range of events to include in analysis
# (for my data, use precipitation end date-time)
# date_start <- as.Date("2016-01-01")
# date_end <- as.Date("2018-12-31")

# Based on LSC stuff
# Tweak based on starting the analysis from 2014, after SCMs implemented?
#date_start <- min(flow_data_to_sample_screened$RainEnd_UTC)
date_start <- ymd_hms("2014-01-01 00:00:00", tz = "UTC")

date_end <- max(flow_data_to_sample_screened$RainEnd_UTC)

# Read in rainfall-runoff data
# And Make similar columns to below.
# Make our data similar to Clarkesburg
data.df <- data.frame(Precip_EndTime = flow_data_to_sample_screened$RainEnd_UTC, Site = flow_data_to_sample_screened$new_sites_name,
                      Precip_Total_mm = flow_data_to_sample_screened$event_rain_mm, Q_Runoff_mm = flow_data_to_sample_screened$new_QF_mm)

#Make into a tibble?
data.df.tibble <- tibble(data.df)

# Read in rainfall-runoff data
#data <- 
 # read_csv("Input/Clarksburg_StormEvents_2004_2018.csv", 
  #         col_types = cols(Precip_Duration_hrs = col_skip(), 
   #                         Precip_EndTime = 
    #                          col_datetime(format = "%m/%d/%Y %H:%M"), 
     #                       Precip_PriorDry_hrs = col_skip(), 
      #                      Precip_StartTime = col_skip(), 
       #                     Q_Duration_hrs = col_skip(), 
        #                    Q_EndTime = col_skip(), 
         #                   Q_PeakTime = col_skip(), 
          #                  Q_Peak_cfs = col_skip(), 
           #                 Q_Peak_cms.sqkm = col_skip(), 
            #                Q_PriorDry_hrs = col_skip(), 
             #               Q_RiseRate_cms.hr = col_skip(),  
              #              Q_StartTime = col_skip(), 
               #             Q_TimeToFlowPeak_hrs = col_skip(), 
                #            Q_Total_PrecipRatio = col_skip())) %>%
 #

data <- data.df.tibble %>%  
  filter(Precip_EndTime >= date_start & Precip_EndTime <= date_end) %>%
   filter(Precip_Total_mm >= precip_min & Precip_Total_mm <= precip_max) %>%
    select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
  
#head(data)
#I think ok to export this to Charlie
save(data, file = "~/Documents/Git/lsc_flow/DCI/lsc_events.Rdata")


### Step 01: Plot Rainfall-Runoff Data -----------------------------------------

plot_list <- list() #MJB to plot rainfall-runoff
#Edited y scale so that it plots up to max of runoff in mm

for (i in unique(data$Site)) {
  plot <- data %>%
    filter(Site == i) %>%
    ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point() +
    scale_x_continuous(limits = c(0, max(data$Precip_Total_mm))) +
    scale_y_continuous(limits = c(0, max(data$Q_Runoff_mm))) +
    labs(title = i) +
    theme_bw()
  plot_list[[i]] <- plot
  #print(plot)
}
rm(i, plot)

#Try plotting, MJB, rainfall-runoff stuff
grid.arrange(grobs=plot_list,ncol=3)



### Step 02: Prelim. Ordinary Least Squares to Remove Outliers -----------------

# Assign threshold for outlier detection (little rain, high runoff... 
# suggesting that rainfall or runoff data may be incorrect)
outlier_std_resid <- -2

# Loop through sites and fit ordinary least squares to filter out outliers
events_all <- tibble()
events_nonzero <- tibble()

for (i in unique(data$Site)) {
  
  # All rainfall events, including events with 0 mm runoff
  data_all_i <- data %>% filter(Site == i)
  ols_all_i <- lm(Q_Runoff_mm ~ Precip_Total_mm, data = data_all_i)
  ols_all_i_aug <- augment(ols_all_i, data = data_all_i)
  events_all_i <- ols_all_i_aug %>%
    filter(.std.resid > outlier_std_resid) %>%
    select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
  events_all <- bind_rows(events_all, events_all_i)
  
  # Only rainfall events that produced runoff > 0 mm
  data_nonzero_i <- data %>%
    filter(Site == i) %>%
    filter(Q_Runoff_mm > 0)
  ols_nonzero_i <- lm(Q_Runoff_mm ~ Precip_Total_mm, data = data_nonzero_i)
  ols_nonzero_i_aug <- augment(ols_nonzero_i, data = data_nonzero_i)
  events_nonzero_i <- ols_nonzero_i_aug %>%
    filter(.std.resid > outlier_std_resid) %>%
    select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
  events_nonzero <- bind_rows(events_nonzero, events_nonzero_i)
}
rm(i, data_all_i, data_nonzero_i, ols_all_i, ols_nonzero_i, 
   ols_all_i_aug, ols_nonzero_i_aug, events_all_i, events_nonzero_i)



### Step 03: Successive Weighted Linear Regression -----------------------------

# Conduct successive weighted linear regression
wls_results <- tibble()
eia_events <- tibble()

for (i in unique(data$Site)) {
  
  # Select event types to include in analysis
  if (event_type == "all") {
    events <- events_all %>% filter(Site == i)
  } else {
    events <- events_nonzero %>% filter(Site == i)
  }
  
  num_combined_events <- 1     # simply to initialize while loop
  while (num_combined_events > 0) {
    
    # Ordinary least squares, then ordinary least squares on residuals vs precip
    ols <- lm(Q_Runoff_mm ~ Precip_Total_mm, data = events)
    ols <- augment(ols, data = events)
    ols_resid <- lm(.resid^2 ~ 0 + Precip_Total_mm, data = ols)
    ols_resid <- augment(ols_resid, data = ols)
    
    # Compute weights; inverse of fitted values from OLS on residuals vs precip
    weights <- ols_resid %>%
      mutate(Weight = 1/.fitted) %>%
      select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm, Weight)
    
    # Conduct weighted linear regression and find standard error
    wls <- lm(Q_Runoff_mm ~ Precip_Total_mm, 
                data = weights, weights = Weight)
    wls_slope <- wls$coefficients[[2]]
    wls_yint <- wls$coefficients[[1]]
    wls <- augment(wls, data = weights)
    wls_se <- sqrt(sum(wls$Weight * wls$.resid^2) / (nrow(wls) - 3))
    
    # Identify combined events
    combined_events <- wls %>%
      mutate(check = .fitted + max(2*wls_se, 1)) %>%
      mutate(combo_event = if_else(Q_Runoff_mm > check, 1, 0))
    
    # # Make a plot
    # plot <- combined_events %>%
    #   ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm, color = combo_event)) +
    #   geom_point() +
    #   geom_abline(slope = wls_slope, intercept = wls_yint) +
    #   labs(x = "Precipitation Depth (mm)", y = "Runoff Depth (mm)", 
    #        title = i, color = "") +
    #   theme_bw() + theme(legend.position = "none")
    # print(plot)
    
    # Find number of combined events
    num_combined_events = sum(combined_events$combo_event)
    
    # If there are still combined events, remove and repeat
    if (num_combined_events > 0) {
      
      events <- combined_events %>%
        filter(combo_event == 0) %>%
        select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
      rm(ols, ols_resid, weights, wls, wls_slope, wls_yint, wls_se, 
         combined_events)
      
    } else {
      
      # Save results
      wls_results_i <- 
        tibble(Site = i, n = nrow(wls), yint = wls_yint, slope = wls_slope)
      wls_results <- bind_rows(wls_results, wls_results_i)
      eia_events_i <- combined_events %>%
        select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
      eia_events <- bind_rows(eia_events, eia_events_i)
      
    }
  }
}
rm(i, events, num_combined_events, ols, ols_resid, 
   weights, wls, wls_slope, wls_yint, wls_se, 
   combined_events, wls_results_i, eia_events_i)



### Step 04: Plot Results ------------------------------------------------------

plot_list <- list()  #Another modification from MJB
get_all_plotting_data <- NULL

for (i in unique(data$Site)) {
  
  if (event_type == "all") {
    
    # Classify events
    eia_events_i <- eia_events %>% 
      filter(Site == i) %>%
      mutate(combo_event = "EIA Event")
    non_eia_events_i <- events_all %>% 
      filter(Site == i) %>%
      anti_join(eia_events_i, by = c("Site", "Precip_EndTime", 
                                     "Precip_Total_mm", "Q_Runoff_mm")) %>%
      mutate(combo_event = "Combo Event")
    events_combo <- bind_rows(eia_events_i, non_eia_events_i)
    rm(eia_events_i, non_eia_events_i)
    
  } else {
    
    # Classify events
    eia_events_i <- eia_events %>% 
      filter(Site == i) %>%
      mutate(combo_event = "EIA Event")
    non_eia_events_i <- events_nonzero %>% 
      filter(Site == i) %>%
      anti_join(eia_events_i, by = c("Site", "Precip_EndTime", 
                                     "Precip_Total_mm", "Q_Runoff_mm")) %>%
      mutate(combo_event = "Combo Event")
    events_combo <- bind_rows(eia_events_i, non_eia_events_i)
    #Print maybe?
    print(i)
    
    rm(eia_events_i, non_eia_events_i)
    
  }
  
  # Pull weighted linear regression coefficiencts
  wls_results_i <- wls_results %>%
    filter(Site == i)
  
  # Make plot
  plot <- events_combo %>%
    ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm, color = combo_event)) +
    geom_point() +
    geom_abline(slope = wls_results_i$slope, intercept = wls_results_i$yint) +
    labs(x = "Precipitation Depth (mm)", y = "Runoff Depth (mm)", 
         title = i, color = "") +
    theme_bw()
 # print(plot)
  plot_list[[i]] <- plot #MJB to plot them all
  get_all_plotting_data <- rbind(get_all_plotting_data, events_combo)
  
}

rm(i, events_combo, wls_results_i, plot)

#Plotting
grid.arrange(grobs=plot_list,ncol=3)

head(wls_results)


head(flow_data_to_sample_screened)
#Plot each site over time, along with above?
#Further screen the original data.frame based on the date start
flow_data_to_sample_screened_start <- flow_data_to_sample_screened[which(flow_data_to_sample_screened$RainEnd_UTC >= date_start),]


#Main label d.f to match
name_match <- data.frame(Site = c("Br", "D4", "Fe", "L4", "L1", "Ln", "Ls", "Sa"), main = c("Br (control)", "D4 (impact)",
                                                                                            "Fe (control)", "L4 (impact)", "L1 (impact)",
                                                                                            "Ln (impact)", "Ls (impact)", "Sa (reference)"))


pdf("~/Documents/Git/lsc_flow/DCI/dci.pdf",width=10,height=10)
par(mfrow=c(3,3))

for (i in unique(flow_data_to_sample_screened_start$new_sites_name)) {
  i_pos <- which(flow_data_to_sample_screened_start$new_sites_name == i) 
  
  i_slope <- wls_results$slope[which(wls_results$Site == i)]
  title <- name_match$main[which(name_match$Site == i)]
  
  
  #Try plotting
  plot(flow_data_to_sample_screened_start$RainEnd_UTC[i_pos], flow_data_to_sample_screened_start$treatments_ei[i_pos], type = 'l', main = title, ylim = c(0, (1.1*max(i_slope, (flow_data_to_sample_screened_start$treatments_ei[i_pos])))),
       xlab = "Datetime", ylab = "EIA")
  
  lines(flow_data_to_sample_screened_start$RainEnd_UTC[i_pos], flow_data_to_sample_screened_start$treatments_ro[i_pos], col = 'blue', lty = 2)
  
  abline(h = i_slope, col = 'red')
  
  print(paste(i, "-", mean(flow_data_to_sample_screened_start$treatments_ei[i_pos]), "-", mean(flow_data_to_sample_screened_start$treatments_ro[i_pos])))
}

plot(0, type = 'n', axes = FALSE, ann = FALSE)
legend("bottomright", c("Mapped EIA", "Mapped EIA + SCMs", "EIA_Ebrahimian"),
       col = c('black', 'blue', 'red'), lty = c(1,2,1), cex = 1, xjust = 1, yjust = 1, bty="n")

dev.off()

#get_all_plotting_data
#Get the N, total and EIA events
str(get_all_plotting_data)
sites_n <- levels(factor(get_all_plotting_data$Site))

number_events <- NULL
for(i in 1:NROW(sites_n)) {
  site_i <- get_all_plotting_data[which(get_all_plotting_data$Site == sites_n[i]),]
  total_n_i <- NROW(site_i)
  eia_n_i <- NROW(which(site_i$combo_event == "EIA Event"))
  combine <- data.frame(site = sites_n[i], N = total_n_i, eia_N = eia_n_i)
  number_events <- rbind(number_events, combine) }

number_events$per <- (number_events$eia_N/number_events$N)*100
number_events$per_round <- round(((number_events$eia_N/number_events$N)*100), 0)

#Percentage
wls_results$ei_estimate <- round((wls_results$slope * 100),0)

#Initial abstraction
#Set y = 0
wls_results$int <- round((wls_results$yint * -1)/wls_results$slope, 1)




