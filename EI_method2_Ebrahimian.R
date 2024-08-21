## Estimate Effective Impervious Area Using Rainfall-Runoff Data
## Based on Ebrahimian and others (2016) Weighted Linear Regression

# Charlie Stillwell
# August 21, 2024
 


### Step 00: Set Up ------------------------------------------------------------

# Load packages
library(tidyverse)
library(broom)

# Specify max and min precip depths to include in analysis (match data units)
precip_min <- 2.54
precip_max <- 254000

# Choose types of events (only runoff-generating events or all precip events)
event_type <- "all"     # one of "all" or "nonzero"

# Specify date range of events to include in analysis
# (for my data, use precipitation end date-time)
date_start_usa <- as.Date("2016-01-01")
date_end_usa <- as.Date("2018-12-31")

# Read in rainfall-runoff data
load("input/lsc_events.RData")
data_aus <- data %>%
  mutate(region = "Melbourne") %>%
  relocate(region, .before = Site)
rm(data)

data_usa <- 
  read_csv("input/Clarksburg_StormEvents_2004_2018.csv", 
           col_types = cols(Precip_Duration_hrs = col_skip(), 
                            Precip_EndTime = 
                              col_datetime(format = "%m/%d/%Y %H:%M"), 
                            Precip_PriorDry_hrs = col_skip(), 
                            Precip_StartTime = col_skip(), 
                            Q_Duration_hrs = col_skip(), 
                            Q_EndTime = col_skip(), 
                            Q_PeakTime = col_skip(), 
                            Q_Peak_cfs = col_skip(), 
                            Q_Peak_cms.sqkm = col_skip(), 
                            Q_PriorDry_hrs = col_skip(), 
                            Q_RiseRate_cms.hr = col_skip(),  
                            Q_StartTime = col_skip(), 
                            Q_TimeToFlowPeak_hrs = col_skip(),
                            Q_Total_mm = col_skip(),
                            Q_Total_PrecipRatio = col_skip())) %>%
  filter(Precip_EndTime >= date_start_usa & Precip_EndTime <= date_end_usa) %>%
  filter(Precip_Total_mm >= precip_min & Precip_Total_mm <= precip_max) %>%
  select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm) %>%
  mutate(region = "Clarksburg") %>%
  relocate(region, .before = Site)

# Combine datasets
data <- bind_rows(data_aus, data_usa)
rm(data_aus, data_usa)



### Step 01: Plot Rainfall-Runoff Data -----------------------------------------

# for (i in unique(data$Site)) {
#   plot <- data %>%
#     filter(Site == i) %>%
#     ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm)) +
#     geom_abline(slope = 1, intercept = 0) +
#     geom_point() +
#     scale_x_continuous(limits = c(0, max(data$Precip_Total_mm))) +
#     scale_y_continuous(limits = c(0, max(data$Precip_Total_mm))) +
#     labs(title = i) +
#     theme_bw()
#   print(plot)
# }
# rm(i, plot)



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
    # wls_xint <- (-1*wls_yint)/wls_slope
    wls_ci <- confint(wls)
    wls_adjrsquared <- summary(wls)$adj.r.squared
    wls_mse <- mean(summary(wls)$residuals^2)
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
      rm(ols, ols_resid, weights, wls, wls_ci, wls_adjrsquared, wls_mse,
         wls_slope, wls_yint, wls_se, combined_events)
      
    } else {
      
      # Save results
      wls_results_i <- 
        tibble(Site = i, n = nrow(wls), 
               adj_r_squared = wls_adjrsquared, mean_squared_error = wls_mse,
               yint = wls_yint, yint_lb = wls_ci[[1]], yint_ub = wls_ci[[3]], 
               slope = wls_slope, slope_lb = wls_ci[[2]], slope_ub = wls_ci[[4]])
      wls_results <- bind_rows(wls_results, wls_results_i)
      eia_events_i <- combined_events %>%
        select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
      eia_events <- bind_rows(eia_events, eia_events_i)
      
    }
  }
}
rm(i, events, num_combined_events, ols, ols_resid, 
   weights, wls, wls_slope, wls_yint, wls_se, wls_ci, wls_adjrsquared, wls_mse,
   combined_events, wls_results_i, eia_events_i)



### Step 04: Plot Results ------------------------------------------------------

# for (i in unique(data$Site)) {
#   
#   if (event_type == "all") {
#     
#     # Classify events
#     eia_events_i <- eia_events %>% 
#       filter(Site == i) %>%
#       mutate(combo_event = "EIA Event")
#     non_eia_events_i <- events_all %>% 
#       filter(Site == i) %>%
#       anti_join(eia_events_i, by = c("Site", "Precip_EndTime", 
#                                      "Precip_Total_mm", "Q_Runoff_mm")) %>%
#       mutate(combo_event = "Combo Event")
#     events_combo <- bind_rows(eia_events_i, non_eia_events_i)
#     rm(eia_events_i, non_eia_events_i)
#     
#   } else {
#     
#     # Classify events
#     eia_events_i <- eia_events %>% 
#       filter(Site == i) %>%
#       mutate(combo_event = "EIA Event")
#     non_eia_events_i <- events_nonzero %>% 
#       filter(Site == i) %>%
#       anti_join(eia_events_i, by = c("Site", "Precip_EndTime", 
#                                      "Precip_Total_mm", "Q_Runoff_mm")) %>%
#       mutate(combo_event = "Combo Event")
#     events_combo <- bind_rows(eia_events_i, non_eia_events_i)
#     rm(eia_events_i, non_eia_events_i)
#     
#   }
#   
#   # Pull weighted linear regression coefficiencts
#   wls_results_i <- wls_results %>%
#     filter(Site == i)
#   
#   # Make plot
#   plot <- events_combo %>%
#     ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm, color = combo_event)) +
#     geom_point() +
#     geom_abline(slope = wls_results_i$slope, intercept = wls_results_i$yint) +
#     labs(x = "Precipitation Depth (mm)", y = "Runoff Depth (mm)", 
#          title = i, color = "") +
#     theme_bw()
#   print(plot)
#   
# }
# rm(i, events_combo, wls_results_i, plot)

