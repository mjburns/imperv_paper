## Estimate Effective Impervious Area Using Rainfall-Runoff Data
## Based on Ebrahimian and others (2016) Weighted Linear Regression

# Charlie Stillwell
# February 27, 2025
 


### Step 00: Set Up ------------------------------------------------------------

# Load packages
library(tidyverse)
library(broom)
library(cowplot)

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
data <- get(load("Input/lsc_events_no_precip_max.Rdata"))
data_aus <- data %>%
  mutate(region = "Melbourne") %>%
  relocate(region, .before = Site)
rm(data)

data_usa <- 
  read_csv("Input/Clarksburg_StormEvents_2004_2018.csv", 
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
  relocate(region, .before = Site) %>%
  filter(Site != "Forested Control")

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
wls_results_table2 <- tibble()
wls_results_appendix <- tibble()
eia_events <- tibble()

for (i in unique(data$Site)) {
  
  # Select event types to include in analysis
  if (event_type == "all") {
    events <- events_all %>% filter(Site == i)
  } else {
    events <- events_nonzero %>% filter(Site == i)
  }
  
  events_initial <- events
  
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
    wls_xint <- (-1*wls_yint)/wls_slope
    wls_ci <- confint(wls)
    wls_adjrsquared <- summary(wls)$adj.r.squared
    wls_mse <- mean(summary(wls)$residuals^2)
    wls_se_at_xint <- predict(wls, 
                              newdata = tibble(Precip_Total_mm = wls_xint), 
                              se.fit = TRUE)$se.fit
    wls_xint_lb <- (-wls_se_at_xint - wls_yint) / wls_slope
    wls_xint_ub <- (wls_se_at_xint - wls_yint) / wls_slope
    wls <- augment(wls, data = weights)
    wls_se <- sqrt(sum(wls$Weight * wls$.resid^2) / (nrow(wls) - 2))
    
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
         wls_slope, wls_yint, wls_xint, wls_se_at_xint, wls_xint_lb, wls_xint_ub, 
         wls_se, combined_events)
      
    } else {
      
      # Save results
      wls_results_table2_i <- 
        tibble(Site = i, 
               total_events = nrow(events_initial), 
               events_per_year = 
                 if_else(i == "Urban Control" | i == "Treatment 1" | i == "Treatment 2", 
                         round(nrow(events_initial) / 2.75, digits = 1),
                         round(nrow(events_initial) / 3.75, digits = 1)),
               ei_events = nrow(events), 
               ei_events_per_year = 
                 if_else(i == "Urban Control" | i == "Treatment 1" | i == "Treatment 2", 
                         round(nrow(events) / 2.75, digits = 1), 
                         round(nrow(events) / 3.75, digits = 1)),
               pct_ei_events = round(nrow(events)/nrow(events_initial), digits = 3),
               slope = wls_slope, 
               init_abstraction = wls_xint)
      wls_results_appendix_i <- 
        tibble(Site = i, 
               adj_r_squared = wls_adjrsquared, 
               mean_squared_error = wls_mse,
               yint = wls_yint, 
               yint_lb = wls_ci[[1]], 
               yint_ub = wls_ci[[3]], 
               slope = wls_slope, 
               slope_lb = wls_ci[[2]], 
               slope_ub = wls_ci[[4]], 
               xint = wls_xint, 
               xint_lb = wls_xint_lb, 
               xint_ub = wls_xint_ub,
               wls_se = wls_se)
      wls_results_table2 <- bind_rows(wls_results_table2, wls_results_table2_i)
      wls_results_appendix <- bind_rows(wls_results_appendix, wls_results_appendix_i)
      eia_events_i <- combined_events %>%
        select(Site, Precip_EndTime, Precip_Total_mm, Q_Runoff_mm)
      eia_events <- bind_rows(eia_events, eia_events_i)
      
    }
  }
}
rm(i, events, num_combined_events, ols, ols_resid, 
   weights, wls, wls_slope, wls_yint, wls_xint, wls_se_at_xint, wls_xint_lb, wls_xint_ub,
   wls_se, wls_ci, wls_adjrsquared, wls_mse,
   combined_events, wls_results_table2_i, wls_results_appendix_i, eia_events_i)



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
#   # Pull weighted linear regression coefficients
#   wls_results_i <- wls_results_appendix %>%
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



### Step 05: Make Plots for Manuscript -----------------------------------------

# Specify if event is EIA or combo
combo_events <- 
  anti_join(events_all, eia_events, 
            by = c("Site", "Precip_EndTime", 
                   "Precip_Total_mm", "Q_Runoff_mm")) %>%
  mutate(Event_Type = "Combo")
eia_events <- mutate(eia_events, Event_Type = "EI")
events_all <- bind_rows(eia_events, combo_events)
rm(combo_events, eia_events)

# Prepare data for plots
events_all_plots <- events_all %>%
  filter(Site != "Urban Control") %>%
  mutate(Event_Type = case_when(Event_Type == "Combo" ~ "Combination", 
                                Event_Type == "EI" ~ "Effective Impervious")) %>%
  mutate(Site = case_when(Site == "D4" ~ "Melbourne - D4", 
                          Site == "L4" ~ "Melbourne - L4", 
                          Site == "L1" ~ "Melbourne - L1", 
                          Site == "Ln" ~ "Melbourne - Ln", 
                          Site == "Ls" ~ "Melbourne - Ls", 
                          Site == "Treatment 1" ~ "Clarksburg - Treatment 1", 
                          Site == "Treatment 2" ~ "Clarksburg - Treatment 2"))

# Figure 2
figure2 <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Melbourne - Ls" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Melbourne - Ls" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, color = "black") +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", color = "Event Type") +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title.x = 
          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), 
                       size = 10, color = "black"), 
        axis.title.y = 
          element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                       size = 10, color = "black"), 
        axis.text = element_text(size = 9, color = "black"))
figure2_legend <- events_all_plots %>%
  filter(Site == "Melbourne - Ls") %>%
  ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm, color = Event_Type)) +
  geom_point(shape = 1) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", color = "Event Type") +
  theme_bw() +
  theme(legend.position = "right", 
        legend.title = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 9, color = "black"))
figure2_legend <- get_legend(figure2_legend) %>% suppressWarnings()
figure2 <- plot_grid(figure2, figure2_legend, rel_heights = c(3,1),
                      align = "hv", axis = "tblr", nrow = 2)
figure2
ggsave("figure2.png", figure2, width = 3, height = 4, units = "in")

# Figure 3
plot_D4 <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - D4" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - D4" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - D4" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - D4" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Melbourne - D4" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Melbourne - D4" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Melbourne - D4") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(), 
        axis.title.y = 
          element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                       size = 10, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_L4 <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - L4" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - L4" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - L4" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - L4" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Melbourne - L4" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Melbourne - L4" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Melbourne - L4") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_L1 <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - L1" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - L1" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - L1" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - L1" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Melbourne - L1" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Melbourne - L1" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Melbourne - L1") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(), 
        axis.title.y = 
          element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                       size = 10, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_Ln <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - Ln" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - Ln" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - Ln" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - Ln" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Melbourne - Ln" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Melbourne - Ln" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Melbourne - Ln") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(), 
        axis.title.y = 
          element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                       size = 10, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_Ls <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Melbourne - Ls" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Melbourne - Ls" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Melbourne - Ls" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Melbourne - Ls") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_Treat1 <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 1" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 1" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 1" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 1" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Clarksburg - Treatment 1" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Clarksburg - Treatment 1" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Clarksburg - Treatment 1") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = 
          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), 
                       size = 10, color = "black"),
        axis.title.y = 
          element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), 
                       size = 10, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_Treat2 <- ggplot() +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 2" & 
                              Event_Type == "Effective Impervious")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 2" & 
                              Event_Type == "Effective Impervious")$Q_Runoff_mm),
             shape = 1, color = "#00BFC4") +
  geom_point(aes(x = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 2" & 
                              Event_Type == "Combination")$Precip_Total_mm,
                 y = filter(events_all_plots, 
                            Site == "Clarksburg - Treatment 2" & 
                              Event_Type == "Combination")$Q_Runoff_mm),
             shape = 1, color = "#F8766D") +
  geom_smooth(aes(x = filter(events_all_plots,
                             Site == "Clarksburg - Treatment 2" &
                               Event_Type == "Effective Impervious")$Precip_Total_mm,
                  y = filter(events_all_plots,
                             Site == "Clarksburg - Treatment 2" &
                               Event_Type == "Effective Impervious")$Q_Runoff_mm),
              method = "lm", formula = y ~ x, se = TRUE, 
              color = "black", linewidth = 0.5) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", 
       title = "Clarksburg - Treatment 2") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.title.x = 
          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), 
                       size = 10, color = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_legend <- events_all_plots %>%
  filter(Site == "Melbourne - D4") %>%
  ggplot(aes(x = Precip_Total_mm, y = Q_Runoff_mm, color = Event_Type)) +
  geom_point(shape = 1) +
  labs(x = "Precipitation (mm)", y = "Quickflow (mm)", color = "Event Type") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 10)) +
  theme_bw() +
  theme(legend.position = "right", 
        legend.title = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 9, color = "black"), 
        plot.title = element_text(size = 11),
        plot.margin = unit(c(0, 0, 0, 0), "in"))
plot_legend <- get_legend(plot_legend) %>% suppressWarnings()
figure3 <- plot_grid(plot_D4, plot_legend, 
                     plot_L1, plot_L4, 
                     plot_Ln, plot_Ls, 
                     plot_Treat1, plot_Treat2, 
                     align = "hv", axis = "tblr", nrow = 4, ncol = 2)
figure3
ggsave("figure3.png", figure3, height = 8, width = 6.5, units = "in")

# Figure 4
fig4_data <- tibble(site = c("Melbourne - D4", "Melbourne - L4", 
                             "Melbourne - L1", "Melbourne - Ln", 
                             "Melbourne - Ls", "Clarksburg - Treatment 1",
                             "Clarksburg - Treatment 2"), 
                    method_1 = c(1.8, 6.1, 21.3, 1.6, 6.4, 7.7, 12.2), 
                    method_2 = c(2.8, 7.7, 12.1, 4.1, 12.2, 5.5, 9.9), 
                    total_imp = c(6.7, 14.1, 27.8, 8.3, 22.3, 33, 44))
figure4_plot <- fig4_data %>%
  pivot_longer(cols = c(method_1, method_2, total_imp), 
               names_to = "Method", values_to = "Imperviousness") %>%
  mutate(Method = case_when(Method == "method_1" ~ "EI, Flow Disturbance Frequency", 
                            Method == "method_2" ~ "EI, Rainfall-Runoff Regression", 
                            Method == "total_imp" ~ "Total Imperviousness")) %>%
  ggplot(aes(x = Imperviousness/100, y = reorder(site, Imperviousness/100, FUN = max))) +
  geom_line() +
  geom_point(aes(color = Method), size = 2) +
  scale_x_continuous(limits = c(0, 0.5), labels = scales::percent) +
  labs(x = "Imperviousness", y = "", color = "") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), 
                                    size = 10, color = "black"), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 9, color = "black"))
figure4_legend <- fig4_data %>%
  pivot_longer(cols = c(method_1, method_2, total_imp), 
               names_to = "Method", values_to = "Imperviousness") %>%
  mutate(Method = case_when(Method == "method_1" ~ "EI, Flow Disturbance Frequency", 
                            Method == "method_2" ~ "EI, Rainfall-Runoff Regression", 
                            Method == "total_imp" ~ "Total Imperviousness")) %>%
  ggplot(aes(x = Imperviousness/100, y = reorder(site, Imperviousness/100, FUN = max))) +
  geom_line() +
  geom_point(aes(color = Method), size = 2) +
  scale_x_continuous(limits = c(0, 0.5), labels = scales::percent) +
  labs(x = "Imperviousness", y = "", color = "") +
  theme_bw() +
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 9, color = "black"))
figure4_legend <- get_legend(figure4_legend) %>% suppressWarnings()
figure4 <- plot_grid(figure4_plot, figure4_legend, 
                     align = "v", axis = "r", nrow = 1, 
                     rel_widths = c(2,1))
figure4
ggsave("figure4.png", figure4, height = 3, width = 6.5, units = "in")
