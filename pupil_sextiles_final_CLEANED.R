# ============================================================================
# Pupil Sextile Analysis: Ripple Rate Modulation by Pupil Size
# ============================================================================
# Analysis of ripple activity across pupil size sextiles in eyes-open periods
# Including EOEC, main pupil analyses, blink analyses, and modulation indices
# Last updated: April 2026, June 2026
# ============================================================================

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(lme4)
library(car)
library(lmerTest)
library(ggResidpanel)
library(ggeffects)
library(emmeans)
library(gridExtra)
library(tidyr)


# ============================================================================
# 2. DATA LOADING AND ORGANIZATION
# ============================================================================

# Load main pupil array with all physiological data
pupil_array <- read_csv("data/pupil_array_1_3.csv", col_names = FALSE)
colnames(pupil_array) <- c(
  "sub_num", "emublock", "matblock", "pup_time", 
  "raw_pup", "raw_sextile", "trial_num", 
  "EOEC", "trial_reject", "ripple", "ripple_count", "amy_count", "ant_hipp_count", "post_hipp_count",
  "electrode_num",
  "hipp_or_amy", "A_or_P", "keep_or_rej", "duration", "amplitude", "freq",
  "ripple_time", "ripple_idx", "z_pup", "z_sextile", "blink", "eye")

# Remove subject 56 (only amygdala contacts that cannot be reliably kept)
pupil_array <- pupil_array %>%
  filter(sub_num != 56)

# Remove hipp_or_amy from dataframe (will be added from electrode position data)
pupil_array <- pupil_array %>% select(-hipp_or_amy)

# Filter array for the eye to use for each subject (one eye per subject)
pupil_array_filter <- pupil_array %>%
  filter(case_when(
    sub_num %in% c(38, 48, 53, 54, 58, 69) ~ eye == 0,
    TRUE ~ eye == 1
  ))

# Filter out rejected ripples
pupil_array_filter2 <- pupil_array_filter %>%
  mutate(ripple = ifelse(keep_or_rej == 1, ripple, 0),
         ripple_count = ifelse(keep_or_rej == 1, ripple_count, 0))

# Load electrode position information
sub_elec_positions <- read_csv("data/sub_elec_positions_allinfo_2025.csv", col_names = FALSE)
colnames(sub_elec_positions) <- c(
  "sub_num", "electrode_num", "left_or_right", "hipp_or_amy", "AP", "RL", "SI", "AP_uncal", "apex_coord",
  "elec_coord", "coord_diff", "close_flag", "MNI_coord", "post_50", "post_35", "post_MNI",
  "overlap_50", "overlap_MNI", "overlap_35"
)

# Load sleep electrode information with modulation scores
sleep_elecs <- read_csv("data/sleep_electrode_array_updated_final_9242025.csv", col_names = FALSE)
colnames(sleep_elecs) <- c(
  "sub_num", "electrode_num", "keep_or_reject", "hipp_or_amy", "sleep_mod_score", 
  "sleep_amp_score", "sleep_dur_score", "sleep_freq_score"
)

# Join electrode position information to pupil array
pupil_array_final <- pupil_array_filter2 %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, hipp_or_amy, left_or_right, AP, AP_uncal), 
            by = c("sub_num", "electrode_num"))

# Join sleep electrode information (including modulation scores)
pupil_array_final <- pupil_array_final %>%
  left_join(sleep_elecs %>%
              select(sub_num, electrode_num, sleep_mod_score, keep_or_reject), 
            by = c("sub_num", "electrode_num"))

# Select only good electrodes (keep_or_reject == 0 or NA)
pupil_array_final <- filter(pupil_array_final, keep_or_reject == 0 | is.na(keep_or_reject))

# Remove specific electrodes that were removed from sleep analysis
pupil_array_final <- pupil_array_final %>%
  filter(!(sub_num == 49 & electrode_num == 38)) %>%
  filter(!(sub_num == 62 & electrode_num == 134))

# Check for artifactual ripples (amplitude >= 300 microvolts)
bad_ripple <- pupil_array_final %>% 
  filter(amplitude >= 300)

# Create dataset subsets
# Main clean dataset: all good trials
pupil_array_final_clean <- pupil_array_final %>%
  filter(trial_reject == 0)

# Eyes-open subset: EOEC == 1
pupil_array_EO <- pupil_array_final %>%
  filter(EOEC == 1)

# Eyes-open clean subset: EOEC == 1 and good trials
pupil_array_EO_clean <- pupil_array_EO %>%
  filter(trial_reject == 0)

# Summary information: unique emublocks and trials per subject
unique_emublocks <- pupil_array_EO_clean %>%
  group_by(sub_num) %>%
  summarise(n_unique_emublocks = n_distinct(emublock), .groups = "drop")

unique_emublocks_trials <- pupil_array_EO_clean %>%
  group_by(sub_num, emublock) %>%
  summarise(n_unique_trials = n_distinct(trial_num), .groups = "drop") %>%
  group_by(sub_num) %>%
  summarise(n_trials = sum(n_unique_trials), .groups = "drop")

# Summary statistics
median(unique_emublocks$n_unique_emublocks)
quantile(unique_emublocks$n_unique_emublocks, 0.25)  # Q1
quantile(unique_emublocks$n_unique_emublocks, 0.75)  # Q3

mean_rows_per_sub_num <- pupil_array_EO_clean %>%
  group_by(sub_num) %>%
  summarise(n_rows = n(), .groups = "drop") %>%
  mutate(mean_duration_sec = n_rows / 200) %>%
  mutate(mean_duration_min = mean_duration_sec / 60)

mean(mean_rows_per_sub_num$mean_duration_sec) / 60
sd(mean_rows_per_sub_num$mean_duration_sec) / 60

# ============================================================================
# 3. EOEC ANALYSES: Eyes-Open vs Eyes-Closed Comparison
# ============================================================================

# 3.1 Data Organization for EOEC
# ----

# Time points per EOEC condition
ripple_sextile_timesEOEC <- pupil_array_final_clean %>%
  group_by(sub_num, EOEC) %>%
  summarise(total_timepoints = n())

# Ripple counts per electrode and EOEC condition
ripple_sextile_countsEOEC <- pupil_array_final_clean %>%
  group_by(sub_num, EOEC, electrode_num) %>%
  summarise(
    total_events = sum(ifelse(hipp_or_amy == 2, amy_count, 
                             ifelse(AP_uncal == 0, ant_hipp_count, post_hipp_count)))) %>%
  ungroup() %>%
  group_by(sub_num) %>%
  complete(EOEC = 1:2,
           electrode_num = unique(electrode_num), 
           fill = list(total_events = 0)) %>%
  ungroup()

# Remove electrode 0 and any NaN rows
ripple_sextile_countsEOEC <- ripple_sextile_countsEOEC %>% filter(electrode_num != 0)
ripple_sextile_countsEOEC <- ripple_sextile_countsEOEC[complete.cases(ripple_sextile_countsEOEC), ]
ripple_sextile_timesEOEC <- ripple_sextile_timesEOEC[complete.cases(ripple_sextile_timesEOEC), ]

# Add region and modulation information
ripple_sextile_countsEOEC <- ripple_sextile_countsEOEC %>%
  left_join(sleep_elecs %>%
              select(sub_num, electrode_num, hipp_or_amy, sleep_mod_score), 
            by = c("sub_num", "electrode_num")) %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, AP_uncal),
            by = c("sub_num", "electrode_num"))

# Merge times and counts
ripple_rate_resultEOEC <- merge(ripple_sextile_timesEOEC, ripple_sextile_countsEOEC, 
                                 by = c("sub_num", "EOEC"), all = TRUE)
ripple_rate_resultEOEC <- ripple_rate_resultEOEC %>%
  select(sub_num, EOEC, total_timepoints, electrode_num, hipp_or_amy, total_events, 
         sleep_mod_score, AP_uncal)

# Calculate ripple rates (pupil sampled at 200Hz, so multiply timepoints by 5ms)
ripple_rate_finalEOEC <- ripple_rate_resultEOEC %>%
  mutate(ripple_rate = total_events / ((total_timepoints * 5) / 1000))

# Get long format data
long_data_rateEOEC <- ripple_rate_finalEOEC %>%
  select(sub_num, EOEC, electrode_num, hipp_or_amy, AP_uncal, total_events, ripple_rate) %>%
  filter(!is.nan(EOEC))

# Relabel subject numbers (42 -> 49 for consistency)
long_data_rateEOEC$sub_num[long_data_rateEOEC$sub_num == 42] <- 49

# Calculate subject and aggregate averages
subject_avg_dataEOEC <- long_data_rateEOEC %>%
  group_by(sub_num, EOEC, hipp_or_amy) %>%
  summarise(ripple_rate = mean(ripple_rate), .groups = 'drop')

agg_data_simpleEOEC <- subject_avg_dataEOEC %>%
  group_by(EOEC, hipp_or_amy) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Separate for hippocampus and amygdala
hipp_aggEOEC <- agg_data_simpleEOEC %>% filter(hipp_or_amy == 1)
hipp_subjectEOEC <- subject_avg_dataEOEC %>% filter(hipp_or_amy == 1)

amy_aggEOEC <- agg_data_simpleEOEC %>% filter(hipp_or_amy == 2)
amy_subjectEOEC <- subject_avg_dataEOEC %>% filter(hipp_or_amy == 2)

# Create A/P/Amy variable
long_data_rate_APEOEC <- long_data_rateEOEC %>%
  mutate(A_or_P_or_amy = ifelse(hipp_or_amy == 2, 3, 
                                ifelse(AP_uncal == 0, 0, 1)))

long_data_rate_hippEOEC <- long_data_rate_APEOEC %>% filter(hipp_or_amy == 1) 

subject_avg_data_hippEOEC <- long_data_rate_hippEOEC %>%
  group_by(sub_num, EOEC, AP_uncal) %>%
  summarise(ripple_rate = mean(ripple_rate), .groups = 'drop')

agg_data_APEOEC <- subject_avg_data_hippEOEC %>% 
  group_by(EOEC, AP_uncal) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

subject_avg_data_hipp_AEOEC <- subject_avg_data_hippEOEC %>% filter(AP_uncal == 0)
subject_avg_data_hipp_PEOEC <- subject_avg_data_hippEOEC %>% filter(AP_uncal == 1)

# 3.2 EOEC Plotting Function
# ----

plot_EOEC_summary <- function(summary_df, EOEC_data, value_col, mean_col, sem_col, colors, 
                               y_label, y_limits) {
  # Convert inputs to symbols for tidy evaluation
  value_sym <- sym(value_col)
  mean_sym <- sym(mean_col)
  sem_sym <- sym(sem_col)
  
  # Subject-level summary
  subject_avg <- EOEC_data %>%
    group_by(sub_num, EOEC) %>%
    summarise(value = mean(!!value_sym, na.rm = TRUE), .groups = "drop")
  
  ggplot() +
    geom_line(data = subject_avg,
              aes(x = factor(EOEC), y = value, group = sub_num),
              color = "lightgrey", alpha = 0.6) +
    geom_point(data = subject_avg,
               aes(x = factor(EOEC), y = value, group = sub_num, color = factor(EOEC)),
               size = 3, alpha = 0.2) +
    geom_line(data = summary_df,
              aes(x = factor(EOEC), y = !!mean_sym, group = 1),
              color = "black", size = 1) +
    geom_errorbar(data = summary_df,
                  aes(x = factor(EOEC), ymin = !!mean_sym - !!sem_sym, ymax = !!mean_sym + !!sem_sym),
                  width = 0.2, size = 1, color = 'black') +
    geom_point(data = summary_df,
               aes(x = factor(EOEC), y = !!mean_sym, fill = factor(EOEC)),
               shape = 21, size = 4, stroke = 1.5) +
    scale_x_discrete(limits = c("2", "1"), labels = c("1" = "Eyes-Open", "2" = "Eyes-Closed")) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(y = y_label) +
    theme_minimal() +
    ylim(y_limits) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) 
}

# 3.3 EOEC Plot and Statistics
# ----

rate2 <- plot_EOEC_summary(
  summary_df = hipp_aggEOEC,
  EOEC_data = hipp_subjectEOEC,
  value_col = "ripple_rate",
  mean_col = "ripple_rate_mean",
  sem_col = "ripple_rate_sem",
  colors = c("1" = "#6E8E6E", "2" = "#C5D3C6"),
  y_label = "Mean Rate",
  y_limits = c(0.01, .07)
)
rate2

# Mixed effects model: EOEC effect on hippocampal ripple rate
long_data_rateEOEC_hipp <- long_data_rateEOEC %>% filter(hipp_or_amy == 1)
model.pup1 <- lmer(ripple_rate ~ EOEC + (1 | sub_num/electrode_num), data = long_data_rateEOEC_hipp)
summary(model.pup1)
car::Anova(model.pup1, type = 'III')

pp <- plot(ggpredict(model.pup1, terms = c('EOEC'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp


# ============================================================================
# 4. MAIN PUPIL ANALYSES: Ripple Rate by Pupil Size Sextiles
# ============================================================================

# 4.1 Data Organization for Pupil Sextiles
# ----

# Time points per pupil sextile
ripple_sextile_times <- pupil_array_EO_clean %>%
  group_by(sub_num, z_sextile) %>%
  summarise(total_timepoints = n())

# Ripple counts per electrode and pupil sextile
ripple_sextile_counts <- pupil_array_EO_clean %>%
  group_by(sub_num, z_sextile, electrode_num) %>%
  summarise(
    total_events = sum(ifelse(hipp_or_amy == 2, amy_count, 
                             ifelse(AP_uncal == 0, ant_hipp_count, post_hipp_count)))) %>%
  ungroup() %>%
  group_by(sub_num) %>%
  complete(z_sextile = 1:6,
           electrode_num = unique(electrode_num), 
           fill = list(total_events = 0)) %>%
  ungroup()

# Remove electrode 0 and any NaN rows
ripple_sextile_counts <- ripple_sextile_counts %>% filter(electrode_num != 0)
ripple_sextile_counts <- ripple_sextile_counts[complete.cases(ripple_sextile_counts), ]

# Add region and modulation information
ripple_sextile_counts <- ripple_sextile_counts %>%
  left_join(sleep_elecs %>%
              select(sub_num, electrode_num, hipp_or_amy, sleep_mod_score), 
            by = c("sub_num", "electrode_num")) %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, AP_uncal),
            by = c("sub_num", "electrode_num"))

# Merge times and counts
ripple_rate_result <- merge(ripple_sextile_times, ripple_sextile_counts, 
                            by = c("sub_num", "z_sextile"), all = TRUE)
ripple_rate_result <- ripple_rate_result %>%
  select(sub_num, z_sextile, total_timepoints, electrode_num, hipp_or_amy, total_events, 
         sleep_mod_score, AP_uncal)

# Calculate ripple rates
ripple_rate_final <- ripple_rate_result %>%
  mutate(ripple_rate = total_events / ((total_timepoints * 5) / 1000))

# Get long format data
long_data_rate <- ripple_rate_final %>%
  select(sub_num, z_sextile, electrode_num, hipp_or_amy, AP_uncal, total_events, ripple_rate) %>%
  filter(!is.nan(z_sextile))

# Create modulation index for each electrode
modulation_df_elec <- long_data_rate %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, AP),
            by = c("sub_num", "electrode_num")) %>%
  group_by(sub_num, electrode_num, z_sextile) %>%
  summarise(mean_ripple_rate = mean(ripple_rate, na.rm = TRUE), 
            AP_uncal = AP_uncal[1],
            hipp_or_amy = hipp_or_amy[1],
            AP = AP[1],
            .groups = "drop")

modulation_diff_df_elec <- modulation_df_elec %>%
  group_by(sub_num, electrode_num) %>%
  summarise(
    diff_mod = (mean(mean_ripple_rate[z_sextile %in% c(1, 2)]) -
                  mean(mean_ripple_rate[z_sextile %in% c(5, 6)])) / 
               (mean(mean_ripple_rate[z_sextile %in% c(1, 2)]) +
                mean(mean_ripple_rate[z_sextile %in% c(5, 6)])),
    
    diff_mod_less = (mean(mean_ripple_rate[z_sextile %in% c(1)]) -
                      mean(mean_ripple_rate[z_sextile %in% c(6)])) / 
                    (mean(mean_ripple_rate[z_sextile %in% c(1)]) +
                     mean(mean_ripple_rate[z_sextile %in% c(6)])),
    AP_uncal = AP_uncal[1],
    hipp_or_amy = hipp_or_amy[1],
    AP = AP[1],
    .groups = "drop"
  )

# 4.2 Pupil Sextile Modulation Analysis for Hippocampus
# ----

# Model AP effect on modulation
modulation_diff_df_elec_hipp <- modulation_diff_df_elec %>%
  filter(hipp_or_amy == 1)
modulation_diff_df_elec_hipp$sub_num[modulation_diff_df_elec_hipp$sub_num == 42] <- 49

model.pup1 <- lmer(diff_mod ~ AP + (1 | sub_num), data = modulation_diff_df_elec_hipp)
summary(model.pup1)
# 4.3 Data Preparation for Main Pupil Plots
# ----

# Define color schemes for plots
grey_shade <- "#bfbfc4"
pink_colors <- c(
  rgb(81, 12, 64, maxColorValue = 255),      # red
  rgb(120, 33, 100, maxColorValue = 255),      # green
  rgb(154, 59, 129, maxColorValue = 255),      # blue
  rgb(187, 89, 160, maxColorValue = 255),    # orange
  rgb(213, 123, 178, maxColorValue = 255),    # purple
  rgb(227, 161, 199, maxColorValue = 255)     # teal
)

orange_colors <- c(
  rgb(110, 14, 16, maxColorValue = 255),      # red
  rgb(148, 34, 30, maxColorValue = 255),      # green
  rgb(183, 62, 37, maxColorValue = 255),      # blue
  rgb(213, 92, 39, maxColorValue = 255),    # orange
  rgb(236, 128, 53, maxColorValue = 255),    # purple
  rgb(249, 163, 95, maxColorValue = 255)     # teal
)

blue_purple <- c(
  rgb(19, 27, 71, maxColorValue = 255),      # red
  rgb(34, 53, 106, maxColorValue = 255),      # green
  rgb(64, 81, 141, maxColorValue = 255),      # blue
  rgb(100, 113, 176, maxColorValue = 255),    # orange
  rgb(135, 146, 201, maxColorValue = 255),    # purple
  rgb(173, 179, 218, maxColorValue = 255)     # teal
)

purple_colors2 <- c(
  rgb(36, 16, 65, maxColorValue = 255),      # red
  rgb(53, 33, 84, maxColorValue = 255),      # green
  rgb(87, 63, 117, maxColorValue = 255),      # blue
  rgb(121, 95, 151, maxColorValue = 255),    # orange
  rgb(155, 128, 187, maxColorValue = 255),    # purple
  rgb(191, 163, 206, maxColorValue = 255)     # teal
)

# Calculate subject averages
subject_avg_data <- long_data_rate %>%
  group_by(sub_num, z_sextile, hipp_or_amy) %>%
  summarise(ripple_rate = mean(ripple_rate), .groups = 'drop')

# Calculate aggregate averages with SEM
agg_data_simple <- subject_avg_data %>%
  group_by(z_sextile, hipp_or_amy) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(color = ifelse(hipp_or_amy == 1, blue_purple[z_sextile], orange_colors[z_sextile]))

# Separate for hippocampus and amygdala
hipp_agg <- agg_data_simple %>% filter(hipp_or_amy == 1)
hipp_subject <- subject_avg_data %>% filter(hipp_or_amy == 1)

amy_agg <- agg_data_simple %>% filter(hipp_or_amy == 2)
amy_subject <- subject_avg_data %>% filter(hipp_or_amy == 2)

# Create A/P/Amy variable
long_data_rate_AP <- long_data_rate %>%
  mutate(A_or_P_or_amy = ifelse(hipp_or_amy == 2, 3, 
                                ifelse(AP_uncal == 0, 0, 1)))

long_data_rate_hipp <- long_data_rate_AP %>% filter(hipp_or_amy == 1)
long_data_rate_amy <- long_data_rate_AP %>% filter(hipp_or_amy == 2)

# Subject averages for anterior/posterior split
subject_avg_data_hipp <- long_data_rate_hipp %>%
  group_by(sub_num, z_sextile, AP_uncal) %>%
  summarise(ripple_rate = mean(ripple_rate), .groups = 'drop')

agg_data_hipp <- subject_avg_data_hipp %>%
  group_by(z_sextile, AP_uncal) %>%
  summarise(ripple_rate = mean(ripple_rate), .groups = 'drop')

agg_data_AP <- subject_avg_data_hipp %>% 
  group_by(z_sextile, AP_uncal) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

agg_data_hipp_A <- agg_data_AP %>% filter(AP_uncal == 0)
agg_data_hipp_P <- agg_data_AP %>% filter(AP_uncal == 1)

subject_avg_data_hipp_A <- subject_avg_data_hipp %>% filter(AP_uncal == 0)
subject_avg_data_hipp_P <- subject_avg_data_hipp %>% filter(AP_uncal == 1)

# 4.4 Plotting Function for Pupil Sextiles
# ----

simple_plot_final_nodots_indiv <- function(agg_data_use, y_limits, color_vector, 
                                           grey_shade, plot_title = "") {
  
  agg_data_use <- agg_data_use %>%
    mutate(color = ifelse(z_sextile2 == 7, grey_shade, color_vector[z_sextile2]))
  
  ggplot(agg_data_use) +
    geom_line(data = agg_data_use %>% filter(z_sextile2 != 7),
              aes(x = z_sextile2, y = ripple_rate_mean, group = 1, color = color),
              size = 1.5) +
    geom_errorbar(aes(x = z_sextile2, ymin = ripple_rate_mean - ripple_rate_sem, 
                      ymax = ripple_rate_mean + ripple_rate_sem),
                  width = 0.2, color = "black", size = 1) +
    geom_point(aes(x = z_sextile2, y = ripple_rate_mean, fill = color),
               color = "black", size = 5, shape = 21, stroke = 1.5) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = unique(agg_data_use$z_sextile2),
      labels = function(x) ifelse(x == 7, "blink", x)
    ) +
    labs(x = "Pupil z-score sextile", y = "Ripple Rate", title = plot_title) +
    theme_light() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white")
    ) +
    ylim(y_limits[1], y_limits[2])
}

# 4.5 Main Pupil Plots: Hippocampus and Amygdala
# ----

hipp_agg$z_sextile2 <- hipp_agg$z_sextile
hipp_final_simple <- simple_plot_final_nodots_indiv(
  agg_data_use = hipp_agg,
  y_limits = c(0.02, 0.06),
  color_vector = blue_purple,
  grey_shade = grey_shade,
  plot_title = "Hippocampus"
)
hipp_final_simple

amy_agg$z_sextile2 <- amy_agg$z_sextile
amy_final_simple <- simple_plot_final_nodots_indiv(
  agg_data_use = amy_agg,
  y_limits = c(0.02, 0.06),
  color_vector = orange_colors,
  grey_shade = grey_shade,
  plot_title = "Amygdala"
)
amy_final_simple


# 4.6 Anterior/Posterior Split Plot
# ----

agg_hipp_A <- agg_data_hipp_A %>%
  mutate(region = "hipp_A",
         line_color = ifelse(z_sextile == 0, grey_shade, purple_colors2[4]),
         dot_color = ifelse(z_sextile == 0, "#dec9e2", purple_colors2[z_sextile]))

agg_hipp_P <- agg_data_hipp_P %>%
  mutate(region = "hipp_P",
         line_color = ifelse(z_sextile == 0, grey_shade, pink_colors[4]),
         dot_color = ifelse(z_sextile == 0, '#efcee3', pink_colors[z_sextile]))

all_data <- bind_rows(agg_hipp_A, agg_hipp_P)

line_types <- c("hipp_A" = "solid", "hipp_P" = "dashed")

plot_all_sextiles_AP <- ggplot() +
  geom_line(data = all_data,
            aes(x = z_sextile, y = ripple_rate_mean, group = region,
                color = line_color, linetype = region),
            linewidth = 1.5, show.legend = FALSE) +
  geom_errorbar(data = all_data,
                aes(x = z_sextile, ymin = ripple_rate_mean - ripple_rate_sem,
                    ymax = ripple_rate_mean + ripple_rate_sem),
                width = 0.2, size = 1, show.legend = FALSE) +
  geom_point(data = all_data,
             aes(x = z_sextile, y = ripple_rate_mean, fill = dot_color),
             shape = 21, size = 5, stroke = 1.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linetype_manual(values = line_types) +
  labs(x = "Pupil z-score sextile", y = "Ripple Rate",
       title = "Ripple Rate by Region and Pupil Sextile") +
  theme_light() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  ) +
  ylim(0.02, 0.06)

plot_all_sextiles_AP

# 4.7 Main Pupil Models: Mixed Effects Analysis
# ----
# Model: Hippocampus only with anterior/posterior split
long_data_rate_hipp <- long_data_rate_AP %>% filter(hipp_or_amy == 1)
model.pup1 <- lmer(ripple_rate ~ z_sextile * AP_uncal + (1 | sub_num/electrode_num), 
                   data = long_data_rate_hipp)
summary(model.pup1)
test(emmeans(model.pup1, ~ z_sextile, at = list(z_sextile = c(1, 6))))
test(emtrends(model.pup1, ~ AP_uncal, var = "z_sextile")) #anterior and posterior

emtrends(model.pup1, pairwise ~ AP_uncal, var = "z_sextile")

test(emmeans(model.pup1, ~ z_sextile | AP_uncal, at = list(z_sextile = c(1, 6))))
#test(emmeans(model.pup1, ~ AP_uncal | z_sextile, at = list(z_sextile = c(1, 6))))
#emm <- emmeans(model.pup1, ~ z_sextile * AP_uncal, at = list(z_sextile = c(1, 6)))
#pairs(emm, by = "z_sextile")

pp <- plot(ggpredict(model.pup1, terms = c('z_sextile', 'AP_uncal'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp


# Model: Hippocampus and Amygdala interaction with pupil sextile
model.pup1 <- lmer(ripple_rate ~ z_sextile * hipp_or_amy + (1 | sub_num/electrode_num), 
                   data = long_data_rate_AP)
summary(model.pup1)
pp <- plot(ggpredict(model.pup1, terms = c('z_sextile', 'hipp_or_amy'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp

test(emtrends(model.pup1, ~ hipp_or_amy, var = "z_sextile"))
emtrends(model.pup1, pairwise ~ hipp_or_amy, var = "z_sextile")

# Model with calibrated AP values

#insert continuous AP then run model
long_data_rate_hipp <- long_data_rate_hipp %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, AP),
            by = c("sub_num", "electrode_num"))

model.pup1 <- lmer(ripple_rate ~ z_sextile * AP + (1 | sub_num/electrode_num), 
                   data = long_data_rate_hipp)
summary(model.pup1)

# ============================================================================
# 6. MOD SCORE ANALYSIS: Sleep Modulation Index by Pupil State
# ============================================================================

# 6.1 Modulation Index Calculation
# ----
long_data_rate_AP <- long_data_rate_AP %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, AP),
            by = c("sub_num", "electrode_num"))
modulation_df <- long_data_rate_AP %>%
  group_by(sub_num, z_sextile, A_or_P_or_amy) %>%
  summarise(mean_ripple_rate = mean(ripple_rate, na.rm = TRUE), .groups = "drop",
            AP = AP[1])

# Calculate symmetry index: (small_pupil_rate - large_pupil_rate) / (small_pupil_rate + large_pupil_rate)
modulation_diff_df <- modulation_df %>%
  group_by(sub_num, A_or_P_or_amy) %>%
  summarise(
    diff_mod = (mean(mean_ripple_rate[z_sextile %in% c(1, 2)]) -
                  mean(mean_ripple_rate[z_sextile %in% c(5, 6)])) / 
               (mean(mean_ripple_rate[z_sextile %in% c(1, 2)]) +
                mean(mean_ripple_rate[z_sextile %in% c(5, 6)])),
    
    diff_mod_simple = (mean(mean_ripple_rate[z_sextile %in% c(1)]) -
                        mean(mean_ripple_rate[z_sextile %in% c(6)])) / 
                      (mean(mean_ripple_rate[z_sextile %in% c(1)]) +
                       mean(mean_ripple_rate[z_sextile %in% c(6)])),
    AP = AP[1],
    .groups = "drop"
  )

modulation_diff_df_filtered <- modulation_diff_df %>%
  filter(A_or_P_or_amy %in% c(0, 1, 3)) %>%
  mutate(
    Region = factor(A_or_P_or_amy,
                    levels = c(0, 1, 3),
                    labels = c("Ant. Hipp", "Post. Hipp", "Amyg.")),
    line_color = case_when(
      A_or_P_or_amy == 0 ~ "#9a81bb",  # Ant Hipp
      A_or_P_or_amy == 1 ~ "#d57bb2",  # Post Hipp
      A_or_P_or_amy == 3 ~ "#ec8035"   # Amyg
    )
  )

# 6.2 Mod Score Summary Statistics
# ----

modulation_summary <- modulation_diff_df_filtered %>%
  group_by(Region, line_color) %>%
  summarise(
    mean_diff_mod = mean(diff_mod, na.rm = TRUE),
    sem_diff_mod = sd(diff_mod, na.rm = TRUE) / sqrt(n()),
    
    mean_diff_modS = mean(diff_mod_simple, na.rm = TRUE),
    sem_diff_modS = sd(diff_mod_simple, na.rm = TRUE) / sqrt(n()),
    AP = AP[1],
    .groups = "drop"
  )

# 6.3 Mod Score Plot
# ----

plot_diff_APamy_mod <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  
  # Lines connecting same subject across regions
  geom_line(
    data = modulation_diff_df_filtered,
    aes(x = Region, y = diff_mod, group = sub_num),
    color = "gray70",
    linewidth = 0.7,
    alpha = 0.3
  ) +
  
  # Subject points
  geom_jitter(
    data = modulation_diff_df_filtered,
    aes(x = Region, y = diff_mod, fill = line_color),
    shape = ifelse(modulation_diff_df_filtered$Region == "Ant. Hipp", 21, 
                   ifelse(modulation_diff_df_filtered$Region == "Post. Hipp", 23, 24)),
    size = 4,
    alpha = 0.3,
    width = 0,
    stroke = .2,
    show.legend = FALSE
  ) +
  
  # Mean ± SEM error bars
  geom_errorbar(
    data = modulation_summary,
    aes(x = Region,
        ymin = mean_diff_mod - sem_diff_mod,
        ymax = mean_diff_mod + sem_diff_mod),
    width = 0.15,
    color = "black",
    size = 1
  ) +
  
  # Mean points
  geom_point(
    data = modulation_summary,
    aes(x = Region, y = mean_diff_mod, fill = line_color),
    shape = ifelse(modulation_summary$Region == "Ant. Hipp", 21, 
                   ifelse(modulation_summary$Region == "Post. Hipp", 23, 24)),
    size = 4,
    stroke = 1.5,
    color = "black"
  ) +
  
  scale_fill_identity() +
  labs(
    x = NULL,
    y = "Ripple rate symmetry index",
    title = "Pupil modulation of ripple rate"
  ) +
  theme_light(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  ylim(-0.5, 1)

plot_diff_APamy_mod
