# ============================================================================
# Sleep ripple analysis - NREM modulation and anatomical variation
# Date last edited: 3/13/2025
# ============================================================================

# Load required packages
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

################################################################################
# SECTION 1: DATA LOADING AND PREPROCESSING
################################################################################

# Load raw ripple data across all channels and nights
ripple_array_sleep <- read_csv("data/all_chan_ripple_final_9242025.csv", col_names = FALSE)
colnames(ripple_array_sleep) <- c(
  "sub_num", "electrode_num", "sleep_night", "block_num", 
  "ripple_index", "ripple_max_time", "keep_index", "ripple_time1", "ripple_time2", 
  "ripple_time3", "duration", "amplitude", "freq", "sleep_stage", "N1_time", "N2_time", "N3_time","R_time","W_time"
)

# Load electrode quality scores and anatomical classifications
sleep_elecs <- read_csv("data/sleep_electrode_array_updated_final_9242025.csv", col_names = FALSE)
colnames(sleep_elecs) <- c(
 "sub_num", "electrode_num", "keep_or_reject", "hipp_or_amy", "sleep_mod_score", "sleep_amp_score", "sleep_dur_score", "sleep_freq_score"
)

# Merge electrode quality information with ripple data
ripple_array_sleep <- ripple_array_sleep %>%
 left_join(sleep_elecs %>%
             select(sub_num, electrode_num, sleep_mod_score, keep_or_reject),
           by = c("sub_num", "electrode_num"))

# Load electrode anatomical positions (MNI coordinates, anterior-posterior location, region)
sub_elec_positions <- read_csv("data/sub_elec_positions_allinfo_2025.csv", col_names = FALSE)
colnames(sub_elec_positions) <- c(
  "sub_num", "electrode_num", "left_or_right", "hipp_or_amy", "AP", "RL", "SI","AP_uncal","apex_coord","elec_coord","coord_diff","close_flag","MNI_coord","post_50","post_35","post_MNI","overlap_50","overlap_MNI","overlap_35"
)

# Merge electrode position data with ripple dataset
ripple_array_sleep <- ripple_array_sleep %>%
  left_join(sub_elec_positions %>%
              select(sub_num, electrode_num, hipp_or_amy, left_or_right, AP, AP_uncal), 
            by = c("sub_num", "electrode_num"))

# Apply subject and session-level exclusions (based on data quality review)
ripple_array_sleep <- ripple_array_sleep %>%
  filter(
    sub_num != 56,
    !(sub_num == 44 & sleep_night == 1),
    !(sub_num == 49 & sleep_night == 1),
    !(sub_num == 54 & sleep_night == 3),
    !(sub_num == 55 & sleep_night == 3),
    !(sub_num == 57 & (sleep_night == 3 | sleep_night == 4)),
    !(sub_num == 64 & sleep_night == 1),
    !(sub_num == 69 & sleep_night == 2),
    !(sub_num == 62 & electrode_num == 134),
    !(sub_num == 49 & electrode_num == 38)
  )

################################################################################
# SECTION 2: DATA FILTERING AND QUALITY CONTROL
################################################################################

# Keep only detected ripple events that passed manual validation
ripple_array_final_sleep <- filter(ripple_array_sleep, keep_index == 1)
# Keep only ripples from high-quality electrodes
ripple_array_final_sleep <- filter(ripple_array_final_sleep, keep_or_reject == 0)

# Remove artifactual ripples with excessively high amplitudes (> 300 µV)
ripple_array_final_sleep <- ripple_array_final_sleep %>%
  filter(amplitude < 300)

# Consolidate sleep stages: combine N1, N2, N3 into single "NREM" category for analysis
# (N1 is brief and unstable; combining with deeper NREM stages for more power)
ripple_array_final_sleep <- ripple_array_final_sleep %>%
  mutate(sleep_stage2 = case_when(
    sleep_stage == -1 ~ NA_real_,     # Convert -1 to NA
    sleep_stage %in% c(-3, -2) ~ -1,  # Convert -3 and -2 to -1 (consolidated NREM)
    TRUE ~ sleep_stage                # Keep REM (0) and Wake (1) unchanged
  ))

################################################################################
# SECTION 3: RIPPLE RATE CALCULATION BY STAGE AND ELECTRODE
################################################################################

# Calculate stage-specific time spent in each sleep stage per session
ripple_array_sleep_rate <- ripple_array_final_sleep %>%
  mutate(stage_time = case_when(
    sleep_stage == -3 ~ N3_time,
    sleep_stage == -2 ~ N2_time,
    sleep_stage == -1 ~ N1_time,
    sleep_stage == 0  ~ R_time,
    sleep_stage == 1  ~ W_time
  ), stage_time2 = case_when(
    sleep_stage2 == -1 ~ N3_time+N2_time,
    sleep_stage2 == 0 ~ R_time,
    sleep_stage2 == 1 ~ W_time))

# Compute ripple rate (events per minute) for each electrode × stage × session
ripple_rate_sleep <- ripple_array_sleep_rate%>%
  group_by(sub_num, electrode_num, sleep_night, sleep_stage) %>%
  summarize(
    hipp_or_amy = first(hipp_or_amy),
    event_count = n(),
    stage_time = stage_time[1],
    mod_score = sleep_mod_score[1],
    AP = AP[1],
    AP_uncal = AP_uncal[1]) %>%
  mutate(ripple_rate = event_count/stage_time,
         A_or_P = ifelse(AP>.5,2,1),
         A_or_P_uncal = ifelse(AP_uncal==1,2,ifelse(AP_uncal==0,1)),
         A_or_P_or_amy = ifelse(is.nan(AP_uncal),2,AP_uncal))

# Compute ripple rate using consolidated NREM stage grouping
ripple_rate_sleep2 <- ripple_array_sleep_rate%>%
  group_by(sub_num, electrode_num, sleep_night, sleep_stage2) %>%
  summarize(
    hipp_or_amy = first(hipp_or_amy),
    event_count2 = n(),
    stage_time2 = stage_time2[1],
    AP = AP[1],
    AP_uncal = AP_uncal[1]) %>%
  mutate(ripple_rate2 = event_count2/stage_time2,
         A_or_P = ifelse(AP>.5,2,1),
         A_or_P_uncal = ifelse(AP_uncal==1,2,ifelse(AP_uncal==0,1)),
         A_or_P_or_amy = ifelse(is.nan(AP_uncal),2,AP_uncal))

# Remove N1-only observations (consolidated into NREM category)
ripple_rate_sleep2 <- filter(ripple_rate_sleep2, !is.na(sleep_stage2))

################################################################################
# SECTION 4: STATISTICAL MODELS - RIPPLE RATE ACROSS SLEEP STAGES
################################################################################

# Account for subject overlap: subjects 42 and 49 have non-overlapping electrodes
# (consolidate as same subject for random effects structure)
ripple_rate_sleep$sub_num[ripple_rate_sleep$sub_num == 42] <- 49
ripple_rate_sleep2$sub_num[ripple_rate_sleep2$sub_num == 42] <- 49

# Extract hippocampal ripples only for anatomical subdivision analyses
ripple_rate_sleep_hipp <- ripple_rate_sleep %>%
  filter(hipp_or_amy == 1)
ripple_rate_sleep2_hipp <- ripple_rate_sleep2 %>%
  filter(hipp_or_amy == 1)

# MODEL 1: Anterior-Posterior dissociation within hippocampus
# Test for differential modulation by sleep stage across anterior and posterior hippocampus
ripple_rate_sleep2_hipp$A_or_P_uncal <- as.factor(ripple_rate_sleep2_hipp$A_or_P_uncal)
ripple_rate_sleep2_hipp$sleep_stage2 <- as.factor(ripple_rate_sleep2_hipp$sleep_stage2)
model.sleep1AP <- lmer(ripple_rate2 ~ sleep_stage2*A_or_P_uncal + (1 | sub_num/electrode_num), REML=F, data = ripple_rate_sleep2_hipp)

# Print model summary and main effects
summary(model.sleep1AP) 
car::Anova(model.sleep1AP, type='III')

# Post-hoc comparisons for overall hippocampus effect
emmeans(model.sleep1AP,  ~ sleep_stage2)
pairs(emmeans(model.sleep1AP,  ~ sleep_stage2))

# Evaluate anterior vs posterior differences
joint_tests(model.sleep1AP, by = "A_or_P_uncal")
emmeans(model.sleep1AP, pairwise ~ sleep_stage2 | A_or_P_uncal, lmer.df='satterthwaite')
emmeans(model.sleep1AP, pairwise ~ A_or_P_uncal | sleep_stage2, lmer.df='satterthwaite')

# MODEL 2: Hippocampus vs Amygdala comparison
# Test for regional differences in ripple rate modulation across sleep stages
ripple_rate_sleep2$hipp_or_amy <- as.factor(ripple_rate_sleep2$hipp_or_amy)
ripple_rate_sleep2$sleep_stage2 <- as.factor(ripple_rate_sleep2$sleep_stage2)
model.sleep1AP <- lmer(ripple_rate2 ~ sleep_stage2*hipp_or_amy + (1 | sub_num/electrode_num), REML=F, data = ripple_rate_sleep2)

summary(model.sleep1AP)
car::Anova(model.sleep1AP, type='III')
emmeans(model.sleep1AP,  ~ sleep_stage2 | hipp_or_amy)
joint_tests(model.sleep1AP, by = "hipp_or_amy")
emmeans(model.sleep1AP, pairwise ~ sleep_stage2 | hipp_or_amy, lmer.df='satterthwaite')
emmeans(model.sleep1AP, pairwise ~ hipp_or_amy | sleep_stage2, lmer.df='satterthwaite')

# MODEL 3: Continuous anterior-posterior position within hippocampus
# Test for graded effect of A-P position (not categorical, but continuous)
model.sleep1 <- lmer(ripple_rate2 ~ sleep_stage2*AP + (1 | sub_num/electrode_num), REML=F, data = ripple_rate_sleep2_hipp)

summary(model.sleep1)
car::Anova(model.sleep1, type='III')
test(emtrends(model.sleep1, ~ sleep_stage2, var = "AP", lmer.df = "satterthwaite"))

################################################################################
# SECTION 5: MODULATION INDEX CALCULATION AND STATS TEST
################################################################################

# Calculate modulation score: compare NREM vs REM ripple rates for each electrode
# Positive values = higher in NREM; negative = higher in REM

mod_scores <- ripple_rate_sleep2_hipp %>%
  filter(sleep_stage2 %in% c(-1, 1)) %>%
  group_by(sub_num, electrode_num, AP, sleep_stage2) %>%
  summarize(
    mean_rate = mean(ripple_rate2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sleep_stage2,
    values_from = mean_rate,
    names_prefix = "rate_"
  ) %>%
  rename(nrem_rate = `rate_-1`, rem_rate = rate_1) %>%
  mutate(
    mod_index = (nrem_rate - rem_rate) / (nrem_rate + rem_rate)
  )

# Fit mixed effects model: does modulation vary by anterior-posterior position?
# Random intercept per subject and electrode (nested structure)
mod_lmer <- lmer(mod_index ~ AP + (1 | sub_num), 
                 data = mod_scores, REML = FALSE)
summary(mod_lmer)

################################################################################
# SECTION 5: RIPPLE RATE DATA PREPARATION (FOR VISUALIZATION PLOTS)
################################################################################

# Prepare NREM-consolidated rate data for hippocampus (separate anterior/posterior)
df_summary_rate_hipp2 <- ripple_rate_sleep2 %>%
  filter(!is.na(sleep_stage2), hipp_or_amy == 1) %>%
  group_by(sub_num, electrode_num, sleep_stage2) %>%
  summarise(avg_rate = mean(ripple_rate2, na.rm = TRUE),
            A_or_P_uncal = A_or_P_uncal[1])

df_subject_rate_hipp2 <- df_summary_rate_hipp2 %>%
  group_by(sub_num, sleep_stage2) %>%
  summarise(avg_rate = mean(avg_rate, na.rm = TRUE), .groups = 'drop')

df_subject_rate_hipp2_AP <- df_summary_rate_hipp2 %>%
  group_by(sub_num, A_or_P_uncal, sleep_stage2) %>%
  summarise(avg_rate = mean(avg_rate, na.rm = TRUE), .groups = 'drop')

df_mean_rate_hipp2 <- df_subject_rate_hipp2 %>%
  group_by(sleep_stage2) %>%
  summarise(
    avg_rate2 = mean(avg_rate, na.rm = TRUE),
    sem = sd(avg_rate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

df_mean_rate_hipp2_AP <- df_subject_rate_hipp2_AP %>%
  group_by(A_or_P_uncal, sleep_stage2) %>%
  summarise(
    avg_rate2 = mean(avg_rate, na.rm = TRUE),
    sem = sd(avg_rate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Prepare full sleep stage rate data (N3, N2, N1, REM, Wake separate)
df_summary_rate_hipp <- ripple_rate_sleep %>%
  filter(!is.nan(sleep_stage), hipp_or_amy == 1) %>%
  group_by(sub_num, electrode_num, A_or_P_uncal, sleep_stage) %>%
  summarise(avg_rate = mean(ripple_rate, na.rm = TRUE), .groups = 'drop')

df_subject_rate_hipp_full <- df_summary_rate_hipp %>%
  group_by(sub_num, sleep_stage, A_or_P_uncal) %>%
  summarise(avg_rate = mean(avg_rate, na.rm = TRUE), .groups = 'drop')

df_mean_rate_hipp_full <- df_subject_rate_hipp_full %>%
  group_by(sleep_stage, A_or_P_uncal) %>%
  summarise(avg_rate2 = mean(avg_rate, na.rm = TRUE),
            sem = sd(avg_rate, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop')

# Prepare amygdala data for full sleep stages
df_summary_rate_amy <- ripple_rate_sleep %>%
  filter(!is.nan(sleep_stage), hipp_or_amy == 2) %>%
  group_by(sub_num, electrode_num, sleep_stage) %>%
  summarise(avg_rate = mean(ripple_rate, na.rm = TRUE), .groups = 'drop')

df_subject_rate_amy <- df_summary_rate_amy %>%
  group_by(sub_num, sleep_stage) %>%
  summarise(avg_rate = mean(avg_rate, na.rm = TRUE), .groups = 'drop')

df_mean_rate_amy <- df_subject_rate_amy %>%
  group_by(sleep_stage) %>%
  summarise(avg_rate2 = mean(avg_rate, na.rm = TRUE),
            sem = sd(avg_rate, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop')

################################################################################
# SECTION 6: PLOTTING FUNCTIONS
################################################################################
plot_sleep_data_errorbar_simple <- function(df_subject, df_mean, 
                                            sleep_stage_info2, sleep_stage_info3,
                                            sleep_stage_levels, sleep_stage_labels, 
                                            y_label2, y_label3,
                                            y_info, sleep_stage_colors, y_text, title_text) {
  ggplot() +
    # Mean line
    geom_line(data = df_mean,
              aes(x = factor({{ sleep_stage_info3 }}, levels = sleep_stage_levels, labels = sleep_stage_labels),
                  y = {{ y_label3 }},
                  group = 1),
              color = "black", linewidth = 1) +
    # SEM error bars
    geom_errorbar(data = df_mean,
                  aes(x = factor({{ sleep_stage_info3 }}, levels = sleep_stage_levels, labels = sleep_stage_labels),
                      ymin = {{ y_label3 }} - sem,
                      ymax = {{ y_label3 }} + sem),
                  #color = factor({{ sleep_stage_info3 }}, levels = sleep_stage_levels)),
                  width = 0.2, size = 1) + # color = "black", size = 1)+
    # Mean points
    geom_point(data = df_mean,
               aes(x = factor({{ sleep_stage_info3 }}, levels = sleep_stage_levels, labels = sleep_stage_labels),
                   y = {{ y_label3 }},
                   #color = factor({{ sleep_stage_info3 }}, levels = sleep_stage_levels),
                   fill = factor({{ sleep_stage_info3 }}, levels = sleep_stage_levels)),size = 4, shape = 21, stroke = 1.5)+
    
    scale_color_manual(values = sleep_stage_colors) +
    scale_fill_manual(values = sleep_stage_colors) +
    labs(x = 'Sleep Stage', y = y_text, title = title_text) +
    theme_minimal() +
    ylim(y_info[1], y_info[2]) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"))
}

plot_sleep_data_combined <- function(df_mean, sleep_stage_info_col,
                                     sleep_stage_levels, sleep_stage_labels, 
                                     y_col, y_limits, fill_colors, y_text, title_text) {
  ggplot(df_mean,
         aes(x = factor({{ sleep_stage_info_col }}, levels = sleep_stage_levels, labels = sleep_stage_labels),
             y = {{ y_col }},
             group = region)) +
    
    # Line with linetype and color
    geom_line(aes(linetype = region,
                  color = ifelse(region == "P", "black", "black")),
              linewidth = 1) +
    
    # Error bars
    geom_errorbar(aes(ymin = {{ y_col }} - sem,
                      ymax = {{ y_col }} + sem,
                      color = ifelse(region == "P", "black", "black")),
                  width = 0.2, size = 1) +
    
    # Points with custom fill
    geom_point(aes(fill = fill_group),
               size = 5, shape = 21, stroke = 1.5) +
    
    scale_fill_manual(values = fill_colors) +
    scale_linetype_manual(values = c("A" = "solid", "P" = "dashed", "amy" = "dotted")) +
    scale_color_identity() +
    
    labs(x = 'Sleep Stage', y = y_text, title = title_text) +
    theme_minimal() +
    ylim(y_limits[1], y_limits[2]) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"))
}


################################################################################
# SECTION 7: VISUALIZATION - RIPPLE RATE ACROSS SLEEP STAGES
################################################################################

# Define sleep stage labels and colors for NREM-consolidated plots
sleep_stage_labels2 <- c('NREM', 'REM', 'W')
sleep_stage_colors2 <- c(
  rgb(32, 120, 180, maxColorValue = 255),     # NREM (blue)
  rgb(97, 173, 205, maxColorValue = 255),     # REM (light blue)
  rgb(191, 191, 196, maxColorValue = 255)     # Wake (grey)
)

# Define sleep stage labels and colors for all individual sleep stages
sleep_stage_labels1 <- c('N3', 'N2', 'N1', 'REM', 'W')
sleep_stage_colors_final <- c(
  rgb(8, 29, 88, maxColorValue = 255),     # N3
  rgb(36, 68, 142, maxColorValue = 255),   # N2
  rgb(70, 118, 224, maxColorValue = 255),  # N1
  
  rgb(104, 193, 176, maxColorValue = 255),   # REM
  rgb(191, 191, 196, maxColorValue = 255)  # W
)

## plot hippocampus simple
hipp_sum_rate_simple <- plot_sleep_data_errorbar_simple(df_subject_rate_hipp2, df_mean_rate_hipp2, 
                                                        sleep_stage2, sleep_stage2,
                                                        c(-1:1), sleep_stage_labels2,
                                                        avg_rate, avg_rate2,
                                                        c(0.02,.07), sleep_stage_colors2, "Rate (hz)", "Hippocampus")
hipp_sum_rate_simple

#### now plot anterior and posterior hippocampus overlaid
df_mean_rate_hipp2_A <- df_mean_rate_hipp2_AP %>% filter(A_or_P_uncal == 1)
df_mean_rate_hipp2_P <- df_mean_rate_hipp2_AP %>% filter(A_or_P_uncal == 2)

df_mean_rate_hipp2_A$region <- "A"
df_mean_rate_hipp2_P$region <- "P"
df_mean_combined <- rbind(df_mean_rate_hipp2_A, df_mean_rate_hipp2_P)
df_mean_combined$fill_group <- interaction(df_mean_combined$region, df_mean_combined$sleep_stage2)
fill_colors <- c(
  setNames(sleep_stage_colors2, paste0("A.", -1:1)),
  setNames(sleep_stage_colors2, paste0("P.", -1:1))
)


hipp_combined_plot <- plot_sleep_data_combined(df_mean_combined, 
                                               sleep_stage_info_col = sleep_stage2,
                                               sleep_stage_levels = c(-1:1),
                                               sleep_stage_labels = sleep_stage_labels2,
                                               y_col = avg_rate2,
                                               y_limits = c(0.02, 0.07),
                                               fill_colors = fill_colors,
                                               y_text = "Rate (hz)",
                                               title_text = "Anterior vs Posterior Hippocampus")

hipp_combined_plot

## ===== PLOT 2: All Three Regions with Full Sleep Stages (N3, N2, N1, REM, Wake) =====
# Prepare combined data for anterior hipp, posterior hipp, and amygdala (all sleep stages)
df_mean_rate_hipp_full_A <- df_mean_rate_hipp_full %>% filter(A_or_P_uncal == 1)
df_mean_rate_hipp_full_P <- df_mean_rate_hipp_full %>% filter(A_or_P_uncal == 2)

df_mean_rate_hipp_full_A$region <- "A"
df_mean_rate_hipp_full_P$region <- "P"
df_mean_rate_amy$region <- "amy"

df_mean_combined_full <- bind_rows(
  df_mean_rate_hipp_full_A,
  df_mean_rate_hipp_full_P,
  df_mean_rate_amy
)

df_mean_combined_full$fill_group <- interaction(df_mean_combined_full$region, df_mean_combined_full$sleep_stage)

fill_colors_full <- c(
  setNames(sleep_stage_colors_final, paste0("A.", -3:1)),
  setNames(sleep_stage_colors_final, paste0("P.", -3:1)),
  setNames(sleep_stage_colors_final, paste0("amy.", -3:1))
)

# Create plot: All three regions across all sleep stages
plot_all_regions_full <- plot_sleep_data_combined(df_mean_combined_full, 
                                                  sleep_stage_info_col = sleep_stage,
                                                  sleep_stage_levels = c(-3, -2, -1, 0, 1),
                                                  sleep_stage_labels = sleep_stage_labels1,
                                                  y_col = avg_rate2,
                                                  y_limits = c(0.02, 0.075),
                                                  fill_colors = fill_colors_full,
                                                  y_text = "Rate (hz)",
                                                  title_text = "Ripple Rate Across Sleep Stages")
plot_all_regions_full

################################################################################
# SECTION 8: VISUALIZATION - NREM MODULATION BY REGION (Subject-level)
################################################################################
# Calculate NREM vs REM modulation directly by subject and region
modulation_diff_df <- ripple_rate_sleep2 %>%
  group_by(sub_num, sleep_stage2, A_or_P_or_amy) %>%
  summarise(mean_ripple_rate = mean(ripple_rate2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sleep_stage2,
    values_from = mean_ripple_rate,
    names_prefix = "rate_"
  ) %>%
  rename(nrem_rate = `rate_-1`, rem_rate = rate_1) %>%
  mutate(
    diff_mod = (nrem_rate - rem_rate) / (nrem_rate + rem_rate)
  ) %>%
  filter(A_or_P_or_amy %in% c(0, 1, 2)) %>%
  mutate(
    Region = factor(A_or_P_or_amy,
                    levels = c(0, 1, 2),
                    labels = c("Ant. Hipp", "Post. Hipp", "Amyg.")),
    line_color = case_when(
      A_or_P_or_amy == 0 ~ "#9a81bb",  # Anterior hippocampus
      A_or_P_or_amy == 1 ~ "#d57bb2",  # Posterior hippocampus
      A_or_P_or_amy == 2 ~ "#ec8035"   # Amygdala – orange
    )
  )

# Calculate summary statistics for means and error bars
modulation_summary <- modulation_diff_df %>%
  group_by(Region, line_color) %>%
  summarise(
    mean_diff_mod = mean(diff_mod, na.rm = TRUE),
    sem_diff_mod = sd(diff_mod, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

plot_diff_APamy_mod <- ggplot() +
  # Individual subject points
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  # Lines connecting the same subject across regions
  geom_line(
    data = modulation_diff_df,
    aes(x = Region, y = diff_mod, group = sub_num),
    color = "gray70",
    linewidth = 0.7,
    alpha = 0.3
  ) +
  geom_jitter(
    data = modulation_diff_df,
    aes(x = Region, y = diff_mod, fill = line_color),
    shape = ifelse(modulation_diff_df$Region == "Ant. Hipp", 21, ifelse(modulation_diff_df$Region == "Post. Hipp", 23, 24)),
    size = 4,
    alpha = 0.3,
    width = 0,
    stroke = .2,
    show.legend = FALSE
  ) +
  # Mean ± SEM
  geom_errorbar(
    data = modulation_summary,
    aes(x = Region,
        ymin = mean_diff_mod - sem_diff_mod,
        ymax = mean_diff_mod + sem_diff_mod),
    width = 0.15,
    color = "black",
    size = 1
  ) +
  geom_point(
    data = modulation_summary,
    aes(x = Region, y = mean_diff_mod, fill = line_color),
    shape = ifelse(modulation_summary$Region == "Ant. Hipp", 21, ifelse(modulation_summary$Region == "Post. Hipp", 23, 24)),
    size = 4,
    stroke = 1.5,
    color = "black"
  ) +
  scale_fill_identity() +
  labs(
    x = NULL,
    y = "Ripple rate symmetry index",
    title = "NREM modulation of ripple rate"
  ) +
  theme_light(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  ylim(-.75,.75)
plot_diff_APamy_mod
