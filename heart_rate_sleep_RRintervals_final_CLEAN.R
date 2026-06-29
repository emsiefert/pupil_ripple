# ============================================================================
# Heart rate analyses for RR intervals during sleep recordings
# mixed effects models (main manuscript), plots (main manuscript), modulation scores (in supplement).
# everything has all sleep (analysis on full recording, across all stages), NREM only, and wake only.
# last updated April/June 2026, ES
# ============================================================================

# LIBRARIES
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
# LOAD DATA
# ============================================================================

heart_rate_sleep_50 <- read.csv("array_for_R_RR_50hz_sleep.csv", header = FALSE)

colnames(heart_rate_sleep_50) <- c(
  "sub_num", "electrode_num", "hipp_or_amy",
  "A_or_P", "sextile", "row_count", "ripple_count", "average_RR",
  "row_countNREM", "ripple_countNREM", "average_RRNREM",
  "row_countW", "ripple_countW", "averageRRW",
  "row_countREM", "ripple_countREM", "averageRRREM"
)

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Set up ripple rate arrays
heart_rate_sleep_50 <- heart_rate_sleep_50 %>% #make values numeric
  mutate(
    RR_sextile = as.numeric(sextile),
    ripple_count = as.numeric(ripple_count),
    row_count = as.numeric(row_count),
    row_countNREM = as.numeric(row_countNREM),
    row_countW = as.numeric(row_countW),
    row_countREM = as.numeric(row_countREM),
    ripple_countNREM = as.numeric(ripple_countNREM),
    ripple_countW = as.numeric(ripple_countW),
    sub_num = as.factor(sub_num)
  )

# Calculate ripple rate, from heart rate that was interpolated to 50 Hz
heart_rate_sleep_50 <- heart_rate_sleep_50 %>%
  mutate(
    ripple_rate = ripple_count / (row_count / 50),
    ripple_rate_NREM = ripple_countNREM / (row_countNREM / 50),
    ripple_rate_W = ripple_countW / (row_countW / 50)
  )

heart_rate_sleep_50 <- heart_rate_sleep_50 %>%
  mutate(A_or_P_or_amy = ifelse(hipp_or_amy == 2, 2, A_or_P))

# Apply exclusions - from sleep data. Removing these electrodes given large artifacts observed during sleep.
heart_rate_sleep_50 <- heart_rate_sleep_50 %>%
  filter(
    !(sub_num == 62 & electrode_num == 134),
    !(sub_num == 49 & electrode_num == 38)
  )

# Load and merge electrode position data
sub_elec_positions <- read_csv("sub_elec_positions_allinfo_2025.csv", col_names = FALSE)
colnames(sub_elec_positions) <- c(
  "sub_num", "electrode_num", "left_or_right", "hipp_or_amy", "AP", "RL", "SI",
  "AP_uncal", "apex_coord", "elec_coord", "coord_diff", "close_flag",
  "MNI_coord", "post_50", "post_35", "post_MNI", "overlap_50", "overlap_MNI", "overlap_35"
)

sub_elec_positions$sub_num <- as.factor(sub_elec_positions$sub_num)

heart_rate_sleep_50 <- heart_rate_sleep_50 %>%
  left_join(sub_elec_positions %>% select(sub_num, electrode_num, AP),
    by = c("sub_num", "electrode_num"))

# ============================================================================
# MIXED EFFECTS MODELS
# ============================================================================

# Fix subject number, two implants but same participant. Non overlapping electrodes, but same participant for nested models
heart_rate_sleep_50$sub_num[heart_rate_sleep_50$sub_num == 49] <- 42

# Filter for hippocampus
heart_rate_sleep_50_hipp <- heart_rate_sleep_50 %>%
  filter(hipp_or_amy == 1)

# Reverse sextiles for model interpretability!!!! This reversal is important for interpreting the reporting in the paper and aligning it with the results printed here.
heart_rate_sleep_50_hipp$RR_sextile_rev <- 7 - heart_rate_sleep_50_hipp$RR_sextile
heart_rate_sleep_50$RR_sextile_rev <- 7 - heart_rate_sleep_50$RR_sextile

# ============================================================================
# ALL SLEEP MODELS - HIPPOCAMPUS
# ============================================================================

# Model 1: All sleep, Hippocampus, A/P as categorical
heart_rate_sleep_50_hipp$A_or_P <- as.factor(heart_rate_sleep_50_hipp$A_or_P)
model.rate.hippAP <- lmer(ripple_rate ~ RR_sextile_rev * A_or_P + (1 | sub_num/electrode_num),
                               REML = F, data = heart_rate_sleep_50_hipp)
summary(model.rate.hippAP)
car::Anova(model.rate.hippAP, type = 'III')
test(emmeans(model.rate.hippAP, ~ RR_sextile_rev, at = list(RR_sextile_rev = c(1, 6)))) #full hipp
emtrends(model.rate.hippAP, ~ A_or_P, var = "RR_sextile_rev") #anterior and posterior
test(emtrends(model.rate.hippAP, ~ A_or_P, var = "RR_sextile_rev"))
emtrends(model.rate.hippAP, pairwise ~ A_or_P, var = "RR_sextile_rev")
test(emmeans(model.rate.hippAP, pairwise ~ A_or_P | RR_sextile_rev,
        at = list(RR_sextile_rev = c(1, 6))))

pp <- plot(ggpredict(model.rate.hippAP, terms = c('RR_sextile_rev', 'A_or_P'))) + #plot model results
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp

# Model 2: All sleep, Hippocampus, AP as continuous
model.EOEC.rate.hippAP <- lmer(ripple_rate ~ RR_sextile_rev * AP + (1 | sub_num/electrode_num),
                               REML = F, data = heart_rate_sleep_50_hipp)
summary(model.EOEC.rate.hippAP)

pp <- plot(ggpredict(model.EOEC.rate.hippAP, terms = c('RR_sextile_rev', 'AP'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp

# Model 3: All sleep, Hipp vs Amy

model.EOEC.rate.hippamy <- lmer(ripple_rate ~ RR_sextile_rev * hipp_or_amy + (1 | sub_num/electrode_num),
                                REML = F, data = heart_rate_sleep_50)
summary(model.EOEC.rate.hippamy)
car::Anova(model.EOEC.rate.hippamy, type = 'III')
emtrends(model.EOEC.rate.hippamy, ~ hipp_or_amy, var = "RR_sextile_rev")
test(emtrends(model.EOEC.rate.hippamy, ~ hipp_or_amy, var = "RR_sextile_rev"))
emtrends(model.EOEC.rate.hippamy, pairwise ~ hipp_or_amy, var = "RR_sextile_rev")
test(emmeans(model.EOEC.rate.hippamy, pairwise ~ hipp_or_amy | RR_sextile_rev,
             at = list(RR_sextile_rev = c(1, 6))))

# ============================================================================
# NREM ONLY MODELS
# ============================================================================

# Model 1: NREM, Hipp vs Amy
model.NREM2 <- lmer(ripple_rate_NREM ~ RR_sextile_rev * hipp_or_amy + (1 | sub_num/electrode_num),
                    REML = F, data = heart_rate_sleep_50)
summary(model.NREM2)
car::Anova(model.NREM2, type = 'III')
pp <- plot(ggpredict(model.NREM2, terms = c('RR_sextile_rev', 'hipp_or_amy'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp
emtrends(model.NREM2, ~ hipp_or_amy, var = "RR_sextile_rev")
test(emtrends(model.NREM2, ~ hipp_or_amy, var = "RR_sextile_rev"))
emtrends(model.NREM2, pairwise ~ hipp_or_amy, var = "RR_sextile_rev")

# Model 2: NREM, Hippocampus, A/P as categorical
model.NREM_hipp <- lmer(ripple_rate_NREM ~ RR_sextile_rev * A_or_P + (1 | sub_num/electrode_num),
                        REML = F, data = heart_rate_sleep_50_hipp)
summary(model.NREM_hipp)
pp <- plot(ggpredict(model.NREM_hipp, terms = c('RR_sextile_rev', 'A_or_P'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp
emtrends(model.NREM_hipp, ~ A_or_P, var = "RR_sextile_rev")
test(emtrends(model.NREM_hipp, ~ A_or_P, var = "RR_sextile_rev"))
emtrends(model.NREM_hipp, pairwise ~ A_or_P, var = "RR_sextile_rev")
test(emmeans(model.NREM_hipp, pairwise ~ A_or_P | RR_sextile_rev,
             at = list(RR_sextile_rev = c(1, 6))))


# ============================================================================
# WAKE ONLY MODELS
# ============================================================================

# Model 1: Wake, Hipp vs Amy
model.wake2 <- lmer(ripple_rate_W ~ RR_sextile_rev * hipp_or_amy + (1 | sub_num/electrode_num),
                    REML = F, data = heart_rate_sleep_50)
summary(model.wake2)
car::Anova(model.wake2, type = 'III')
pp <- plot(ggpredict(model.wake2, terms = c('RR_sextile_rev', 'hipp_or_amy'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp
emtrends(model.wake2, ~ hipp_or_amy, var = "RR_sextile_rev")
test(emtrends(model.wake2, ~ hipp_or_amy, var = "RR_sextile_rev"))
emtrends(model.wake2, pairwise ~ hipp_or_amy, var = "RR_sextile_rev")

# Model 2: Wake, Hippocampus, A/P as categorical
model.wakeHipp <- lmer(ripple_rate_W ~ RR_sextile_rev * A_or_P + (1 | sub_num/electrode_num),
                       REML = F, data = heart_rate_sleep_50_hipp)
summary(model.wakeHipp)
car::Anova(model.wakeHipp, type = 'III')
pp <- plot(ggpredict(model.wakeHipp, terms = c('RR_sextile_rev', 'A_or_P'))) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
pp
emtrends(model.wakeHipp, ~ A_or_P, var = "RR_sextile_rev")
test(emtrends(model.wakeHipp, ~ A_or_P, var = "RR_sextile_rev"))
emtrends(model.wakeHipp, pairwise ~ A_or_P, var = "RR_sextile_rev")
test(emmeans(model.wakeHipp, pairwise ~ A_or_P | RR_sextile_rev,
             at = list(RR_sextile_rev = c(1, 6))))


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

simple_plot_final_nodots_indiv <- function(agg_data_use, y_limits,
                                           color_vector, grey_shade, plot_title = "") {
  agg_data_use <- agg_data_use %>%
    mutate(color = ifelse(sextile == 7, grey_shade, color_vector[sextile]))

  ggplot() +
    geom_line(data = agg_data_use %>% filter(sextile != 7),
              aes(x = sextile, y = ripple_rate_mean, group = 1, color = color),
              size = 1.5) +
    geom_errorbar(data = agg_data_use,
                  aes(x = sextile, ymin = ripple_rate_mean - ripple_rate_sem,
                      ymax = ripple_rate_mean + ripple_rate_sem),
                  width = 0.2, color = "black", size = 1) +
    geom_point(data = agg_data_use,
               aes(x = sextile, y = ripple_rate_mean, fill = color),
               color = "black", size = 5, shape = 21, stroke = 1.5) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_continuous(breaks = unique(agg_data_use$sextile),
                       labels = function(x) ifelse(x == 7, "blink", x)) +
    labs(x = "Heart rate sextile", y = "Ripple Rate", title = plot_title) +
    theme_light() +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_x_reverse(breaks = 1:6) +
    ylim(y_limits[1], y_limits[2])
}

# ============================================================================
# COLOR DEFINITIONS
# ============================================================================

grey_shade <- "#A9A9A9"

purple_colors2 <- c( #anterior hipp
  rgb(191, 163, 206, maxColorValue = 255),
  rgb(155, 128, 187, maxColorValue = 255),
  rgb(121, 95, 151, maxColorValue = 255),
  rgb(87, 63, 117, maxColorValue = 255),
  rgb(53, 33, 84, maxColorValue = 255),
  rgb(36, 16, 65, maxColorValue = 255)

)

blue_purple <- c( #hippocampus
  rgb(173, 179, 218, maxColorValue = 255),
  rgb(135, 146, 201, maxColorValue = 255),
  rgb(100, 113, 176, maxColorValue = 255),
  rgb(64, 81, 141, maxColorValue = 255),
  rgb(34, 53, 106, maxColorValue = 255),
  rgb(19, 27, 71, maxColorValue = 255)
)

pink_colors <- c( #posterior hipp
  rgb(227, 161, 199, maxColorValue = 255),
  rgb(213, 123, 178, maxColorValue = 255),
  rgb(187, 89, 160, maxColorValue = 255),
  rgb(154, 59, 129, maxColorValue = 255),
  rgb(120, 33, 100, maxColorValue = 255),
  rgb(81, 12, 64, maxColorValue = 255)
)

orange_colors <- c( #amygdala
  rgb(249, 163, 95, maxColorValue = 255),
  rgb(236, 128, 53, maxColorValue = 255),
  rgb(213, 92, 39, maxColorValue = 255),
  rgb(183, 62, 37, maxColorValue = 255),
  rgb(148, 34, 30, maxColorValue = 255),
  rgb(110, 14, 16, maxColorValue = 255)
)

# ============================================================================
# AGGREGATE ARRAYS - ALL SLEEP
# ============================================================================

subject_avg_data50 <- heart_rate_sleep_50 %>%
  group_by(sub_num, sextile, hipp_or_amy) %>%
  summarise(ripple_rate2 = mean(ripple_rate),
            sem_ripple = sd(ripple_rate, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop') %>%
  mutate(color = ifelse(sextile == 7,
                        grey_shade,
                        ifelse(hipp_or_amy == 1, blue_purple[sextile], orange_colors[sextile])))

agg_data_simple50 <- subject_avg_data50 %>%
  group_by(sextile, hipp_or_amy) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate2, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(color = ifelse(sextile == 7,
                        grey_shade,
                        ifelse(hipp_or_amy == 1, blue_purple[sextile], orange_colors[sextile])))

hipp_agg50 <- agg_data_simple50 %>% filter(hipp_or_amy == 1)
amy_agg50 <- agg_data_simple50 %>% filter(hipp_or_amy == 2)

ripple_rate_final_hipp50 <- heart_rate_sleep_50 %>% filter(hipp_or_amy == 1)
subject_avg_data_hipp50 <- ripple_rate_final_hipp50 %>%
  group_by(sub_num, sextile, A_or_P) %>%
  summarise(ripple_rate2 = mean(ripple_rate),
            sem_ripple = sd(ripple_rate, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop')

agg_data_AP50 <- subject_avg_data_hipp50 %>%
  group_by(sextile, A_or_P) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate2, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

agg_data_hipp_A50 <- agg_data_AP50 %>% filter(A_or_P == 0)
agg_data_hipp_P50 <- agg_data_AP50 %>% filter(A_or_P == 1)

# ============================================================================
# AGGREGATE ARRAYS - NREM
# ============================================================================

subject_avg_data50NREM <- heart_rate_sleep_50 %>%
  group_by(sub_num, sextile, hipp_or_amy) %>%
  summarise(ripple_rate2 = mean(ripple_rate_NREM), .groups = 'drop')

agg_data_simple50NREM <- subject_avg_data50NREM %>%
  group_by(sextile, hipp_or_amy) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate2, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(color = ifelse(sextile == 7,
                        grey_shade,
                        ifelse(hipp_or_amy == 1, blue_purple[sextile], orange_colors[sextile])))

hipp_agg50NREM <- agg_data_simple50NREM %>% filter(hipp_or_amy == 1)
amy_agg50NREM <- agg_data_simple50NREM %>% filter(hipp_or_amy == 2)

subject_avg_data_hipp50NREM <- ripple_rate_final_hipp50 %>%
  group_by(sub_num, sextile, A_or_P) %>%
  summarise(ripple_rate2 = mean(ripple_rate_NREM), .groups = 'drop')

agg_data_AP50NREM <- subject_avg_data_hipp50NREM %>%
  group_by(sextile, A_or_P) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate2, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

agg_data_hipp_A50NREM <- agg_data_AP50NREM %>% filter(A_or_P == 0)
agg_data_hipp_P50NREM <- agg_data_AP50NREM %>% filter(A_or_P == 1)

# ============================================================================
# AGGREGATE ARRAYS - WAKE
# ============================================================================

subject_avg_data50W <- heart_rate_sleep_50 %>%
  group_by(sub_num, sextile, hipp_or_amy) %>%
  summarise(ripple_rate2 = mean(ripple_rate_W), .groups = 'drop')

agg_data_simple50W <- subject_avg_data50W %>%
  group_by(sextile, hipp_or_amy) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate2, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(color = ifelse(sextile == 7,
                        grey_shade,
                        ifelse(hipp_or_amy == 1, blue_purple[sextile], orange_colors[sextile])))

hipp_agg50W <- agg_data_simple50W %>% filter(hipp_or_amy == 1)
amy_agg50W <- agg_data_simple50W %>% filter(hipp_or_amy == 2)

subject_avg_data_hipp50W <- ripple_rate_final_hipp50 %>%
  group_by(sub_num, sextile, A_or_P) %>%
  summarise(ripple_rate2 = mean(ripple_rate_W), .groups = 'drop')

agg_data_AP50W <- subject_avg_data_hipp50W %>%
  group_by(sextile, A_or_P) %>%
  summarise(
    ripple_rate_mean = mean(ripple_rate2, na.rm = TRUE),
    ripple_rate_sem = sd(ripple_rate2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

agg_data_hipp_A50W <- agg_data_AP50W %>% filter(A_or_P == 0)
agg_data_hipp_P50W <- agg_data_AP50W %>% filter(A_or_P == 1)

# ============================================================================
# RIPPLE RATE PLOTS - ALL SLEEP
# ============================================================================

hipp_final_simple5 <- simple_plot_final_nodots_indiv(
  agg_data_use = hipp_agg50,
  y_limits = c(0.02, 0.06),
  color_vector = blue_purple,
  grey_shade = grey_shade,
  plot_title = "Hippocampus"
)
hipp_final_simple5

amy_final_simple5 <- simple_plot_final_nodots_indiv(
  agg_data_use = amy_agg50,
  y_limits = c(0.02, 0.06),
  color_vector = orange_colors,
  grey_shade = grey_shade,
  plot_title = "Amygdala"
)
amy_final_simple5

hipp_final_simple_A5 <- simple_plot_final_nodots_indiv(
  agg_data_use = agg_data_hipp_A50,
  y_limits = c(0.02, 0.06),
  color_vector = purple_colors2,
  grey_shade = grey_shade,
  plot_title = "Anterior hippocampus"
)
hipp_final_simple_A5

hipp_final_simple_P5 <- simple_plot_final_nodots_indiv(
  agg_data_use = agg_data_hipp_P50,
  y_limits = c(0.02, 0.06),
  color_vector = pink_colors,
  grey_shade = grey_shade,
  plot_title = "Posterior hippocampus"
)
hipp_final_simple_P5

# ============================================================================
# RIPPLE RATE PLOTS - NREM
# ============================================================================

hipp_final_simple5_NREM <- simple_plot_final_nodots_indiv(
  agg_data_use = hipp_agg50NREM,
  y_limits = c(0.02, 0.06),
  color_vector = blue_purple,
  grey_shade = grey_shade,
  plot_title = "Hippocampus"
)
hipp_final_simple5_NREM

amy_final_simple5_NREM <- simple_plot_final_nodots_indiv(
  agg_data_use = amy_agg50NREM,
  y_limits = c(0.02, 0.07),
  color_vector = orange_colors,
  grey_shade = grey_shade,
  plot_title = "Amygdala"
)
amy_final_simple5_NREM

hipp_final_simple_A5_NREM <- simple_plot_final_nodots_indiv(
  agg_data_use = agg_data_hipp_A50NREM,
  y_limits = c(0.02, 0.07),
  color_vector = purple_colors2,
  grey_shade = grey_shade,
  plot_title = "Anterior hippocampus"
)
hipp_final_simple_A5_NREM

hipp_final_simple_P5_NREM <- simple_plot_final_nodots_indiv(
  agg_data_use = agg_data_hipp_P50NREM,
  y_limits = c(0.02, 0.07),
  color_vector = pink_colors,
  grey_shade = grey_shade,
  plot_title = "Posterior hippocampus"
)
hipp_final_simple_P5_NREM

# ============================================================================
# RIPPLE RATE PLOTS - WAKE
# ============================================================================

hipp_final_simple5_W <- simple_plot_final_nodots_indiv(
  agg_data_use = hipp_agg50W,
  y_limits = c(0.02, 0.07),
  color_vector = blue_purple,
  grey_shade = grey_shade,
  plot_title = "Hippocampus"
)
hipp_final_simple5_W

amy_final_simple5_W <- simple_plot_final_nodots_indiv(
  agg_data_use = amy_agg50W,
  y_limits = c(0.02, 0.07),
  color_vector = orange_colors,
  grey_shade = grey_shade,
  plot_title = "Amygdala"
)
amy_final_simple5_W

hipp_final_simple_A5_W <- simple_plot_final_nodots_indiv(
  agg_data_use = agg_data_hipp_A50W,
  y_limits = c(0.02, 0.07),
  color_vector = purple_colors2,
  grey_shade = grey_shade,
  plot_title = "Anterior hippocampus"
)
hipp_final_simple_A5_W

hipp_final_simple_P5_W <- simple_plot_final_nodots_indiv(
  agg_data_use = agg_data_hipp_P50W,
  y_limits = c(0.02, 0.07),
  color_vector = pink_colors,
  grey_shade = grey_shade,
  plot_title = "Posterior hippocampus"
)
hipp_final_simple_P5_W

# ============================================================================
# COMBINED RIPPLE RATE PLOTS (A/P/Amy overlay)
# ============================================================================

agg_hipp_A <- agg_data_hipp_A50 %>%
  mutate(region = "hipp_A",
         line_color = ifelse(sextile == 7, grey_shade, purple_colors2[4]),
         dot_color = ifelse(sextile == 7, "#dec9e2", purple_colors2[sextile]))
agg_hipp_P <- agg_data_hipp_P50 %>%
  mutate(region = "hipp_P",
         line_color = ifelse(sextile == 7, grey_shade, pink_colors[4]),
         dot_color = ifelse(sextile == 7, '#efcee3', pink_colors[sextile]))
agg_amy <- amy_agg50 %>%
  mutate(region = "amy",
         line_color = ifelse(sextile == 7, grey_shade, orange_colors[4]),
         dot_color = ifelse(sextile == 7, '#f9cfae', orange_colors[sextile]))

all_data_APamy <- bind_rows(agg_hipp_A, agg_hipp_P, agg_amy)

line_types_APamy <- c("hipp_A" = "solid", "hipp_P" = "dashed", "amy" = "dotted")
plot_all_sextiles_APamy <- ggplot() +
  geom_line(data = all_data_APamy %>% filter(sextile != 7),
            aes(x = sextile, y = ripple_rate_mean, group = region,
                color = line_color, linetype = region),
            linewidth = 1.5, show.legend = FALSE) +
  geom_errorbar(data = all_data_APamy,
                aes(x = sextile, ymin = ripple_rate_mean - ripple_rate_sem,
                    ymax = ripple_rate_mean + ripple_rate_sem),
                width = 0.2, size = 1, show.legend = FALSE) +
  geom_point(data = all_data_APamy,
             aes(x = sextile, y = ripple_rate_mean, fill = dot_color),
             shape = 21, size = 5, stroke = 1.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linetype_manual(values = line_types_APamy) +
  scale_x_continuous(breaks = unique(all_data_APamy$sextile),
                     labels = function(x) ifelse(x == 7, "blink", x)) +
  labs(x = "RR interval sextile", y = "Ripple Rate (Hz)",
       title = "Ripple Rate by Region and RR Sextile") +
  theme_light() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white")) +
  ylim(0.02, 0.07) +
  scale_x_reverse(breaks = 1:6)

plot_all_sextiles_APamy

# NREM version
agg_hipp_ANREM <- agg_data_hipp_A50NREM %>%
  mutate(region = "hipp_A",
         line_color = ifelse(sextile == 7, grey_shade, purple_colors2[4]),
         dot_color = ifelse(sextile == 7, "#dec9e2", purple_colors2[sextile]))
agg_hipp_PNREM <- agg_data_hipp_P50NREM %>%
  mutate(region = "hipp_P",
         line_color = ifelse(sextile == 7, grey_shade, pink_colors[4]),
         dot_color = ifelse(sextile == 7, '#efcee3', pink_colors[sextile]))
agg_amyNREM <- amy_agg50NREM %>%
  mutate(region = "amy",
         line_color = ifelse(sextile == 7, grey_shade, orange_colors[4]),
         dot_color = ifelse(sextile == 7, '#f9cfae', orange_colors[sextile]))

all_data_APamyNREM <- bind_rows(agg_hipp_ANREM, agg_hipp_PNREM, agg_amyNREM)

plot_all_sextiles_APamyNREM <- ggplot() +
  geom_line(data = all_data_APamyNREM %>% filter(sextile != 7),
            aes(x = sextile, y = ripple_rate_mean, group = region,
                color = line_color, linetype = region),
            linewidth = 1.5, show.legend = FALSE) +
  geom_errorbar(data = all_data_APamyNREM,
                aes(x = sextile, ymin = ripple_rate_mean - ripple_rate_sem,
                    ymax = ripple_rate_mean + ripple_rate_sem),
                width = 0.2, size = 1, show.legend = FALSE) +
  geom_point(data = all_data_APamyNREM,
             aes(x = sextile, y = ripple_rate_mean, fill = dot_color),
             shape = 21, size = 5, stroke = 1.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linetype_manual(values = line_types_APamy) +
  scale_x_continuous(breaks = unique(all_data_APamyNREM$sextile),
                     labels = function(x) ifelse(x == 7, "blink", x)) +
  labs(x = "RR interval sextile", y = "Ripple Rate (Hz)",
       title = "Ripple Rate by Region and RR Sextile") +
  theme_light() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white")) +
  ylim(0.02, 0.07) +
  scale_x_reverse(breaks = 1:6)

plot_all_sextiles_APamyNREM

# Wake version
agg_hipp_AWake <- agg_data_hipp_A50W %>%
  mutate(region = "hipp_A",
         line_color = ifelse(sextile == 7, grey_shade, purple_colors2[4]),
         dot_color = ifelse(sextile == 7, "#dec9e2", purple_colors2[sextile]))
agg_hipp_PWake <- agg_data_hipp_P50W %>%
  mutate(region = "hipp_P",
         line_color = ifelse(sextile == 7, grey_shade, pink_colors[4]),
         dot_color = ifelse(sextile == 7, '#efcee3', pink_colors[sextile]))
agg_amyWake <- amy_agg50W %>%
  mutate(region = "amy",
         line_color = ifelse(sextile == 7, grey_shade, orange_colors[4]),
         dot_color = ifelse(sextile == 7, '#f9cfae', orange_colors[sextile]))

all_data_APamyWake <- bind_rows(agg_hipp_AWake, agg_hipp_PWake, agg_amyWake)

plot_all_sextiles_APamyWake <- ggplot() +
  geom_line(data = all_data_APamyWake %>% filter(sextile != 7),
            aes(x = sextile, y = ripple_rate_mean, group = region,
                color = line_color, linetype = region),
            linewidth = 1.5, show.legend = FALSE) +
  geom_errorbar(data = all_data_APamyWake,
                aes(x = sextile, ymin = ripple_rate_mean - ripple_rate_sem,
                    ymax = ripple_rate_mean + ripple_rate_sem),
                width = 0.2, size = 1, show.legend = FALSE) +
  geom_point(data = all_data_APamyWake,
             aes(x = sextile, y = ripple_rate_mean, fill = dot_color),
             shape = 21, size = 5, stroke = 1.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linetype_manual(values = line_types_APamy) +
  scale_x_continuous(breaks = unique(all_data_APamyWake$sextile),
                     labels = function(x) ifelse(x == 7, "blink", x)) +
  labs(x = "RR interval sextile", y = "Ripple Rate (Hz)",
       title = "Ripple Rate by Region and RR Sextile") +
  theme_light() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white")) +
  ylim(0.02, 0.07) +
  scale_x_reverse(breaks = 1:6)

plot_all_sextiles_APamyWake

# ============================================================================
# CREATE MODULATION SCORE ARRAYS
# ============================================================================

# Array by electrode
modulation_df_elec <- heart_rate_sleep_50 %>%
  group_by(sub_num, electrode_num, RR_sextile) %>%
  summarise(
    mean_ripple_rate = mean(ripple_rate, na.rm = TRUE),
    mean_ripple_rate_NREM = mean(ripple_rate_NREM, na.rm = TRUE),
    mean_ripple_rate_W = mean(ripple_rate_W, na.rm = TRUE),
    AP_uncal = A_or_P[1],
    hipp_or_amy = hipp_or_amy[1],
    A_or_P_or_amy = A_or_P_or_amy[1],
    AP = AP[1],
    .groups = "drop")

# Calculate modulation differences by electrode
modulation_diff_df_elec <- modulation_df_elec %>%
  group_by(sub_num, electrode_num) %>%
  summarise(
    diff_mod = (mean(mean_ripple_rate[RR_sextile %in% c(1, 2)]) -
                  mean(mean_ripple_rate[RR_sextile %in% c(5, 6)])) /
      (mean(mean_ripple_rate[RR_sextile %in% c(1, 2)]) +
         mean(mean_ripple_rate[RR_sextile %in% c(5, 6)])),
    diff_mod_NREM = (mean(mean_ripple_rate_NREM[RR_sextile %in% c(1, 2)]) -
                       mean(mean_ripple_rate_NREM[RR_sextile %in% c(5, 6)])) /
      (mean(mean_ripple_rate_NREM[RR_sextile %in% c(1, 2)]) +
         mean(mean_ripple_rate_NREM[RR_sextile %in% c(5, 6)])),
    diff_mod_Wake = (mean(mean_ripple_rate_W[RR_sextile %in% c(1, 2)]) -
                       mean(mean_ripple_rate_W[RR_sextile %in% c(5, 6)])) /
      (mean(mean_ripple_rate_W[RR_sextile %in% c(1, 2)]) +
         mean(mean_ripple_rate_W[RR_sextile %in% c(5, 6)])),
    AP_uncal = AP_uncal[1],
    hipp_or_amy = hipp_or_amy[1],
    AP = AP[1],
    .groups = "drop"
  )

# Fix participant number. These are two different implants of the same participant. Non-overlapping electrodes, but same participant.
modulation_diff_df_elec$sub_num[modulation_diff_df_elec$sub_num == 49] <- 42
modulation_df_elec$sub_num[modulation_df_elec$sub_num == 49] <- 42

# Array by region
modulation_diff_df_region <- modulation_df_elec %>%
  group_by(sub_num, A_or_P_or_amy) %>%
  summarise(
    diff_mod = (mean(mean_ripple_rate[RR_sextile %in% c(1, 2)]) -
                  mean(mean_ripple_rate[RR_sextile %in% c(5, 6)])) /
      (mean(mean_ripple_rate[RR_sextile %in% c(1, 2)]) +
         mean(mean_ripple_rate[RR_sextile %in% c(5, 6)])),
    diff_mod_NREM = (mean(mean_ripple_rate_NREM[RR_sextile %in% c(1, 2)]) -
                       mean(mean_ripple_rate_NREM[RR_sextile %in% c(5, 6)])) /
      (mean(mean_ripple_rate_NREM[RR_sextile %in% c(1, 2)]) +
         mean(mean_ripple_rate_NREM[RR_sextile %in% c(5, 6)])),
    diff_mod_Wake = (mean(mean_ripple_rate_W[RR_sextile %in% c(1, 2)]) -
                       mean(mean_ripple_rate_W[RR_sextile %in% c(5, 6)])) /
      (mean(mean_ripple_rate_W[RR_sextile %in% c(1, 2)]) +
         mean(mean_ripple_rate_W[RR_sextile %in% c(5, 6)])),
    AP_uncal = AP_uncal[1],
    hipp_or_amy = hipp_or_amy[1],
    .groups = "drop"
  )

# Filter and prepare for plotting
modulation_diff_df_filtered <- modulation_diff_df_region %>%
  filter(A_or_P_or_amy %in% c(0, 1, 2)) %>%
  mutate(
    Region = factor(A_or_P_or_amy,
                    levels = c(0, 1, 2),
                    labels = c("Ant. Hipp", "Post. Hipp", "Amyg.")),
    line_color = case_when(
      A_or_P_or_amy == 0 ~ "#9a81bb",
      A_or_P_or_amy == 1 ~ "#d57bb2",
      A_or_P_or_amy == 2 ~ "#ec8035"
    )
  )

# Summary statistics
modulation_summary <- modulation_diff_df_filtered %>%
  group_by(Region, line_color) %>%
  summarise(
    mean_diff_mod = mean(diff_mod, na.rm = TRUE),
    sem_diff_mod = sd(diff_mod, na.rm = TRUE) / sqrt(n()),
    mean_diff_mod_NREM = mean(diff_mod_NREM, na.rm = TRUE),
    sem_diff_mod_NREM = sd(diff_mod_NREM, na.rm = TRUE) / sqrt(n()),
    mean_diff_mod_W = mean(diff_mod_Wake, na.rm = TRUE),
    sem_diff_mod_W = sd(diff_mod_Wake, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ============================================================================
# PLOT 1: RR modulation of ripple rate - ALL SLEEP - SUPPLEMENT
# ============================================================================

plot_diff_APamy_mod_allsleep <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  geom_line(
    data = modulation_diff_df_filtered,
    aes(x = Region, y = diff_mod, group = sub_num),
    color = "gray70",
    linewidth = 0.7,
    alpha = 0.3
  ) +
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
    shape = ifelse(modulation_summary$Region == "Ant. Hipp", 21,
                   ifelse(modulation_summary$Region == "Post. Hipp", 23, 24)),
    size = 4,
    stroke = 1.5,
    color = "black"
  ) +
  scale_fill_identity() +
  labs(
    x = NULL,
    y = "Ripple rate selectivity index",
    title = "RR modulation of ripple rate"
  ) +
  theme_light(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_reverse(limits = c(-.5, .6))
plot_diff_APamy_mod_allsleep


# Model 4: All sleep, modulation score
modulation_diff_df_elec_hipp <- modulation_diff_df_elec %>% filter(hipp_or_amy == 1)
modulation_diff_df_elec_hipp$diff_mod_rev <- 0 - modulation_diff_df_elec_hipp$diff_mod
model.df <- lmer(diff_mod_rev ~ AP + (1 | sub_num), REML = F,
                 data = modulation_diff_df_elec_hipp)
summary(model.df)

# ============================================================================
# PLOT 2: RR modulation of ripple rate - NREM ONLY - SUPPLEMENT
# ============================================================================

plot_diff_APamy_mod_NREM <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  geom_line(data = modulation_diff_df_filtered,
            aes(x = Region, y = diff_mod_NREM, group = sub_num),
            color = "gray70", linewidth = 0.7, alpha = 0.3) +
  geom_jitter(data = modulation_diff_df_filtered,
              aes(x = Region, y = diff_mod_NREM, fill = line_color),
              shape = ifelse(modulation_diff_df_filtered$Region == "Ant. Hipp", 21,
                             ifelse(modulation_diff_df_filtered$Region == "Post. Hipp", 23, 24)),
              size = 4, alpha = 0.3, width = 0, stroke = .2, show.legend = FALSE) +
  geom_errorbar(data = modulation_summary,
                aes(x = Region,
                    ymin = mean_diff_mod_NREM - sem_diff_mod_NREM,
                    ymax = mean_diff_mod_NREM + sem_diff_mod_NREM),
                width = 0.15, color = "black", size = 1) +
  geom_point(data = modulation_summary,
             aes(x = Region, y = mean_diff_mod_NREM, fill = line_color),
             shape = ifelse(modulation_summary$Region == "Ant. Hipp", 21,
                            ifelse(modulation_summary$Region == "Post. Hipp", 23, 24)),
             size = 4, stroke = 1.5, color = "black") +
  scale_fill_identity() +
  labs(x = NULL, y = "Ripple rate selectivity index",
       title = "RR modulation of ripple rate") +
  theme_light(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_reverse(limits = c(-.5, .6))
plot_diff_APamy_mod_NREM

# ============================================================================
# PLOT 3: RR modulation of ripple rate - WAKE ONLY - SUPPLEMENT
# ============================================================================

plot_diff_APamy_mod_Wake <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  geom_line(data = modulation_diff_df_filtered,
            aes(x = Region, y = diff_mod_Wake, group = sub_num),
            color = "gray70", linewidth = 0.7, alpha = 0.3) +
  geom_jitter(data = modulation_diff_df_filtered,
              aes(x = Region, y = diff_mod_Wake, fill = line_color),
              shape = ifelse(modulation_diff_df_filtered$Region == "Ant. Hipp", 21,
                             ifelse(modulation_diff_df_filtered$Region == "Post. Hipp", 23, 24)),
              size = 4, alpha = 0.3, width = 0, stroke = .2, show.legend = FALSE) +
  geom_errorbar(data = modulation_summary,
                aes(x = Region,
                    ymin = mean_diff_mod_W - sem_diff_mod_W,
                    ymax = mean_diff_mod_W + sem_diff_mod_W),
                width = 0.15, color = "black", size = 1) +
  geom_point(data = modulation_summary,
             aes(x = Region, y = mean_diff_mod_W, fill = line_color),
             shape = ifelse(modulation_summary$Region == "Ant. Hipp", 21,
                            ifelse(modulation_summary$Region == "Post. Hipp", 23, 24)),
             size = 4, stroke = 1.5, color = "black") +
  scale_fill_identity() +
  labs(x = NULL, y = "Ripple rate selectivity index",
       title = "RR modulation of ripple rate") +
  theme_light(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_reverse(limits = c(-.5, .6))
plot_diff_APamy_mod_Wake




