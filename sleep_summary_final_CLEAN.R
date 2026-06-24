# Sleep Stage Analysis and Visualization
# This script loads sleep data, filters for quality records, and creates 
# stacked bar plots showing sleep stage duration across participants
# last edited 3/24/2025
# Load required libraries
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
# SECTION 1: LOAD AND PREPARE DATA
# ============================================================================

# Load raw sleep data from CSV file
sleep_array <- read_csv("data/stage_table_9242025.csv", col_names = TRUE)

# Rename columns to meaningful variable names for consistency
colnames(sleep_array) <- c(
  "sub_num",        # Participant ID
  "session",        # Session number
  "time_W",         # Wake time (minutes)
  "time_R",         # REM sleep time (minutes)
  "time_N1",        # Stage N1 sleep time (minutes)
  "time_N2",        # Stage N2 sleep time (minutes)
  "time_N3",        # Stage N3 (deep) sleep time (minutes)
  "time_sleep",     # Total sleep time (minutes)
  "time_NREM",      # Total NREM sleep time (minutes)
  "time_other",     # Other time (minutes)
  "frag_idx",       # Sleep fragmentation index
  "frag_2",         # Alternative fragmentation metric
  "stage_conf"      # Sleep stage confidence score
)

# ============================================================================
# SECTION 2: FILTER DATA - REMOVE PROBLEMATIC RECORDS
# ============================================================================

# Filter out participants and sessions with data quality issues
# These exclusions are based on quality control checks
final_sleep_array <- sleep_array %>%
  mutate(
    sub_num = as.character(sub_num),
    session = as.character(session)
  ) %>%
  dplyr::filter(
    sub_num != 'PAV056',                                      # Exclude participant
    !(sub_num == 'PAV044' & session == 'session1'),          # Exclude specific session
    !(sub_num == 'PAV049' & session == 'session1'),
    !(sub_num == 'PAV054' & session == 'session3'),
    !(sub_num == 'PAV055' & session == 'session3'),
    !(sub_num == 'PAV057' & session %in% c('session3', 'session4')),
    !(sub_num == 'PAV064' & session == 'session1'),
    !(sub_num == 'PAV069' & session == 'session2')
  )

# ============================================================================
# SECTION 3: CALCULATE DESCRIPTIVE STATISTICS
# ============================================================================

# Calculate summary statistics for key sleep variables
# Results are in minutes (divide by 60 to convert to hours)

# Mean total sleep time in hours
mean_sleep_hours <- mean(final_sleep_array$time_sleep) / 60  # ~7.5 hours

# Standard deviation of total sleep time
sd_sleep_hours <- sd(final_sleep_array$time_sleep) / 60

# Minimum values for quality checking
min_other_hours <- min(final_sleep_array$time_other) / 60
min_NREM_minutes <- min(final_sleep_array$time_NREM)
min_wake_minutes <- min(final_sleep_array$time_W)

# ============================================================================
# SECTION 4: CREATE SUMMARY TABLE
# ============================================================================

# Generate comprehensive summary statistics for all time variables
summary_table <- final_sleep_array %>%
  select(starts_with("time_")) %>%
  pivot_longer(everything(), 
               names_to = "time_variable", 
               values_to = "value") %>%
  group_by(time_variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    std = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# Optional: view and save the summary table
# View summary_table using: print(summary_table)
# Save to file using: write.csv(summary_table, "tables/summary_table.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: CREATE SUBJECT-TO-PAPER CODE MAPPING
# ============================================================================

# Create lookup table mapping EMU participant codes to paper codes
subject_mapping <- tibble(
  sub_num = c("PAV037", "PAV038", "PAV040", "PAV042", "PAV043", "PAV044", 
              "PAV045", "PAV048", "PAV049", "PAV051", "PAV053", "PAV054", 
              "PAV055", "PAV057", "PAV058", "PAV060", "PAV062", "PAV064", 
              "PAV065", "PAV067", "PAV069"),
  paper_code = c("1", "2", "3", "4-1", "5", "6", 
                 "7", "8", "4-2", "9", "10", "11", 
                 "12", "13", "14", "15", "16", "17", 
                 "18", "19", "20")
)
# ============================================================================
# SECTION 6: PREPARE DATA FOR VISUALIZATION
# ============================================================================

# Create a row index and calculate total sleep for sorting
final_sleep_array <- final_sleep_array %>%
  left_join(subject_mapping, by = "sub_num") %>%  # Add paper codes
  mutate(
    row_index = row_number(),
    total_sleep = time_N1 + time_N2 + time_N3 + time_R + time_W,
    # Re-order by total sleep time and create factor for plotting
    row_index = row_number()
  ) %>%
  arrange(sub_num, session) %>%
  mutate(
    sub_label = paper_code,                         # Display label for x-axis (paper code)
    row_index = row_number(),                       # Numeric position in plot
    row_index = factor(row_index, levels = row_index) # Convert to factor for discrete axis
  )

# ============================================================================
# SECTION 7: DEFINE COLOR PALETTES
# ============================================================================

# Color palette for 3-stage visualization (NREM, REM, Wake)
# Uses RGB values normalized to 0-255 range
sleep_stage_colors_simple <- c(
  "time_NREM" = rgb(32, 120, 180, maxColorValue = 255),    # Blue (N1-N3 combined)
  "time_R"    = rgb(104, 193, 176, maxColorValue = 255),   # Teal (REM)
  "time_W"    = rgb(191, 191, 196, maxColorValue = 255)    # Gray (Wake)
)

# Color palette for 5-stage visualization (individual NREM stages, REM, Wake)
# Provides more detailed breakdown of sleep architecture
sleep_stage_colors_detailed <- c(
  "time_N1" = rgb(70, 118, 224, maxColorValue = 255),     # Light blue (N1)
  "time_N2" = rgb(36, 68, 142, maxColorValue = 255),      # Medium blue (N2)
  "time_N3" = rgb(8, 29, 88, maxColorValue = 255),        # Dark blue (N3 - deep sleep)
  "time_R"  = rgb(104, 193, 176, maxColorValue = 255),    # Teal (REM)
  "time_W"  = rgb(191, 191, 196, maxColorValue = 255)     # Gray (Wake)
)

# ============================================================================
# SECTION 8: PREPARE DATA FOR PLOTTING
# ============================================================================

# Create long-format dataset for 3-stage plot (NREM, REM, Wake)
sleep_long_simple <- final_sleep_array %>%
  pivot_longer(cols = c(time_NREM, time_R, time_W),
               names_to = "stage_var",
               values_to = "time_asleep")

# Create long-format dataset for 5-stage plot (individual NREM stages, REM, Wake)
sleep_long_detailed <- final_sleep_array %>%
  pivot_longer(cols = c(time_N3, time_N2, time_N1, time_R, time_W),
               names_to = "stage_var",
               values_to = "time_asleep")

# ============================================================================
# SECTION 9: CREATE VISUALIZATIONS
# ============================================================================

# Plot 1: Simple 3-stage sleep visualization (NREM combined, REM, Wake)
plot_simple <- ggplot(sleep_long_simple, 
                      aes(x = row_index, y = time_asleep / 60, fill = stage_var)) +
  geom_col() +
  scale_fill_manual(values = sleep_stage_colors_simple,
                    labels = c("time_NREM" = "NREM", 
                               "time_R" = "REM", 
                               "time_W" = "Wake")) +
  scale_x_discrete(labels = final_sleep_array$sub_label) +
  labs(title = "Sleep Stage Duration - Combined NREM",
       x = "Sleep Session (Participant)",
       y = "Time (hours)",
       fill = "Sleep Stage") +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Detailed 5-stage sleep visualization (N1, N2, N3, REM, Wake)
# This provides a more complete view of sleep architecture
plot_detailed <- ggplot(sleep_long_detailed, 
                        aes(x = row_index, y = time_asleep / 60, fill = stage_var)) +
  geom_col() +
  scale_fill_manual(values = sleep_stage_colors_detailed,
                    labels = c("time_N1" = "N1", 
                               "time_N2" = "N2", 
                               "time_N3" = "N3", 
                               "time_R" = "REM", 
                               "time_W" = "Wake")) +
  scale_x_discrete(labels = final_sleep_array$sub_label) +
  labs(title = "Sleep Stage Duration - Individual NREM Stages",
       x = "Sleep Session (Participant)",
       y = "Time (hours)",
       fill = "Sleep Stage") +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
print(plot_simple)
print(plot_detailed)
