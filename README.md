# pupil_ripple
Data and code for main figs and analyses in pupil ripple paper

sleep_ripple_plots_final_CLEAN.R creates main figures and runs mixed-effects models for sleep-ripple rate analyses
This piece of code relies on datafiles:
- all_chan_ripple_final_9242025.csv (main data file)
- sleep_electrode_array_updated_final_9242025.csv (electrode information)
- sub_elec_positions_allinfo_2025.csv (additional electrode information)

sleep_summary_final_CLEAN.R generates sleep duration information across participants, and associated supplemental figures
This piece of code relies on datafiles:
- stage_table_9242025.csv (stage duration data)

pupil_sextiles_final_CLEANED.R creates main figures and mixed-effects models for pupil-ripple rate analyses
Relies on:
- pupil_array_1_3.csv (main data)
- sub_elec_positions_allinfo_2025.csv (electrode information)
- sleep_electrode_array_updated_final_9242025.csv (electrode information)
