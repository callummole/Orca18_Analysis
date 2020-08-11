Analysis repository for Silent Failures in Automation. Experiment Repository: https://github.com/callummole/SP_18-19

Preprint: Mole et al., 2020, Predicting Takeover Response to silent automated failures.

Data information about raw and processed data can be found in the Data folder README.

The key analysis files are in the folder _manuscript_analysis_. Also in this folder are the saved model fits.

The folder _Processing_ contains code for eye-tracking, most of which has not been applied to the current dataset so is not worth looking at. The folder _Post-processing_ contains a range of part-baked scripts of different analysis adventures - also not worth looking at.  

The analysis pipeline from the Raw_Data folder hosted on the OSF (https://osf.io/aw8kp/) to the manuscript figures is as follows:

- Extract the Raw_Data folder into the local repo _Data_ folder.

To generate the condition onset times and steering angle biases, and generate the simulated time-to-line-crossings (all relevant csvs are also found in the Raw Data folder so you can skip these steps), do the following:

1) run TrackSimulation.py to output 'simulated_roadcrossing.csv'.
2) run TrackSimulation_sobol.py to obtain the relevant steering angle biases for the random conditions. This file outputs 'SimResults_samplesobol_onsettimes.csv'.
3) run plot_failures_perspective.py to plot all the failures and output 'simulated_ttlcs.csv', used in the analysis.

To run the analysis and output figures found in the manuscript, do the following:

1) run processing_steering_only.py to output a collated steering csv in the _Data_ folder.
2) _optional_, run save_as_rds.R to save time for loading csvs into R.
3) run manuscript_figures_ttlc.rmd for TTLC results (put REFIT = TRUE to refit model).
4) run manuscript_figures_swa.rmd for steering results (put REFIT = TRUE to refit model).
5) run cogtask_performance.rmd for the cognitive task analysis.

