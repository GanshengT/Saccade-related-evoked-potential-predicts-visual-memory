# Description

This repository contains MATLAB pipelines, figure notebooks, and summary tables for the manuscript:
["Mind's eye: Saccade-related evoked potentials support visual encoding in humans"](https://pubmed.ncbi.nlm.nih.gov/41292659/).

## Repository structure

### MATLAB pipelines

- `prep_oscillation_gaze.m`: main analysis pipeline for SREP extraction, trial clustering, oscillation features, and gaze/behavior summaries.
- `prep_erp_gaze_control.m`: control analyses for evoked responses.
- `identify_ch_with_epileptic_discharge.m`: detection and flagging of channels with epileptic discharges.
- `co_registration_seeg_electrodes.m`: CT-MRI co-registration and anatomical localization of SEEG contacts.
- `viz_seeg_coverage.m`: visualization of anatomical SEEG coverage.

### Figure notebooks

- `Fig1_coverage.ipynb`
- `Fig1_saccade_behavior.ipynb`
- `Fig2_bimodal_distribution.ipynb`
- `Fig2_consistent_SREP.ipynb`
- `Fig3_SREP_present_in_fixation_cross_presentation.ipynb`
- `Fig3_srep_is_not_EMG.ipynb`
- `Fig4_polarity_oscillation_phase.ipynb`
- `Fig4_polarity_saccade_dir.ipynb`
- `Fig5_overall_spatial_dist_SREP_onset_latency.ipynb`
- `Fig5_spatial_dist_SREP_onset_latency.ipynb`
- `Fig5_spatial_dist_SREP_peak_amp.ipynb`
- `Fig6_SRNDs_predict_memory.ipynb`
- `Fig6_rf_cv_learns_the_same_mapping.ipynb`

### Data and outputs

- `data/`: tracked summary tables used by the notebooks, including:
  - `sum_responsive_erp.csv`
  - `sum_df_evoked_potential.csv`
  - `sum_df_evoked_potential_up_to_3OCs_before_SREP.csv`
  - `sum_control_responsive_erp.csv`
  - `sum_control_df_evoked_potential.csv`
  - `sum_subclassification_ks.csv`
  - `all_subj_task_gaze_properties*.csv`
  - `rf_performance_*.csv`
  - `logreg_performance_saccade_based.csv`
  - `norm_contact_data_table.csv`
  - `cv_similarity/`
- `result/`: saved subject-level outputs currently present in this checkout:
  - `BJH025_BLAES_study_twosource_gaze_properties.csv`
  - `BJH025_epileptic_activity_results_all_trials.csv`
  - `control_df_evoked_potential_BJH025_BLAES_study_twosource.csv`
  - `control_responsive_erp_BJH025_BLAES_study_twosource.csv`
- `freesurfer_file/`: FreeSurfer segmentation resources used by the coverage scripts.
- `CircStat2012a/`: CircStat dependency used by the MATLAB circular-statistics analyses.

## Download and setup

### Clone the tracked repository files

Some tracked notebooks, `.mat` files, and `data/*.csv` files are stored with Git LFS.

```bash
git clone https://github.com/GanshengT/Saccade-related-evoked-potential-predicts-visual-memory.git
cd Saccade-related-evoked-potential-predicts-visual-memory
git lfs install
git lfs pull
```

### Download `data/BJH025`

`data/BJH025` is not a tracked GitHub folder in this repository. The `.mat` files in that directory are ignored by git, so they will not appear when you clone the repository on GitHub.

To download the raw/preprocessed data files used by the MATLAB pipelines, go to (https://github.com/GanshengT/Saccade-related-evoked-potential-predicts-visual-memory/commits/main/), and search the commit on Feb 17, 2026, for example, "extended data - oscillation before SREP"

`data/BJH025` contains 772 `.mat` files and occupies > 2 GB. The expected filenames include session-level files such as:

- `BJH025_BLAES_study_twosource_session1_prep_saccade_event.mat`
- `BJH025_BLAES_study_twosource_session1_saccade_table.mat`
- `BJH025_BLAES_study_twosource_session1_channel*_saccade_onset_not_image_onset_prep_signal.mat`


## Notes

- The tracked summary table `data/sum_df_evoked_potential_up_to_3OCs_before_SREP.csv` contains the extended oscillation features before SREP.
- Additional local outputs such as `data/cv_similarity/perm_*`, are not part of the tracked repository state in this checkout.
