# Description

This repository contains MATLAB and Python scripts for SEEG data processing and analysis for the manuscript:
["Mind's eye: Saccade-related evoked potentials support visual encoding in humans"](https://pubmed.ncbi.nlm.nih.gov/41292659/).

## Folder overview

- `prep_oscillation_gaze.m`: main pipeline for SREP extraction, trial clustering, oscillation features, and gaze/behavior summaries.The following were included in `prep_oscillation_gaze.m` and correspond to Figure2 A-E: 
    - Description and illustration workflow for calculating SREP.
    - Validation of hierarchical clustering of trial responses.
    - Quantification of distribution separation between clustered responses (including cluster separativity metrics and ROC/AUC-based comparisons).
- `prep_erp_gaze_control.m`: control analyses for evoked responses.
- `identify_ch_with_epileptic_discharge.m`: detection and flagging of channels with epileptic discharges.
- `co_registration_seeg_electrodes.m`: CT-MRI co-registration and anatomical localization of SEEG contacts.
- `viz_seeg_coverage.m`: MATLAB visualization for anatomical contact coverage.
- `Fig1_coverage.ipynb`: Figure 1 coverage visualization workflow.
- `Fig1_saccade_behavior.ipynb`: Figure 1 saccade-behavior analysis/visualization workflow.
- `data/`: input data (including subject/session-level SEEG, gaze, and metadata files).
- `result/`: saved analysis outputs (tables and MAT files).
- `freesurfer_file/`: FreeSurfer segmentation resources used by visualization scripts.
- `CircStat2012a/`: circular statistics toolbox dependency.




## Key outputs

- `result/*_gaze_properties.csv`: gaze-event level summaries.
- `result/responsive_erp_*.mat`: channel-level SREP and clustering metrics.
- `result/df_evoked_potential_*.mat`: trial-level evoked and oscillation features.
