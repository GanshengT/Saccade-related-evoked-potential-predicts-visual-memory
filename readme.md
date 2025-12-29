# Description

This repository contains MATLAB / Python scripts for SEEG data processing and analysis for the manuscript titled ["Mind's eye: Saccade-related evoked potentials support visual encoding in humans"](https://pubmed.ncbi.nlm.nih.gov/41292659/)

---

## Structure
root/
│
├── functions/
│   ├── identify_ch_with_epileptic_discharge.m
│   ├── co_registration_seeg.m
│   └── …
│
├── data/
│   ├── example subject (BJH025)
        ├── Local field potential (LFP) timelocked to saccade onset during image presentation
        ├── LFP timelocked to saccade onset during fixation cross
        ├── EOG signals
│   ├── electrode contact information
    └── …
│
│
├── results/
│   ├── Saccade-related evoked potential (SREP) characterization (single-trial)
    ├── Averaged SREP
│
│
└── README.md

---

## Functions

### `identify_ch_with_epileptic_discharge`

**Description**  
Identifies SEEG channels exhibiting epileptic discharges 

**Requirements**
- MATLAB R2021b 
- Signal Processing Toolbox

**Input Data**
-  SEEG recordings 
---

### `co_registration_seeg_electrodes`

**Description**  
Performs CT–MRI co-registration and maps SEEG electrodes to a referenced anatomical space.

**Requirements**
- MATLAB
- FreeSurfer
- ANTs
- Gansheng Tan's imaging processing [toolbox] (https://github.com/GanshengT/intracranial_contact_loc)


## Data 


## Result

Analysis outputs are saved in `/result/`.

