# Leveraging a Genetic Proxy to Investigate the Effects of Lifelong Cardiac Sodium Channel Blockade

This repository contains the R code and analysis pipeline used for the study **"Leveraging a Genetic Proxy to Investigate the Effects of Lifelong Cardiac Sodium Channel Blockade"**.

## Overview

In this study, we investigated the Finnish-enriched *SCN5A* missense variant **T220I (rs45620037)** across three large cohorts: **FinnGen**, **UK Biobank**, and **Health 2000**. We utilized survival analysis, competing risk models, and polygenic risk score (PRS) interactions to demonstrate that T220I acts as a genetic proxy for lifelong cardiac sodium channel blockade, conferring protection against atrial fibrillation while increasing the risk for specific bradyarrhythmias.

## Citation

If you use the code or insights from this repository, please cite our paper:

## Repository Structure

The analysis is split across three main R scripts:

  * **`ana.R`**: The **main execution script**. It loads the data, calls the cleaning functions, and systematically generates the Main Figures (1, 2, 3) and Supplementary Figures found in the manuscript.
  * **`new_methods.R`**: Contains the core statistical and plotting functions, including:
      * `get_clean_data()`: Data preprocessing and feature engineering.
      * `run_and_plot_analysis()`: Wrapper for Cox Proportional Hazards models.
      * `fit_crr_model()` & `plot_truncated_cif()`: Competing risk analyses and cumulative incidence functions.
      * `create_figure_X()`: Specific functions to compose the multi-panel figures.
  * **`Methods.R`**: Contains utility functions, data loading helpers, and legacy code for PheWAS and basic Kaplan-Meier plotting.

## Data Availability

**Important:** The individual-level genotype and phenotype data used in this study (FinnGen, UK Biobank, Health 2000) are **not** public.

  * **FinnGen:** Data may be accessed through Finnish Biobanks’ FinBB portal (www.finbb.fi).
  * **UK Biobank:** Data is available via application at [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/).
  * **Health 2000:** Access procedures are described [here](https://thl.fi/en/research-and-development/thl-biobank/for-researchers/application-process).

The code in this repository is provided to ensure transparency of the methods and to allow reproduction by researchers with approved access to these datasets.

## Dependencies

The analysis was performed using **R version 4.3.0**. To run the scripts, you will need the following R packages:

```r
install.packages(c(
  "tidyverse",   # includes dplyr, ggplot2, purrr, tibble, forcats, readr, tidyr
  "survival",
  "survminer",
  "cmprsk",
  "tidycmprsk",
  "fastcmprsk",  # For fast competing risks
  "cowplot",     # For figure arrangement
  "scales",
  "broom",
  "data.table",
  "prodlim",
  "gridExtra",
  "RColorBrewer",
  "rlang"
))
```

## Usage

1.  **Setup Environment:**
    Ensure all dependencies are installed. Note that `new_ana.R` contains hardcoded paths (e.g., `/home/ivm/SCN5A`) specific to the analysis environment. You must update `wdPath` and source paths to match your local directory structure.

2.  **Data Preparation:**

3.  **Running the Analysis:**
    Execute the `new_ana.R` script. 


## Statistical Methods

The code implements the following statistical approaches detailed in the Methods section of the paper:

  * **Cox Proportional Hazards:** Used for time-to-event analysis of incident disease and mortality.
  * **Competing Risks:** Fine-Gray subdistribution hazard models were used to account for the competing risk of death when analyzing non-fatal cardiac events.
  * **Left Truncation:** Survival analyses account for delayed entry (age at recruitment/DNA donation).
  * **PRS Interaction:** Interaction models between the *SCN5A* genotype and atrial fibrillation Polygenic Risk Scores (PRS-CS).

## Contact

For questions regarding the code or analysis, please contact:

  * **Julian S. Wanner**: julian.wanner@hpi.de
  * **Henrike O. Heyne**: henrike.heyne@hpi.de
