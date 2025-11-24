
# ==============================================================================
# ANALYSIS_PIPELINE.R
# Purpose: Execute survival analysis and generate Figures 1, 2, and 3.
#
# This script:
# 1. Loads dependencies and local functions.
# 2. Loads and preprocesses the phenotype data.
# 3. Defines endpoint logic.
# 4. Generates and saves the primary figures.
# ==============================================================================

# --- 1. Setup ---
# setwd("/home/ivm/SCN5A") # Set your working directory here
source("Methods.R") # Load helper functions

# Ensure output directory exists
dir.create("output", showWarnings = FALSE)

# Define Covariate List (as used in Cox Models)
covlist <- "byear + SEX_IMPUTED + batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

# --- 2. Load and Preprocess Data ---
message("Loading and cleaning data...")
# Note: Ensure "df_clean_R12_new.txt" is in your WD or provide full path
df_raw <- get_clean_data(pheno_path = "df_clean_R12_new.txt")

# Prepare Final Dataframe
# (Add Competing Risks and Pacemaker logic)
df_final <- df_raw %>%
  # Example placeholder for Pacemaker (PM) logic based on provided source code
  mutate(
    PM_STATUS = if_else(!is.na(tryCatch(FIRST_EVENT_AGE, error=function(e) NA)), 1L, 0L),
    PM_TIME_TO_EVENT = coalesce(tryCatch(FIRST_EVENT_AGE, error=function(e) NA), AGE_AT_DEATH_OR_END_OF_FOLLOWUP)
  )

# Apply Competing Risk Logic (Recodes events if Death occurs first)
df_final <- add_competing_risk(df_final)

# --- 3. Figure Generation ---

#' Create Figure 1: Arrhythmia Endpoints
#' Generates a 2x2 grid for AF, Sick Sinus, Ventricular Premature Depol, and AV Block.
create_final_figure1 <- function(do_trunc = FALSE) {
  
  endpoints <- list(
    list(col="I48",  ylab="Incidence of Atrial Fibrillation (I48)"),
    list(col="I495", ylab="Incidence of Sick Sinus Syndrome (I49.5)"),
    list(col="I493", ylab="Incidence of Ventricular Premature Depol. (I49.3)"),
    list(col="I4",   ylab="Incidence of AV Block (I44.0 - I44.3)")
  )
  
  plot_list <- list()
  
  for(ep in endpoints) {
    # Censor at 80 years for Arrhythmia endpoints
    data_censored <- censor_survival_data(df_final, 
                                          paste0(ep$col, "_TIME_TO_EVENT"), 
                                          paste0(ep$col, "_STATUS"), 80) %>%
      filter(Genotype != "Homozygous") %>%
      mutate(Genotype = fct_drop(Genotype), Genotype = fct_relevel(Genotype, "Wild Type"))
    
    p <- run_and_plot_analysis(
      data = data_censored,
      model_predictors = paste0("Genotype + ", covlist),
      time_col = paste0(ep$col, "_TIME_TO_EVENT_censored"),
      status_col = paste0(ep$col, "_STATUS_censored"),
      group_col = "Genotype",
      hr_term = "GenotypeHeterozygous",
      hr_label = "Heterozygotes",
      plot_title = "",
      ylab = ep$ylab,
      legend_title = "T220 carrier status",
      left_trunc = do_trunc
    )
    plot_list[[ep$col]] <- style_ggsurvplot(p)
  }
  
  combine_survplots(
    styled_list = plot_list,
    ylabels = sapply(endpoints, function(x) x$ylab),
    file = "output/Figure1.pdf",
    width = 10, height = 10
  )
}

#' Create Figure 2: Mortality Endpoints
#' Compares Cardiac Arrhythmia Death vs All-Cause Mortality.
create_final_figure2 <- function(do_trunc = FALSE) {
  
  # Panel A: Cardiac Death (Censor at 100)
  data_ca <- censor_survival_data(df_final, "DEATH_CA_TIME_TO_EVENT", "DEATH_CA_STATUS", 100) %>%
    filter(Genotype != "Homozygous") %>%
    mutate(Genotype = fct_drop(Genotype))
  
  plotCADeath <- run_and_plot_analysis(
    data = data_ca,
    model_predictors = paste0("Genotype + ", covlist),
    time_col = "DEATH_CA_TIME_TO_EVENT_censored",
    status_col = "DEATH_CA_STATUS_censored",
    group_col = "Genotype",
    hr_term = "GenotypeHeterozygous",
    hr_label = "Heterozygotes",
    plot_title = "",
    ylab = "Incidence of death due to cardiac arrhythmias",
    legend_title = "T220 carrier status",
    left_trunc = do_trunc
  )
  
  # Panel B: All-Cause Death
  data_all <- censor_survival_data(df_final, "DEATH_TIME_TO_EVENT", "DEATH_STATUS", 100) %>%
    filter(Genotype != "Homozygous")
  
  plotDeath <- run_and_plot_analysis(
    data = data_all,
    model_predictors = paste0("Genotype + ", covlist),
    time_col = "DEATH_TIME_TO_EVENT_censored",
    status_col = "DEATH_STATUS_censored",
    group_col = "Genotype",
    hr_term = "GenotypeHeterozygous",
    hr_label = "Heterozygotes",
    plot_title = "",
    ylab = "Incidence of all-cause death",
    legend_title = "T220 carrier status",
    left_trunc = do_trunc
  )
  
  combine_survplots(
    styled_list = list(CA=style_ggsurvplot(plotCADeath), All=style_ggsurvplot(plotDeath)),
    ylabels = c("Cumulative incidence of \n death (Cardiac Arrhythmias)",
                "Cumulative incidence of \n all-cause mortality"),
    file = "output/Figure2.pdf",
    width = 9, height = 6
  )
}

#' Create Figure 3: PGS x Genotype Interaction
#' Shows incidence of AF based on Genotype and PGS strata.
create_final_figure3 <- function(do_trunc = FALSE) {
  
  data_pgs <- censor_survival_data(df_final, "I48_TIME_TO_EVENT", "I48_STATUS", 80) %>%
    filter(Genotype != "Homozygous")
  
  # Custom Dynamic Palette for PGS groups
  pal <- get_dynamic_palette(total_colors = 3, n_shades = 3, lighten_amt = 0.3, darken_amt = 0.3)[c(6,3,5,2,4,1)]
  
  plotPGS <- run_and_plot_analysis(
    data = data_pgs,
    model_predictors = paste0("Genotype + PRSnorm + ", covlist),
    time_col = "I48_TIME_TO_EVENT_censored",
    status_col = "I48_STATUS_censored",
    group_col = "pgs_group_GT_mapped",
    hr_term = "PRSnorm",
    hr_label = "Per PGS SD",
    plot_title = "",
    ylab = "Incidence of AF (I48)",
    legend_title = "T220I carrier and PGS Status",
    palette_name = pal,
    left_trunc = do_trunc
  )
  
  # Save directly (single large plot)
  combine_survplots(
    styled_list = list(PGS=style_ggsurvplot(plotPGS)),
    ylabels = "Cumulative incidence of AF (I48)",
    ncol_grid = 1,
    file = "output/Figure3.pdf",
    width = 8, height = 7
  )
}

# --- 4. Execution ---
message("Generating Figures...")
create_final_figure1()
create_final_figure2()
create_final_figure3()
message("Done. Figures saved to output/ directory.")
