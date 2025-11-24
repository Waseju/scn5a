# ==============================================================================
# EXTENDED_ANALYSIS.R
# Purpose: Advanced analysis including Competing Risks (Fine-Gray), 
#          Supplementary Truncated CIF plots, and Forest Plots.
# Dependencies: functions.R, tidycmprsk, broom, cowplot
# ==============================================================================

source("Methods.R")
library(tidycmprsk)
library(broom)

# ==============================================================================
# 1. Fine-Gray Competing Risks (sHR) Analysis
# ==============================================================================

#' Fit Fine-Gray Model via Tidycmprsk
#' 
#' @param df Data frame with censored time/status columns.
#' @param fail_label Label for the primary failure event.
#' @param time_col Name of time column (default "time_r").
#' @param status_col Name of status column (default "status").
fit_tidycmprsk <- function(df, fail_label, time_col = "time_r", status_col = "status",
                           covariates = c("Genotype", "batch", "byear", "SEX_IMPUTED", paste0("PC", 1:10))) {
  
  f_crr <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", 
                             paste(covariates, collapse = " + ")))
  
  list(
    # Cumulative Incidence Object
    ci = tidycmprsk::cuminc(as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ Genotype")), data = df),
    # Fine-Gray Model Object
    fit_crr = tidycmprsk::crr(f_crr, data = df, failcode = fail_label),
    # Tidy Table Output
    tidy = broom::tidy(
      tidycmprsk::crr(f_crr, data = df, failcode = fail_label),
      exponentiate = TRUE, 
      conf.int = TRUE
    ) %>% mutate(failcode = fail_label)
  )
}

#' Prepare Data for Competing Risks
#' Rounds times and sets status to 2 (Death) if death occurred before event.
make_comp_data <- function(df, time_col, status_col, cutoff, event_label, round_by = 1) {
  
  df_cens <- censor_survival_data(df, time_col, status_col, cutoff)
  
  status_c <- paste0(status_col, "_censored")
  time_c   <- paste0(time_col, "_censored")
  
  df_cens %>%
    mutate(
      # Check if Death occurred (1) before cutoff and event status is 0
      !!status_c := as.integer(if_else(
        DEATH_STATUS == 1 & .data[[status_c]] == 0 & DEATH_TIME_TO_EVENT <= cutoff,
        2L, .data[[status_c]]
      )),
      # Round time for computational efficiency
      time_r = pmax(1e-9, round(.data[[time_c]] / round_by) * round_by),
      status = factor(.data[[status_c]], levels = c(0, 1, 2),
                      labels = c("Censor", event_label, "DEATH_COMPETING"))
    )
}

#' Generate sHR Table
#' Runs Fine-Gray models for Cardiac Death and AF.
create_sHRs <- function(df_final) {
  
  endpoints <- tibble::tribble(
    ~time_col,                ~status_col,      ~cutoff, ~label,
    "DEATH_CA_TIME_TO_EVENT", "DEATH_CA_STATUS", 100,     "Cardiac Death",
    "I48_TIME_TO_EVENT",      "I48_STATUS",       80,     "Atrial Fibrillation"
  )
  
  results <- endpoints %>%
    mutate(
      data = pmap(list(time_col, status_col, cutoff, label),
                  ~ make_comp_data(df_final, ..1, ..2, ..3, ..4)),
      fit  = map2(data, label, ~ fit_tidycmprsk(.x, .y))
    )
  
  # Return tidy table
  shr_tbl <- results %>%
    transmute(label, tidy = map(fit, "tidy")) %>%
    tidyr::unnest(tidy)
  
  return(shr_tbl)
}

# ==============================================================================
# 2. Supplementary Figure 10: Truncated CIF Grid
# ==============================================================================

#' Create Supplementary Figure 10
#' Generates a multi-panel plot of Left-Truncated Cumulative Incidence Functions.
create_sup_figure10 <- function(df_final) {
  
  # Define Palette
  my_cols <- c(
    "Heterozygous" = RColorBrewer::brewer.pal(3,"Dark2")[2],
    "Wild Type"    = RColorBrewer::brewer.pal(3,"Dark2")[1]
  )
  
  # Helper to prepare truncated data
  prep_trunc <- function(df, t_col, s_col, max_age) {
    censor_survival_data(df, t_col, s_col, max_age) %>%
      filter(Genotype != "Homozygous", BL_AGE <= max_age, BL_AGE <= .data[[paste0(t_col, "_censored")]])
  }
  
  # 1. AFib
  p1 <- plot_truncated_cif(
    prep_trunc(df_final, "I48_TIME_TO_EVENT", "I48_STATUS", 80),
    "BL_AGE", "I48_TIME_TO_EVENT_censored", "I48_STATUS_censored", "Genotype",
    "Incidence of AF (I48)", palette = my_cols
  )
  
  # 2. Sick Sinus
  p2 <- plot_truncated_cif(
    prep_trunc(df_final, "I495_TIME_TO_EVENT", "I495_STATUS", 80),
    "BL_AGE", "I495_TIME_TO_EVENT_censored", "I495_STATUS_censored", "Genotype",
    "Incidence of SSS (I49.5)", palette = my_cols
  )
  
  # 3. Ventricular Premature Depol
  p3 <- plot_truncated_cif(
    prep_trunc(df_final, "I493_TIME_TO_EVENT", "I493_STATUS", 80),
    "BL_AGE", "I493_TIME_TO_EVENT_censored", "I493_STATUS_censored", "Genotype",
    "Incidence of VPD (I49.3)", palette = my_cols
  )
  
  # 4. AV Block
  p4 <- plot_truncated_cif(
    prep_trunc(df_final, "I4_TIME_TO_EVENT", "I4_STATUS", 80),
    "BL_AGE", "I4_TIME_TO_EVENT_censored", "I4_STATUS_censored", "Genotype",
    "Incidence of AV Block", palette = my_cols
  )
  
  # 5. Cardiac Death
  p5 <- plot_truncated_cif(
    prep_trunc(df_final, "DEATH_CA_TIME_TO_EVENT", "DEATH_CA_STATUS", 100),
    "BL_AGE", "DEATH_CA_TIME_TO_EVENT_censored", "DEATH_CA_STATUS_censored", "Genotype",
    "Incidence of Cardiac Death", palette = my_cols
  )
  
  # 6. Pacemaker (PM)
  p6 <- plot_truncated_cif(
    prep_trunc(df_final, "PM_TIME_TO_EVENT", "PM_STATUS", 80),
    "BL_AGE", "PM_TIME_TO_EVENT_censored", "PM_STATUS_censored", "Genotype",
    "Incidence of Pacemaker", palette = my_cols
  )
  
  # Combine Grid
  grid_main <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
  
  # Add Legend
  shared_legend <- cowplot::get_legend(p2 + theme(legend.position = "bottom"))
  
  final_plot <- cowplot::plot_grid(grid_main, shared_legend, ncol = 1, rel_heights = c(0.95, 0.05))
  
  ggsave("output/Supp_Figure10_Truncated.pdf", final_plot, width = 15, height = 22.5)
}

# ==============================================================================
# 3. Forest Plot Utilities (PheWAS)
# ==============================================================================

#' Generate Forest Plot
#' Creates a forest plot from model results.
generate_forest_plot <- function(data, estimate_col="OR", study_col="name", 
                                 lower_ci_col="lower_ci", upper_ci_col="upper_ci", 
                                 plot_color="black", xlab="Odds Ratio") {
  
  ggplot(data, aes(x = .data[[estimate_col]], y = fct_rev(as.factor(.data[[study_col]])))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = .data[[lower_ci_col]], xmax = .data[[upper_ci_col]]), 
                   color = plot_color, height = 0.2, linewidth = 0.8) +
    geom_point(color = plot_color, shape = 15, size = 3) +
    labs(x = xlab, y = NULL) + 
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.y = element_text(color = "black")
    )
}
