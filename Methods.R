# ==============================================================================
# Methods.R
# Purpose: Core utility functions for Survival Analysis, Data Cleaning, and Plotting.
#
# This script contains reusable functions for:
# 1. Data Cleaning (get_clean_data, add_competing_risk, censor_survival_data)
# 2. Survival Analysis Wrappers (run_and_plot_analysis, format_hr_text)
# 3. Plotting Utilities (style_ggsurvplot, combine_survplots, plot_truncated_cif)
# ==============================================================================

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)
library(tidyr)
library(forcats)
library(stringr)
library(prodlim) 

# ==============================================================================
# 1. Data Cleaning and Preprocessing
# ==============================================================================

#' Read and Clean Phenotype Data
#'
#' Reads the raw phenotype file, filters out problematic batches (zero events),
#' and generates mapped columns for PGS x Genotype interactions.
#'
#' @param pheno_path Character string. Path to the raw phenotype text file.
#' @return A processed data frame (tibble).
get_clean_data <- function(pheno_path) {
  
  # --- 1. Load Data ---
  if (!file.exists(pheno_path)) stop("Phenotype file not found: ", pheno_path)
  df <- data.table::fread(pheno_path) 
  
  # --- 2. Filter Batches (Optional: removes batches with 0 competing events) ---
  if ("DEATH_COMPETING_STATUS" %in% names(df)) {
    problem_batches <- df %>%
      dplyr::count(batch) %>%
      dplyr::right_join(dplyr::distinct(df, batch), by = "batch") %>%
      dplyr::filter(is.na(n)) %>%
      dplyr::pull(batch)
    
    if (length(problem_batches) > 0) {
      message("Info: Removing batches with zero events: ", paste(problem_batches, collapse = ", "))
      # Uncomment to apply filtering:
      # df <- df %>% filter(!batch %in% problem_batches)
    }
  }

  # --- 3. Define PGS x Genotype Mapping ---
  # Explicit ordering for plotting levels
  level_order <- c(
    "Non-carriers (High PGS)", "T220I carriers (High PGS)",
    "Non-carriers (Medium PGS)", "T220I carriers (Medium PGS)",
    "Non-carriers (Low PGS)", "T220I carriers (Low PGS)",
    "T220I homozygous (High PGS)", "T220I homozygous (Low PGS)","T220I homozygous (Medium PGS)",
    "Other"
  )
  
  # Map raw "prs_group + Genotype" strings to clean labels
  pgs_gt_map <- c(
    "High Wild Type"        = "Non-carriers (High PGS)",
    "Medium Wild Type"      = "Non-carriers (Medium PGS)",
    "Low Wild Type"         = "Non-carriers (Low PGS)",
    "High Heterozygous"     = "T220I carriers (High PGS)",
    "Medium Heterozygous"   = "T220I carriers (Medium PGS)",
    "Low Heterozygous"      = "T220I carriers (Low PGS)",
    "High Homozygous"       = "T220I homozygous (High PGS)",
    "Medium Homozygous"     = "T220I homozygous (Medium PGS)",
    "Low Homozygous"        = "T220I homozygous (Low PGS)"
  )
  
  # --- 4. Apply Transformations ---
  df_clean <- df %>%
    mutate(
      Genotype = factor(Genotype, levels = c("Wild Type", "Heterozygous", "Homozygous")),
      SEX_IMPUTED = factor(SEX_IMPUTED),
      batch = factor(batch)
    ) %>%
    tidyr::drop_na(batch) %>%
    dplyr::mutate(
      pgs_group_GT_mapped = dplyr::recode(
        paste(prs_group, Genotype),
        !!!pgs_gt_map,
        .default = "Other"
      ),
      pgs_group_GT_mapped = factor(pgs_group_GT_mapped, levels = level_order)
    ) %>%
    tidyr::drop_na(batch, PC1) # Ensure PC1 exists for adjustment
  
  return(df_clean)
}

#' Add Competing Risk Logic
#'
#' Modifies status columns to reflect competing risks (death). 
#' If a patient dies before the event of interest, the status is recoded to 2.
#'
#' @param data Data frame containing survival data.
#' @param death_status_col Name of the column indicating death status (0/1).
#' @param death_time_col Name of the column indicating time to death.
#' @param status_suffix Suffix used to identify status columns (default "_STATUS").
#' @param time_suffix Suffix used to identify time columns (default "_TIME_TO_EVENT").
#' @param tie_rule Character. "death_wins" (death <= event) or "event_wins" (death < event).
#' @return Data frame with updated status columns (0=Censor, 1=Event, 2=Competing Death).
add_competing_risk <- function(data,
                               death_status_col = "DEATH_STATUS",
                               death_time_col   = "DEATH_TIME_TO_EVENT",
                               status_suffix    = "_STATUS",
                               time_suffix      = "_TIME_TO_EVENT",
                               tie_rule         = "death_wins") {
  
  # Identify all status columns matching the suffix
  status_cols <- grep(paste0(status_suffix, "$"), names(data), value = TRUE)
  status_cols <- setdiff(status_cols, death_status_col) # Exclude the death column itself
  
  for (s_col in status_cols) {
    t_col <- sub(paste0(status_suffix, "$"), time_suffix, s_col)
    
    if (t_col %in% names(data)) {
      # Vectors for comparison
      ev_stat <- as.numeric(data[[s_col]])
      ev_time <- data[[t_col]]
      d_stat  <- data[[death_status_col]]
      d_time  <- data[[death_time_col]]
      
      # Determine if death happened before or at the same time as the event
      cmp <- if (tie_rule == "death_wins") (d_time <= ev_time) else (d_time < ev_time)
      
      # Handle cases where event time is missing but death occurred (assume death precluded event)
      cmp[is.na(ev_time)] <- TRUE 
      
      # Identify rows to recode: Death occurred (1), Death time known, and Death <= Event
      is_competing <- (d_stat == 1) & !is.na(d_time) & cmp
      
      # Apply recode: Change status to 2, update time to death time
      to_change <- is_competing & (is.na(ev_stat) | ev_stat != 1 | (ev_stat == 1 & cmp))
      
      ev_stat[to_change] <- 2
      data[[s_col]] <- ev_stat
      data[[t_col]][to_change] <- d_time[to_change]
    }
  }
  return(data)
}

#' Censor Survival Data
#'
#' Caps follow-up time at a specific threshold (e.g., age 80).
#'
#' @param df Data frame.
#' @param time_col Name of the time-to-event column.
#' @param status_col Name of the status column.
#' @param max_time Numeric. The cut-off time (e.g., 80 or 100).
#' @return Data frame with new `_censored` columns.
censor_survival_data <- function(df, time_col, status_col, max_time) {
  time_censored_col <- paste0(time_col, "_censored")
  status_censored_col <- paste0(status_col, "_censored")
  
  df %>%
    dplyr::mutate(
      !!time_censored_col := ifelse(.data[[time_col]] > max_time, max_time, .data[[time_col]]),
      !!status_censored_col := ifelse(.data[[time_col]] > max_time, 0, .data[[status_col]])
    )
}

# ==============================================================================
# 2. Plotting Utilities
# ==============================================================================

#' Format HR Text for Plots
#'
#' Extracts Hazard Ratio, 95% CI, and P-value from a Cox model object.
#'
#' @param model A coxph object.
#' @param term_name The name of the term in the model (e.g., "GenotypeHeterozygous").
#' @param label The label to display (e.g., "Heterozygotes").
#' @return A formatted string string.
format_hr_text <- function(model, term_name, label) {
  summ <- summary(model)
  coefs <- summ$coefficients
  conf <- summ$conf.int
  
  if (term_name %in% rownames(coefs)) {
    hr <- conf[term_name, "exp(coef)"]
    lower <- conf[term_name, "lower .95"]
    upper <- conf[term_name, "upper .95"]
    pval <- coefs[term_name, "Pr(>|z|)"]
    
    pval_fmt <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3)))
    
    return(sprintf("%s HR: %.2f (%.2f-%.2f)\n%s", label, hr, lower, upper, pval_fmt))
  } else {
    return("")
  }
}

#' Generate Dynamic Color Palette
#'
#' Creates a gradient palette based on RColorBrewer base colors.
#' Useful for plots with sub-groupings (e.g., PGS + Genotype).
get_dynamic_palette <- function(total_colors = 3, n_shades = 3, lighten_amt=0.2, darken_amt=0.2) {
  requireNamespace("colorspace")
  requireNamespace("RColorBrewer")
  
  base_colors <- RColorBrewer::brewer.pal(max(3, total_colors), "Set2")[1:total_colors]
  
  dynamic_palette <- lapply(base_colors, function(col) {
    dark_col  <- colorspace::darken(col, darken_amt)
    light_col <- colorspace::lighten(col, lighten_amt)
    ramp_func <- colorRampPalette(c(dark_col, col, light_col))
    ramp_func(n_shades)
  })
  
  rev(unlist(dynamic_palette))
}

#' Plot Truncated Cumulative Incidence (CIF)
#'
#' Uses `prodlim` to calculate and plot CIFs with specific time truncation.
#'
#' @param df Data frame.
#' @param entry_age_col Column for entry age.
#' @param exit_age_col Column for exit age.
#' @param status_col Column for status.
#' @param group_col Column for grouping factor.
#' @param ylab Y-axis label.
#' @param palette Optional color palette.
#' @return A ggplot object.
plot_truncated_cif <- function(df, entry_age_col, exit_age_col, status_col, group_col, 
                               ylab, palette = NULL) {
  
  # Define Formula: Hist(entry, exit, status) ~ Group
  fmla <- as.formula(paste0("Hist(entry=", entry_age_col, ", time=", exit_age_col, 
                            ", event=", status_col, ") ~ ", group_col))
  
  fit <- prodlim::prodlim(fmla, data = df)
  
  # Extract data for ggplot (cause = 1 is the primary event)
  cif_df <- as.data.frame(prodlim::as.data.frame.prodlim(fit, cause = 1))
  
  p <- ggplot(cif_df, aes(x = time, y = absolute_risk, colour = get(group_col))) +
    geom_step(linewidth = 1) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic(base_size = 12) +
    labs(x = "Age (years)", y = ylab, color = group_col) +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  
  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }
  
  return(p)
}

# ==============================================================================
# 3. Wrapper: Run Analysis & Plot
# ==============================================================================

#' Run Cox Analysis and Generate Survival Plot
#'
#' This is the primary workhorse function. It fits a Cox model to get HRs
#' and generates a publication-ready Kaplan-Meier-like plot.
#'
#' @param data Data frame.
#' @param model_predictors String. RHS of the formula (e.g. "Genotype + sex").
#' @param time_col String. Time column name.
#' @param status_col String. Status column name.
#' @param group_col String. Grouping variable for the plot lines.
#' @param hr_term String. The specific term to extract HR for.
#' @param hr_label String. Label for the HR annotation.
#' @param plot_title String. Plot title.
#' @param ylab String. Y-axis label.
#' @param legend_title String. Legend title.
#' @param palette_name String or Vector. Color palette to use.
#' @param left_trunc Logical. If TRUE, uses `BL_AGE` for left truncation.
#' @param skip_cox Logical. If TRUE, skips HR calculation.
#'
#' @return A `ggsurvplot` object.
run_and_plot_analysis <- function(data, model_predictors, time_col, status_col,
                                  group_col, hr_term, hr_label, plot_title,
                                  ylab, legend_title, palette_name = "brewer_dark2",
                                  left_trunc = FALSE, skip_cox = FALSE) {
  
  # 1. Run Cox Model for Annotation
  hr_text_string <- ""
  if (!skip_cox) {
    if (left_trunc) {
      fmla_str <- paste0("Surv(BL_AGE, ", time_col, ", ", status_col, "==1) ~ ", model_predictors)
    } else {
      fmla_str <- paste0("Surv(", time_col, ", ", status_col, "==1) ~ ", model_predictors)
    }
    
    # Try fitting model, handle errors gracefully
    cox_model <- tryCatch({
      coxph(as.formula(fmla_str), data = data)
    }, error = function(e) NULL)
    
    if (!is.null(cox_model)) {
      hr_text_string <- format_hr_text(cox_model, hr_term, hr_label)
    }
  }
  
  # 2. Fit Survival Curve for Plotting
  if (left_trunc) {
    surv_obj <- Surv(data$BL_AGE, data[[time_col]], data[[status_col]])
  } else {
    surv_obj <- Surv(data[[time_col]], data[[status_col]])
  }
  fit_formula <- as.formula(paste("surv_obj ~", group_col))
  fit <- surv_fit(fit_formula, data = data)
  
  # 3. Resolve Palette
  if(length(palette_name) > 1) {
    cols <- palette_name 
  } else {
    cols <- c("#1B9E77", "#D95F02", "#7570B3") # Default safe colors
  }

  # 4. Generate Plot
  p <- ggsurvplot(
    fit,
    data = data,
    title = plot_title,
    xlab = "Age in Years",
    ylab = ylab,
    palette = cols,
    risk.table = TRUE,
    censor = FALSE,
    conf.int = FALSE,
    legend.title = legend_title,
    ggtheme = theme_classic(base_size = 14)
  )
  
  # 5. Add HR Annotation
  if (hr_text_string != "") {
    p$plot <- p$plot + 
      annotate("text", x = min(data[[time_col]], na.rm=TRUE), 
               y = 0.15, label = hr_text_string, hjust = 0, size = 3.5)
  }
  
  return(p)
}

#' Style ggsurvplot
#' 
#' Applies standard styling: adds border, removes risk table titles, cleans margins.
#' @param ggsurv_object A ggsurvplot object.
#' @return A styled ggsurvplot object.
style_ggsurvplot <- function(ggsurv_object) {
  
  ggsurv_object$plot <- ggsurv_object$plot + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
    )
  
  ggsurv_object$table <- ggsurv_object$table + 
    theme(
      legend.position = "none",
      plot.title = element_blank(), # Remove "Number at risk" title
      plot.margin = margin(t = 0, r = 0, b = 10, l = 0)
    )
  
  return(ggsurv_object)
}

#' Combine Survival Plots
#'
#' Arranges multiple styled ggsurvplots into a grid with a shared legend.
#'
#' @param styled_list List of styled ggsurvplot objects.
#' @param ylabels Vector of Y-axis labels corresponding to the plots.
#' @param ncol_grid Integer. Number of columns in grid.
#' @param labels Vector. Labels for panels (e.g. "AUTO").
#' @param file Optional string. Filename to save PDF.
#' @param width Numeric. Width in inches.
#' @param height Numeric. Height in inches.
#' @return A cowplot object.
combine_survplots <- function(styled_list, ylabels, ncol_grid = 2, labels = "AUTO",
                              file = NULL, width = 8, height = 10) {
  
  # Extract Legend from first plot (assuming shared legend)
  first_plot_w_lgd <- styled_list[[1]]$plot + theme(legend.position = "bottom")
  shared_legend <- cowplot::get_legend(first_plot_w_lgd)
  
  # Combine Plot + Table for each item in list
  combos <- Map(function(sp, y_lab) {
    cowplot::plot_grid(
      sp$plot + labs(y = y_lab) + theme(legend.position = "none"),
      sp$table,
      ncol = 1,
      align = "v",
      rel_heights = c(0.75, 0.25)
    )
  }, styled_list, ylabels)
  
  # Main Grid
  grid_main <- cowplot::plot_grid(plotlist = combos, ncol = ncol_grid, labels = labels)
  
  # Combine Grid with Legend
  final_plot <- cowplot::plot_grid(
    grid_main,
    shared_legend,
    ncol = 1,
    rel_heights = c(0.9, 0.1)
  )
  
  # Save Output
  if (!is.null(file)) {
    ggsave(file, final_plot, width = width, height = height)
    message("Saved figure to: ", file)
  }
  
  return(final_plot)
}
