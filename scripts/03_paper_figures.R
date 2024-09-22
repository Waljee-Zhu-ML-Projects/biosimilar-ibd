rm(list = ls())
DATA_DIR <- here::here("data")
RESULTS_DIR <- here::here("results")
SURV_RESULTS_DIR <- here::here(RESULTS_DIR, "survival_models")
CLASS_RESULTS_DIR <- here::here(RESULTS_DIR, "classification_models")
FIGURES_DIR <- here::here(RESULTS_DIR, "paper_figures")

###################### Helper Functions ###################### 
rename_survival_methods <- function(x) {
  dplyr::case_when(
    x == "ltrc" ~ "LTRC RF",
    x == "cox" ~ "Cox",
    x == "cox_lasso" ~ "Cox + Lasso",
    x == "cox_ridge" ~ "Cox + Ridge"
  ) |>
    factor(levels = c("LTRC RF", "Cox + Lasso", "Cox + Ridge", "Cox"))
}

rename_classification_methods <- function(x) {
  dplyr::case_when(
    x == "autogluon_medium" ~ "Autogluon",
    x == "elnet" ~ "Elastic Net",
    x == "gb" ~ "GBDT",
    x == "lasso" ~ "Lasso",
    x == "logistic" ~ "Logistic",
    x == "randf" ~ "RF",
    x == "rfplus" ~ "RF+",
    x == "ridge" ~ "Ridge",
    TRUE ~ x
  ) |>
    factor(
      levels = c(
        "RF+", "RF", "GBDT", 
        "Lasso", "Ridge", "Elastic Net", "Logistic", "Autogluon"
      )
    )
}

rename_varset_modes <- function(x) {
  dplyr::case_when(
    x == "no_labs" ~ "No Lab Info",
    x == "recent_labs" ~ "Recent Labs Only",
    x == "recent_mean_max_labs" ~ "Mean, Max, Recent Labs Only",
    x == "recent_mean_max_recentolsslope_labs" ~ "Mean, Max, Slope, Recent Labs Only",
    x == "recent_accrual_labs" ~ "All Lab Info",
    TRUE ~ x
  ) |>
    factor(
      levels = c(
        "Mean, Max, Recent Labs Only",
        "Mean, Max, Slope, Recent Labs Only",
        "Recent Labs Only",
        "All Lab Info",
        "No Lab Info"
      )
    )
}

rename_variables <- function(x) {
  x <- stringr::str_remove(x, "\\.") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_remove("lab value ") |>
    stringr::str_remove("accrual ")
  x <- dplyr::case_when(
    x == "n flares preifx" ~ "# Flares Pre-IFX",
    x == "n prior flare" ~ "# Prior Flares",
    x == "count anyenc" ~ "# Any Encounters",
    x == "count IBDenc" ~ "# IBD Encounters",
    stringr::str_detect(x, "hgb") ~ 
      stringr::str_replace(x, "hgb", "Hemoglobin") |>
      stringr::str_to_title(),
    stringr::str_detect(x, "wbc") ~ 
      stringr::str_replace(x, "wbc", "White Blood Cells") |>
      stringr::str_to_title(),
    x == "treatment" ~ "Biosimilar Switch",
    x == "treatmentbiosimilar" ~ "Biosimilar Switch",
    x == "TimeYrs IBD IFX" ~ "IBD Duration",
    x == "priority" ~ "Enrollment Priority",
    x == "concomitant med" ~ "Concomitant Medications",
    x == "Age IFX" ~ "Age at IFX Index",
    stringr::str_detect(x, "co2") ~ 
      stringr::str_to_title(x) |>
      stringr::str_replace("Co2", "CO2"),
    stringr::str_detect(x, "rurality") ~ "Rurality",
    stringr::str_detect(x, "GenderF") ~ "Gender",
    stringr::str_detect(x, "Ethnicity$") ~ "Ethnicity",
    stringr::str_detect(x, "priority") ~ "Enrollment Priority",
    stringr::str_detect(x, "ProviderType gp") ~ "Provider Type",
    stringr::str_detect(x, "IBD type c") ~ "IBD Type",
    stringr::str_detect(x, "race gp$") ~ 
      stringr::str_replace(x, "race gp", "Race"),
    stringr::str_detect(x, "MS cat$") ~ 
      stringr::str_replace(x, "MS cat", "Marital Status"),
    stringr::str_detect(x, "COMORBIDITY cat$") ~ 
      stringr::str_replace(x, "COMORBIDITY cat", "Comorbidity"),
    stringr::str_detect(x, "Ethnicity") ~ 
      stringr::str_replace(x, "Ethnicity", "Ethnicity: ") |>
      stringr::str_to_title(),
    stringr::str_detect(x, "race gp") ~ 
      stringr::str_replace(x, "race gp", "Race: ") |>
      stringr::str_to_title(),
    stringr::str_detect(x, "MS cat") ~ 
      stringr::str_replace(x, "MS cat", "Marital Status: ") |>
      stringr::str_to_title(),
    stringr::str_detect(x, "COMORBIDITY cat") ~ 
      stringr::str_replace(x, "COMORBIDITY cat", "Comorbidity: ") |>
      stringr::str_to_title(),
    TRUE ~ stringr::str_to_title(x)
  )
  return(x)
}

load_classification_results <- function(valid_fname = NULL, test_fname = NULL) {
  valid_df <- NULL
  test_df <- NULL
  if (!is.null(valid_fname)) {
    valid_df <- data.table::fread(
      file.path(CLASS_RESULTS_DIR, valid_fname)
    ) |>
      dplyr::mutate(
        split_mode = "validation"
      )
  }
  if (!is.null(test_fname)) {
    test_df <- data.table::fread(
      file.path(CLASS_RESULTS_DIR, test_fname)
    ) |>
      dplyr::mutate(
        split_mode = "test"
      )
  }
  results_df <- dplyr::bind_rows(valid_df, test_df) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      accrual_name = sprintf("Accrual = %s year", accrual) |>
        forcats::fct_inorder(),
      forecast_name = sprintf("Forecast = %s months", forecast) |>
        forcats::fct_inorder(),
      method = rename_classification_methods(method),
      varset_mode = rename_varset_modes(varset_mode)
    )
  return(results_df)
}

###################### Survival Analysis Figures ###################### 

#### C-index results ####
cindex_df <- dplyr::bind_rows(
  readRDS(
    file.path(SURV_RESULTS_DIR, "validation_cindexs.rds")
  ),
  readRDS(
    file.path(SURV_RESULTS_DIR, "test_cindexs.rds")
  )
)

plt_df <- cindex_df |>
  dplyr::filter(
    varset_mode %in% c("no_labs", "labs_rf_imp"),
    censor_switcher == FALSE
  ) |>
  tidyr::pivot_longer(
    cols = c(ltrc, cox, cox_lasso, cox_ridge),
    names_to = "method",
    values_to = "cindex"
  ) |>
  dplyr::mutate(
    method = rename_survival_methods(method),
    split_mode = stringr::str_to_title(split_mode) |>
      factor(levels = c("Validation", "Test"))
  ) |>
  dplyr::mutate(
    time = time / 365.25,
    varset_mode = dplyr::case_when(
      varset_mode == "no_labs" ~ "No Labs",
      varset_mode == "labs_rf_imp" ~ "Include Labs"
    )
  ) |>
  dplyr::group_by(
    varset_mode, time, forecast, forecast_name, censor_switcher, split_mode, 
    method
  ) |>
  dplyr::summarise(
    cindex_mean = mean(cindex, na.rm = TRUE),
    cindex_sd = sd(cindex, na.rm = TRUE),
    cindex_se = sd(cindex, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

plt <- plt_df |>
  dplyr::filter(
    time > 1/12, time <= 4,
    split_mode == "Validation"
  ) |>
  ggplot2::ggplot() +
  ggplot2::aes(
    x = time, 
    y = cindex_mean,
    linetype = varset_mode
  ) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(method ~ forecast_name) +
  ggplot2::labs(
    x = "Accrual Period (yrs)",
    y = "C-index",
    linetype = ""
  ) +
  vthemes::theme_vmodern(size_preset = "medium")
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, "survival_cindex_by_variables_used.pdf"
  ),
  width = 10,
  height = 6
)

plt <- plt_df |>
  dplyr::filter(
    time > 1/12, time <= 4,
    split_mode == "Validation",
    varset_mode == "Include Labs"
  ) |>
  ggplot2::ggplot() +
  ggplot2::geom_line(
    ggplot2::aes(
      x = time,
      y = cindex_mean,
      color = method
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      x = time,
      ymin = cindex_mean - cindex_se,
      ymax = cindex_mean + cindex_se,
      fill = method
    ),
    alpha = 0.1
  ) +
  ggplot2::facet_grid(~ forecast_name) +
  vthemes::scale_color_vmodern(discrete = TRUE) +
  vthemes::scale_fill_vmodern(discrete = TRUE) +
  ggplot2::labs(
    x = "Accrual Period (yrs)",
    y = "C-index",
    color = "Method",
    fill = "Method",
    linetype = "Data Split"
  ) +
  vthemes::theme_vmodern(size_preset = "medium")
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, "survival_cindex_by_method.pdf"
  ),
  width = 10,
  height = 3.5
)

plt <- plt_df |>
  dplyr::filter(
    time > 1/12, time <= 4,
    varset_mode == "Include Labs"
  ) |>
  ggplot2::ggplot() +
  ggplot2::geom_line(
    ggplot2::aes(
      x = time,
      y = cindex_mean,
      color = method
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      x = time,
      ymin = cindex_mean - cindex_se,
      ymax = cindex_mean + cindex_se,
      fill = method
    ),
    alpha = 0.1
  ) +
  ggplot2::facet_grid(split_mode ~ forecast_name) +
  vthemes::scale_color_vmodern(discrete = TRUE) +
  vthemes::scale_fill_vmodern(discrete = TRUE) +
  ggplot2::labs(
    x = "Accrual Period (yrs)",
    y = "C-index",
    color = "Method",
    fill = "Method"
  ) +
  vthemes::theme_vmodern(size_preset = "medium")
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, "survival_cindex_by_method_and_split.pdf"
  ),
  width = 10,
  height = 5
)

#### Kaplan-Meier plots ####
test_km_df <- readRDS(
  file.path(SURV_RESULTS_DIR, "test_km.rds")
)

plt_df <- test_km_df |>
  dplyr::filter(
    censor_switcher == FALSE,
    varset_mode == "labs_rf_imp",
    !is.na(group),
    accrual == "0.5"
  ) |>
  dplyr::mutate(
    method = rename_survival_methods(method),
    strata = stringr::str_to_title(strata)
  ) |>
  # interpolate survival curves to have the same time points
  dplyr::group_by(
    strata, group, rep, varset_mode, censor_switcher,
    forecast, forecast_name, accrual, accrual_name, split_mode, method
  ) |>
  dplyr::summarise(
    dplyr::across(
      c(surv, std.err, upper, lower),
      ~ approx(
        x = time, y = .x, xout = seq(0, max(time), by = 0.01)
      ) |>
        tibble::as_tibble() |>
        setNames(c("time", dplyr::cur_column())) |>
        dplyr::select(-time) |>
        list()
    ),
    time = tibble::tibble(time = seq(0, max(time), by = 0.01)) |>
      list(),
    .groups = "drop"
  ) |>
  tidyr::unnest(
    c(time, surv, std.err, upper, lower)
  ) |>
  # aggregate survival curves
  dplyr::group_by(
    time, strata, group, varset_mode, censor_switcher,
    forecast, forecast_name, accrual, accrual_name, split_mode, method
  ) |>
  dplyr::summarise(
    dplyr::across(
      c(surv, std.err, upper, lower),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

plt <- plt_df |>
  dplyr::filter(time <= 4) |>
  ggplot2::ggplot() +
  ggplot2::geom_line(
    ggplot2::aes(
      x = time, y = surv, color = strata
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      x = time, ymin = lower, ymax = upper, fill = strata
    ),
    alpha = 0.2
  ) +
  ggplot2::facet_grid(
    forecast_name ~ method
  ) +
  ggplot2::labs(
    x = "Time (yrs)", y = "Flare-free Survival",, 
    color = "Risk Strata", fill = "Risk Strata"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(reverse = TRUE),
    fill = ggplot2::guide_legend(reverse = TRUE)
  ) +
  vthemes::theme_vmodern(size_preset = "medium")
ggplot2::ggsave(
  plt,
  filename = file.path(FIGURES_DIR, "survival_km_by_risk.pdf"),
  width = 10,
  height = 6
)

#### feature importances ####
vimps_df <- readRDS(
  file.path(SURV_RESULTS_DIR, "test_importances.rds")
)

plt_df_all <- vimps_df |>
  dplyr::mutate(
    method = rename_survival_methods(method),
    variable = rename_variables(variable)
  ) |>
  dplyr::filter(
    censor_switcher == FALSE,
    varset_mode == "labs_rf_imp"
  ) |>
  dplyr::group_by(
    forecast, forecast_name, censor_switcher, varset_mode, method, variable
  ) |>
  dplyr::summarise(
    dplyr::across(
      c(importance, estimate, statistic, p.value, HR),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        se = ~ sd(.x, na.rm = TRUE) / sqrt(dplyr::n())
      )
    ),
    .groups = "drop"
  )

plt_ls <- list()
for (method_name in levels(plt_df_all$method)) {
  plt_df <- plt_df_all |>
    dplyr::filter(
      method == !!method_name
    )
  var_order <- plt_df |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      mean_importance = mean(importance_mean)
    ) |>
    dplyr::arrange(-mean_importance) |>
    dplyr::pull(variable)
  plt <- plt_df |>
    dplyr::mutate(
      variable = factor(variable, levels = rev(var_order))
    ) |>
    ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(
        x = importance_mean, y = variable
      ),
      size = 0.75
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = importance_mean - importance_se, 
        xmax = importance_mean + importance_se,
        y = variable
      ),
      width = 0.1
    ) +
    ggplot2::facet_grid(method ~ forecast_name) +
    ggplot2::labs(
      x = "Importance", y = "Variable"
    ) +
    vthemes::theme_vmodern(size_preset = "medium")
  if (method_name != unique(plt_df_all$method)[1]) {
    plt <- plt +
      ggplot2::theme(
        strip.background.x = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank()
      )
  }
  plt_ls[[method_name]] <- plt
}

plt <- patchwork::wrap_plots(
  plt_ls, ncol = 1, heights = c(1, 1.15, 1.15, 1.15)
) +
  patchwork::plot_layout(axes = "collect")
ggplot2::ggsave(
  plt,
  filename = file.path(FIGURES_DIR, "survival_importance.pdf"),
  width = 12,
  height = 17
)

###################### Classification Analysis Figures ###################### 

#### class proportions ####
class_prop_ls <- list()
for (forecast in c(3, 6, 12)) {
  for (accrual in c(1, 2)) {
    data_fname <- file.path(
      DATA_DIR,
      "classification_data",
      sprintf(
        "classification_accrual_%syr_forecast_%smo_rf_imputed_rep1.csv",
        accrual, forecast
      )
    )
    data_df <- data.table::fread(data_fname)
    class_prop_ls <- c(
      class_prop_ls, 
      list(
        tibble::tibble(
          forecast = !!forecast,
          accrual = !!accrual,
          n_events = sum(data_df$event),
          n_total = nrow(data_df)
        )
      )
    )
  }
}

class_prop_df <- dplyr::bind_rows(class_prop_ls) |>
  dplyr::mutate(
    class_prop = n_events / n_total
  )
saveRDS(
  class_prop_df, 
  file = file.path(CLASS_RESULTS_DIR, "class_proportions.rds")
)

#### prediction performance ####
metrics <- c("auroc", "auprc")
rm_metrics <- c("f1", "accuracy")
rm_varset_modes <- c("Mean, Max, Slope, Recent Labs Only")
errs_df <- load_classification_results(
  valid_fname = "validation_errors.csv",
  test_fname = "test_errors.csv"
) |>
  dplyr::select(
    -tidyselect::all_of(rm_metrics)
  )

errs_df_long <- errs_df |>
  tidyr::pivot_longer(
    cols = tidyselect::all_of(metrics),
    names_to = "metric",
    values_to = "estimate"
  )

for (split_mode in c("validation", "test")) {
  for (metric in metrics) {
    metric_name <- dplyr::case_when(
      metric == "auroc" ~ "AUROC",
      metric == "auprc" ~ "AUPRC",
      TRUE ~ stringr::str_to_title(metric)
    )
    
    # prediction error plot
    plt_df <- errs_df_long |>
      dplyr::filter(
        split_mode == !!split_mode,
        is.na(censor_switcher),
        impute_mode == "rf",
        !(varset_mode %in% !!rm_varset_modes),
        metric == !!metric
      ) |>
      dplyr::group_by(
        method, censor_switcher, impute_mode,
        forecast, forecast_name, accrual, accrual_name, varset_mode, metric
      ) |>
      dplyr::summarise(
        mean_estimate = mean(estimate),
        sd_estimate = sd(estimate),
        se_estimate = sd(estimate) / sqrt(dplyr::n()),
        n_estimate = dplyr::n(),
        .groups = "drop"
      )
    plt <- plt_df |>
      dplyr::mutate(
        forecast = forcats::fct_inseq(as.factor(forecast))
      ) |>
      ggplot2::ggplot() +
      ggplot2::geom_line(
        ggplot2::aes(
          x = forecast, 
          y = mean_estimate, 
          color = varset_mode,
          group = varset_mode
        )
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          x = forecast,
          ymin = mean_estimate - se_estimate,
          ymax = mean_estimate + se_estimate,
          fill = varset_mode,
          group = varset_mode
        ),
        alpha = 0.1
      ) +
      ggplot2::facet_grid(accrual_name ~ method, scales = "free") +
      ggplot2::labs(
        x = "Forecast Period (months)",
        y = metric_name, 
        color = "Variables Used",
        fill = "Variables Used"
      ) +
      vthemes::scale_color_vmodern(discrete = TRUE) +
      vthemes::scale_fill_vmodern(discrete = TRUE) +
      vthemes::theme_vmodern(size_preset = "medium")
    ggplot2::ggsave(
      plt,
      filename = file.path(
        FIGURES_DIR, 
        sprintf("classification_%s_%s.pdf", split_mode, metric)
      ),
      width = 16,
      height = 5
    )
    
    # compare RF+ to logistic
    if (split_mode == "validation") {
      plt_df <- errs_df_long |>
        dplyr::filter(
          split_mode == !!split_mode,
          is.na(censor_switcher),
          impute_mode == "rf",
          varset_mode == "Mean, Max, Recent Labs Only",
          metric == !!metric
        ) |>
        tidyr::pivot_wider(
          names_from = method,
          values_from = estimate
        ) |>
        dplyr::mutate(
          dplyr::across(
            tidyselect::all_of(unique(errs_df_long$method)),
            ~ `RF+` - .x
          )
        ) |>
        tidyr::pivot_longer(
          cols = tidyselect::all_of(unique(errs_df_long$method)),
          names_to = "method",
          values_to = "diff"
        ) |>
        dplyr::filter(
          method != "RF+"
        ) |>
        dplyr::mutate(
          method = rename_classification_methods(method)
        )
      plt_pvals_df <- plt_df |>
        dplyr::group_by(
          censor_switcher, impute_mode, varset_mode, split_mode, metric,
          forecast, forecast_name, accrual, accrual_name, method
        ) |>
        dplyr::summarise(
          perm_pval = purrr::map_lgl(
            1:1000,
            ~ mean(
              sample(c(1, -1), size = length(diff), replace = TRUE) * diff
            ) >
              mean(diff)
          ) |> 
            mean(),
          t_pval = t.test(
            diff, mu = 0, alternative = "greater"
          )$p.value,
          wilcox_pval = wilcox.test(
            diff, mu = 0, alternative = "greater"
          )$p.value,
          .groups = "drop"
        ) |>
        dplyr::mutate(
          pval_label = dplyr::case_when(
            perm_pval < 0.001 ~ "***",
            perm_pval < 0.01 ~ "**",
            perm_pval < 0.05 ~ "*",
            TRUE ~ ""
          )
        )
      plt <-  plt_df |>
        ggplot2::ggplot() +
        ggplot2::geom_boxplot(
          ggplot2::aes(
            x = method, y = diff 
          ),
          outlier.size = 0.75
        ) +
        ggplot2::geom_text(
          ggplot2::aes(
            x = method, y = Inf, label = pval_label
          ),
          angle = 90, hjust = 1, vjust = 0.75, size = 7,
          data = plt_pvals_df
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::facet_grid(accrual_name ~ forecast_name, scales = "free_y") +
        ggplot2::labs(
          x = "Method", 
          y = sprintf("RF+ %s - Model %s", metric_name, metric_name)
        ) +
        vthemes::theme_vmodern(
          size_preset = "medium"
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            size = 12, angle = 90, hjust = 1, vjust = 0.5
          )
        )
      ggplot2::ggsave(
        plt,
        filename = file.path(
          FIGURES_DIR, 
          sprintf("classification_%s_%s_comparison.pdf", split_mode, metric)
        ),
        width = 10,
        height = 6
      )
    }
  }
}

# look specifically at mean, max, recent labs validation AUROC
plt_df <- errs_df_long |>
  dplyr::filter(
    split_mode == "validation",
    is.na(censor_switcher),
    impute_mode == "rf",
    varset_mode == "Mean, Max, Recent Labs Only",
    metric == "auroc"
  ) |>
  dplyr::group_by(
    method, censor_switcher, impute_mode,
    forecast, forecast_name, accrual, accrual_name, varset_mode, metric
  ) |>
  dplyr::summarise(
    mean_estimate = mean(estimate),
    sd_estimate = sd(estimate),
    se_estimate = sd(estimate) / sqrt(dplyr::n()),
    n_estimate = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    method_type = dplyr::case_when(
      method %in% c("RF", "GBDT") ~ "tree",
      method %in% c("Lasso", "Ridge", "Elastic Net", "Logistic") ~ "linear",
      TRUE ~ method
    )
  )
plt <- plt_df |>
  ggplot2::ggplot() +
  ggplot2::geom_point(
    ggplot2::aes(x = method, y = mean_estimate, color = method_type)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      x = method,
      ymin = mean_estimate - se_estimate,
      ymax = mean_estimate + se_estimate,
      color = method_type
    ),
    width = 0
  ) +
  ggplot2::facet_grid(accrual_name ~ forecast_name, scales = "free") +
  ggplot2::labs(
    x = "Method",
    y = "AUROC"
  ) +
  ggplot2::scale_color_manual(
    values = c("grey70", "#F8766D", "black", "#14CED3")
  ) +
  vthemes::theme_vmodern(size_preset = "medium") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 12, angle = 90, hjust = 1, vjust = 0.5
    )
  ) +
  ggplot2::guides(color = "none")
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, 
    "classification_validation_recent_mean_max_labs_auroc.pdf"
  ),
  width = 7,
  height = 5
)

#### sensitivity/specificity at best cutoffs ####
cutoffs_df <- load_classification_results(
  valid_fname = "validation_best_cutoff_results.csv",
  test_fname = "test_best_cutoff_results.csv"
)
errs_df <- load_classification_results(
  valid_fname = "validation_errors.csv",
  test_fname = "test_errors.csv"
) |>
  dplyr::left_join(
    cutoffs_df,
    by = c(
      "seed", "method", "censor_switcher", "impute_mode", "split_mode",
      "forecast", "forecast_name", "accrual", "accrual_name", "varset_mode"
    )
  )

test_errs_df <- errs_df |>
  dplyr::filter(
    split_mode == "test",
    method == "RF+",
    is.na(censor_switcher),
    impute_mode == "rf",
    varset_mode == "Mean, Max, Recent Labs Only"
  ) |>
  dplyr::group_by(
    method, censor_switcher, impute_mode, split_mode,
    forecast, forecast_name, accrual, accrual_name, varset_mode
  ) |>
  dplyr::summarise(
    dplyr::across(
      c(auroc, auprc, sensitivity, specificity, thr),
      list(
        mean = ~ sprintf("%.2f (%.2f)", mean(.x), sd(.x)),
        ci = ~ sprintf("(%.2f, %.2f)", quantile(.x, 0.025), quantile(.x, 0.975))
      )
    ),
    .groups = "drop"
  ) |>
  dplyr::select(
    `Forecast (months)` = forecast,
    `Accrual (years)` = accrual,
    auroc_mean:thr_ci
  )
class_prop_df <- readRDS(file.path(CLASS_RESULTS_DIR, "class_proportions.rds"))

kab_df <- dplyr::left_join(
  class_prop_df |>
    dplyr::mutate(
      forecast = as.integer(forecast),
      accrual = as.integer(accrual),
      class_prop = sprintf("%.2f", class_prop)
    ) |>
    dplyr::select(
      `Forecast (months)` = forecast, `Accrual (years)` = accrual, class_prop
    ),
  test_errs_df,
  by = c("Forecast (months)", "Accrual (years)")
)
kab <- vthemes::pretty_kable(
  kab_df,
  col.names = c(
    "Forecast<br>(months)", "Accrual<br>(years)", "Prop.<br>Flares",
    rep(c("Mean (SD)", "95% CI"), 5)
  )
) |>
  kableExtra::add_header_above(
    c(
      " " = 3, 
      "AUROC" = 2, 
      "AUPRC" = 2, 
      "Sensitivity" = 2, 
      "Specificity" = 2, 
      "Threshold" = 2
    )
  ) |>
  kableExtra::kable_paper() #> exported as image (1250 x 453)
saveRDS(
  kab,
  file = file.path(FIGURES_DIR, "test_errors_table.rds")
)

#### calibration ####
calibration_df <- load_classification_results(
  test_fname = "test_calibration_results.csv"
) |>
  dplyr::filter(
    method == "RF+",
    is.na(censor_switcher),
    impute_mode == "rf",
    varset_mode == "Mean, Max, Recent Labs Only",
    split_mode == "test"
  )

plt_df <- calibration_df |>
  dplyr::group_by(
    method, quantile, censor_switcher, impute_mode,
    forecast, forecast_name, accrual, accrual_name, varset_mode, split_mode
  ) |>
  dplyr::summarise(
    num_flares = mean(num_flares),
    predicted_flares = mean(predicted_flares),
    .groups = "drop"
  )

plt <- plt_df |>
  ggplot2::ggplot() +
  ggplot2::aes(
    x = num_flares, 
    y = predicted_flares, 
    color = as.factor(forecast), 
    shape = as.factor(accrual)
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(
    x = "Number of Observed Flares", 
    y = "Number of Predicted Flares",
    color = "Forecast (mo)",
    shape = "Accrual (yr)"
  ) +
  vthemes::theme_vmodern(size_preset = "medium")

ggplot2::ggsave(
  plt,
  filename = file.path(FIGURES_DIR, "calibration_plot.pdf"),
  height = 4,
  width = 6
)

#### decision curve analysis ####
dca_df <- load_classification_results(test_fname = "test_dca_results.csv") |>
  dplyr::filter(
    is.na(censor_switcher), 
    impute_mode == "rf",
    varset_mode == "Mean, Max, Recent Labs Only"
  )

best_cutoffs <- load_classification_results(
  test_fname = "test_best_cutoff_results.csv"
) |>
  dplyr::filter(
    is.na(censor_switcher), 
    impute_mode == "rf",
    varset_mode == "Mean, Max, Recent Labs Only",
    method == "RF+"
  ) |>
  dplyr::group_by(
    method, censor_switcher, impute_mode,
    forecast, forecast_name, accrual, accrual_name, varset_mode
  ) |>
  dplyr::summarise(
    sensitivity = mean(sensitivity),
    specificity = mean(specificity),
    thr = mean(thr),
    .groups = "drop"
  ) |>
  dplyr::rename(model = method)

plt_df <- dca_df |>
  dplyr::group_by(
    model, threshold, method, censor_switcher, impute_mode, 
    forecast, forecast_name, accrual, accrual_name, varset_mode
  ) |>
  dplyr::summarise(
    net_benefit = mean(net_benefit),
    net_benefit_sd = sd(net_benefit),
    net_benefit_se = sd(net_benefit) / sqrt(dplyr::n()),
    net_benefit_n = dplyr::n(),
    .groups = "drop"
  )

plt <- plt_df |>
  dplyr::filter(
    !(method == "Logistic" & model != "pred"),
    threshold <= 0.35,
    net_benefit >= -0.1
  ) |>
  dplyr::mutate(
    model = dplyr::case_when(
      model == "all" ~ "Treat All",
      model == "none" ~ "Treat None",
      method == "Logistic" ~ "Logistic Model",
      method == "RF+" ~ "RF+ Model"
    ) |>
      factor(
        levels = c("RF+ Model", "Logistic Model", "Treat All", "Treat None")
      )
  ) |>
  ggplot2::ggplot() +
  ggplot2::geom_line(
    ggplot2::aes(
      x = threshold, y = net_benefit, color = model
    )
  ) +
  ggplot2::geom_vline(
    ggplot2::aes(xintercept = thr),
    data = best_cutoffs,
    linetype = "dotted"
  ) +
  ggplot2::facet_grid(accrual_name ~ forecast_name) + 
  ggplot2::scale_color_manual(
    values = c(
      "black", "#F8766D", "#619CFF", "#00BA38"
    )
  ) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::labs(
    x = "Threshold Probability",
    y = "Net Benefit",
    color = ""
  ) +
  vthemes::theme_vmodern(size_preset = "medium")
ggplot2::ggsave(
  plt, 
  filename = file.path(FIGURES_DIR, "dca_all.pdf"),
  width = 8,
  height = 5
)

#### feature importances ####
vimps_df <- load_classification_results(
  test_fname = "test_importances.csv"
) |>
  dplyr::filter(
    is.na(censor_switcher),
    impute_mode == "rf",
    varset_mode == "Mean, Max, Recent Labs Only"
  ) |>
  dplyr::mutate(
    var = rename_variables(var)
  )
rfplus_vimps_df <- vimps_df |>
  dplyr::filter(
    method == "RF+"
  ) |>
  tidyr::pivot_longer(
    cols = c(importance, importance0, importance1),
    names_to = "metric",
    values_to = "importance"
  )

# RF+ Importance
plt_ls <- list()
for (forecast in c(3, 6, 12)) {
  for (accrual in c(1, 2)) {
    plt_name <- sprintf(
      "Accrual = %syr, Forecast = %smo", accrual, forecast
    )
    plt_df <- rfplus_vimps_df |>
      dplyr::filter(
        forecast == !!forecast,
        accrual == !!accrual
      ) |>
      dplyr::group_by(
        forecast_name, accrual_name, var, metric
      ) |>
      dplyr::summarise(
        mean_importance = mean(importance),
        sd_importance = sd(importance),
        se_importance = sd(importance) / sqrt(dplyr::n()),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        importance_mode = dplyr::case_when(
          is.na(stringr::str_extract(metric, "[0-9]+")) ~ "Total",
          stringr::str_extract(metric, "[0-9]+") == "0" ~ "Nonlinear",
          stringr::str_extract(metric, "[0-9]+") == "1" ~ "Linear"
        ) |>
          factor(levels = c("Total", "Linear", "Nonlinear"))
      )
    var_order <- plt_df |>
      dplyr::filter(importance_mode == "Total") |>
      dplyr::arrange(-mean_importance) |>
      dplyr::pull(var)
    plt_ls[[plt_name]] <- plt_df |>
      dplyr::mutate(
        var = factor(var, levels = rev(var_order[1:15]))
      ) |>
      dplyr::filter(
        !is.na(var)
      ) |>
      ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(
          x = var,
          y = mean_importance,
          color = importance_mode
        )
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          x = var,
          ymin = mean_importance - se_importance,
          ymax = mean_importance + se_importance,
          color = importance_mode
        ),
        width = 0
      ) +
      ggplot2::scale_color_manual(
        values = c("black", "#F8766D", "#14CED3")
      ) +
      ggplot2::facet_grid(accrual_name ~ forecast_name) +
      ggplot2::labs(
        x = "Variable", 
        y = "MDI+ Importance", 
        color = "Importance"
      ) +
      ggplot2::coord_flip() +
      vthemes::theme_vmodern(
        size_preset = "large"
      )
    if (accrual == 2) {
      plt_ls[[plt_name]] <- plt_ls[[plt_name]] +
        ggplot2::theme(
          strip.text.x = ggplot2::element_blank(),
          strip.background.x = ggplot2::element_blank()
        )
    }
    if (forecast != 12) {
      plt_ls[[plt_name]] <- plt_ls[[plt_name]] +
        ggplot2::theme(
          strip.text.y = ggplot2::element_blank(),
          strip.background.y = ggplot2::element_blank()
        )
    }
  }
}

plt <- patchwork::wrap_plots(
  plt_ls, nrow = 2, ncol = 3, byrow = FALSE, guides = "collect"
) +
  patchwork::plot_layout(axis_titles = "collect")
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, "rfplus_importances_linear_v_nonlinear.pdf"
  ),
  height = 8,
  width = 20
)

# Importance across methods
vimps_rank_df <- vimps_df |>
  dplyr::group_by(
    method, censor_switcher, impute_mode, varset_mode,
    forecast, forecast_name, accrual, accrual_name, split_mode, var
  ) |>
  dplyr::summarise(
    importance_mean = mean(importance),
    .groups = "drop"
  ) |>
  dplyr::group_by(
    method, censor_switcher, impute_mode, varset_mode,
    forecast, forecast_name, accrual, accrual_name, split_mode
  ) |>
  dplyr::mutate(
    importance_rank = rank(-importance_mean)
  )

plt_ls <- list()
for (forecast in c(3, 6, 12)) {
  for (accrual in c(1, 2)) {
    plt_name <- sprintf(
      "Accrual = %syr, Forecast = %smo", accrual, forecast
    )
    plt_df <- vimps_rank_df |>
      dplyr::filter(
        forecast == !!forecast,
        accrual == !!accrual
      )
    var_order <- plt_df |>
      dplyr::filter(method == "RF+") |>
      dplyr::arrange(importance_rank) |>
      dplyr::pull(var)
    plt_ls[[plt_name]] <- plt_df |>
      dplyr::mutate(
        var = factor(var, levels = rev(var_order[1:10]))
      ) |>
      dplyr::filter(
        !is.na(var)
      ) |>
      ggplot2::ggplot() +
      ggplot2::aes(
        x = method, y = var, fill = importance_rank, label = importance_rank
      ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(size = 3) +
      ggplot2::facet_grid(accrual_name ~ forecast_name) +
      ggplot2::labs(
        x = "Method", 
        y = "Variable", 
        fill = "Variable\nImportance\nRank"
      ) +
      viridis::scale_fill_viridis(
        option = "magma", direction = -1, limits = c(1, 50),
        guide = ggplot2::guide_colorbar(reverse = TRUE),
        begin = 0.2
      ) +
      ggplot2::coord_cartesian(expand = FALSE) +
      vthemes::theme_vmodern(
        size_preset = "large"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          size = 12, angle = 90, hjust = 1, vjust = 0.5
        ),
        panel.border = ggplot2::element_rect(
          color = "black", fill = NA, linewidth = 1
        )
      )
    if (accrual == 2) {
      plt_ls[[plt_name]] <- plt_ls[[plt_name]] +
        ggplot2::theme(
          strip.text.x = ggplot2::element_blank(),
          strip.background.x = ggplot2::element_blank()
        )
    }
    if (forecast != 12) {
      plt_ls[[plt_name]] <- plt_ls[[plt_name]] +
        ggplot2::theme(
          strip.text.y = ggplot2::element_blank(),
          strip.background.y = ggplot2::element_blank()
        )
    }
  }
}

plt <- patchwork::wrap_plots(
  plt_ls, nrow = 2, ncol = 3, byrow = FALSE, guides = "collect"
) +
  patchwork::plot_layout(axis_titles = "collect")
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, "importances_summary.pdf"
  ),
  height = 9,
  width = 18
)

#### partial dependence plots ####
varset_mode <- "recent_mean_max_labs"
censor_switcher <- ""
impute_mode <- "rf"
seeds <- 1:50
features <- c(
  "# Prior Flares" = "n_prior_flare",
  "# Any Encounters" = "count_anyenc",
  "# IBD Encounters" = "count_IBDenc",
  "Max White Blood Cells" = "max_lab_value_accrual_wbc",
  "Recent Hemoglobin" = "recent_lab_value_hgb",
  "Mean Hemoglobin" = "mean_lab_value_accrual_hgb",
  "Max Platelet" = "max_lab_value_accrual_platelet",
  "# Flares Pre-IFX" = "n_flares_preifx",
  "Mean CO2" = "mean_lab_value_accrual_co2"
)
plt_ls_all <- list()
for (forecast in c(3, 6, 12)) {
  for (accrual in c(1, 2)) {
    plt_name <- sprintf(
      "Accrual = %syr, Forecast = %smo", accrual, forecast
    )
    var_order <- vimps_rank_df |>
      dplyr::filter(
        method == "RF+",
        accrual == !!accrual,
        forecast == !!forecast
      ) |>
      dplyr::arrange(importance_rank) |>
      dplyr::pull(var)
    keep_vars <- var_order[1:6]
    plt_ls <- list()
    for (feature_name in keep_vars) {
      feature <- features[[feature_name]]
      Xs_ls <- list()
      plt_df_ls <- list()
      for (seed in seeds) {
        pdp_df <- data.table::fread(
          file.path(
            CLASS_RESULTS_DIR,
            sprintf(
              "classification_accrual_%syr_forecast_%smo%s_%s_imputed_rep%s",
              accrual, forecast, censor_switcher, impute_mode, seed
            ),
            sprintf("test_%s_pdp_%s.csv", varset_mode, feature)
          )
        ) |>
          tibble::as_tibble()
        
        X <- data.table::fread(
          file.path(
            DATA_DIR,
            "classification_data",
            sprintf(
              "classification_accrual_%syr_forecast_%smo%s_%s_imputed_rep%s.csv",
              accrual, forecast, censor_switcher, impute_mode, seed
            )
          )
        ) |>
          tibble::as_tibble()
        feature_data <- X[[feature]]
        x_mean <- mean(feature_data)
        x_sd <- sd(feature_data)
        
        plt_df_ls[[seed]] <- pdp_df |>
          dplyr::mutate(
            orig_grid_values = grid_values * x_sd + x_mean
          )
        Xs_ls[[seed]] <- X |>
          dplyr::select(tidyselect::all_of(feature))
      }
      
      plt_df <- plt_df_ls |>
        dplyr::bind_rows(.id = "seed")
      common_grid_values <- seq(
        quantile(plt_df$orig_grid_values, 0.01),
        quantile(plt_df$orig_grid_values, 0.99),
        length.out = 100
      )
      plt_df <- plt_df |>
        dplyr::group_by(seed) |>
        dplyr::summarise(
          imputed_pdp = approx(
            x = orig_grid_values, 
            y = average, 
            xout = common_grid_values
          ) |>
            tibble::as_tibble() |>
            list(),
          .groups = "drop"
        ) |>
        tidyr::unnest(cols = imputed_pdp) |>
        dplyr::group_by(x) |>
        dplyr::summarise(
          mean = mean(y, na.rm = TRUE),
          lo = quantile(y, 0.025, na.rm = TRUE),
          hi = quantile(y, 0.975, na.rm = TRUE)
        )
      
      plt_ls[[feature]] <- plt_df |>
        ggplot2::ggplot() +
        ggplot2::aes(
          x = x,
          y = mean
        ) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = x, ymin = lo, ymax = hi),
          alpha = 0.2
        ) +
        ggplot2::geom_rug(
          ggplot2::aes(
            x = .data[[feature]]
          ),
          inherit.aes = FALSE,
          data = dplyr::bind_rows(Xs_ls),
          alpha = 1 / length(seeds) * 0.2
        ) +
        ggplot2::labs(
          x = feature_name, 
          y = "Predicted Risk of Flare"
        ) +
        vthemes::theme_vmodern()
    }
    
    plt <- patchwork::wrap_plots(plt_ls, nrow = 1, ncol = length(plt_ls)) +
      patchwork::plot_layout(axis_titles = "collect_y") +
      patchwork::plot_annotation(
        title = sprintf("Accrual = %syr, Forecast = %smo", accrual, forecast),
        theme = ggplot2::theme(
          title = ggplot2::element_text(face = "bold")
        )
      )
    ggplot2::ggsave(
      plt,
      filename = file.path(
        FIGURES_DIR, sprintf("pdp_forecast%s_accrual%s.pdf", forecast, accrual)
      ),
      width = 12,
      height = 3
    )
    plt_ls_all[[plt_name]] <- patchwork::wrap_elements(plt)
  }
}

plt <- patchwork::wrap_plots(plt_ls_all, nrow = length(plt_ls_all), ncol = 1)
ggplot2::ggsave(
  plt,
  filename = file.path(
    FIGURES_DIR, "pdp_plots.pdf"
  ),
  width = 12,
  height = 18
)


