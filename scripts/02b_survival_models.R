rm(list = ls())
library(future)
library(survival)
library(prodlim)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
forecast_i <- as.integer(args[2])
censor_switcher_i <- as.integer(args[3])

# initialize directory variables and load in helper functions
USE_CACHED <- TRUE
DATA_DIR <- here::here("data", "survival_data")
RESULTS_DIR <- here::here("results", "survival_models")
FUNCTIONS_DIR <- here::here("functions")
for (fname in list.files(FUNCTIONS_DIR, pattern = "\\.R$", full.names = TRUE)) {
  source(fname, chdir = TRUE)
}

# helper variables
yr <- 365.25
month <- 30
forecasts <- c("3" = 3 * month, "6" = 0.5 * yr, "12" = 1 * yr)
accruals <- c("0.5" = 0.5, "1" = 1, "1.5" = 1.5, "2" = 2, "2.5" = 2.5) * yr
censor_switchers <- c(FALSE, "_censor_switcher" = TRUE)

censor_switcher <- censor_switchers[censor_switcher_i]
censor_switcher_name <- names(censor_switchers)[censor_switcher_i]
forecast <- forecasts[forecast_i]
forecast_name <- names(forecasts)[forecast_i]
data_fname <- file.path(
  DATA_DIR,
  sprintf("survival_forecast_%smo%s.RData", forecast_name, censor_switcher_name)
)
OUT_DIR <- file.path(
  RESULTS_DIR, sprintf("seed%s", seed), sprintf("forecast%s", forecast_name)
)
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# set up clusters
B <- 20
# n_cores <- 10
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
plan(multicore, workers = n_cores)

# load in data
load(data_fname)

# remove labs (too much missing data or irrelevant)
rm_labs <- c("mch", "mchc", "mpv", "rbc", "hematocrit", "neutrophils", "FeCal", "crp")
lab_df_orig <- load_lab_data()
lab_names <- setdiff(unique(lab_df_orig$LabName), rm_labs)

# helper variables for survival modeling
formula <- "Surv(event_censor_date, event) ~ ."
formula_timevar <- "Surv(start_date, end_date, event) ~ ."
eval_formula <- Surv(event_censor_date, event) ~ 1
keep_vars <- c(
  "switch_factor", "COMORBIDITY_cat", "Age_IFX", "Gender", "race_gp",
  "Ethnicity", "MS_cat", "priority", "concomitant_med", "TimeYrs_IBD_IFX",
  "IBD_type_c", "rurality", "count_anyenc", "count_IBDenc", "ProviderType_gp",
  "n_flares_preifx"
)

# set seed
set.seed(317 + seed)

# do data splitting
deids <- unique(survival_df$deId)
split_idx <- sample(
  c("train", "valid", "test"), 
  size = length(deids), 
  replace = TRUE, 
  prob = c(0.6, 0.15, 0.25)
)
train_deids <- deids[split_idx == "train"]
valid_deids <- deids[split_idx == "valid"]
test_deids <- deids[split_idx == "test"]
save(
  train_deids, valid_deids, test_deids,
  file = file.path(
    OUT_DIR, 
    sprintf(
      "split_indices_forecast_%smo%s.RData", 
      forecast_name, censor_switcher_name
    )
  )
)

# get common evaluation times for all data splits
eval_times_fname <- file.path(
  RESULTS_DIR, 
  sprintf("eval_times_forecast%s%s.rds", forecast_name, censor_switcher_name)
)
if (file.exists(eval_times_fname)) {
  # read eval_times from cache
  eval_times <- readRDS(eval_times_fname)
} else {
  # get eval_times from one run
  base_fit_df <- survival_df |>
    dplyr::select(tidyselect::all_of(keep_vars), event_censor_date, event)
  train_df <- base_fit_df[survival_df$deId %in% train_deids, ]
  test_df <- base_fit_df[survival_df$deId %in% valid_deids, ]
  
  # fit survival forest to get reasonable eval times
  base_model <- fit_rfsrc(
    formula = formula,
    train_df = train_df, 
    test_df = test_df
  )
  # retrieve unique death times for evaluation
  eval_times <- sort(c(base_model$fit$time.interest, accruals))
  saveRDS(
    eval_times,
    file = eval_times_fname
  )
}

# time-varying data frame without labs
tfit_df <- survival_df |>
  dplyr::select(
    tidyselect::all_of(keep_vars),
    start_date, end_date, event, deId
  )

# non-time-varying data frame without labs used solely for permutation importance
cov_df <- tfit_df |>
  dplyr::group_by(deId) |>
  dplyr::mutate(
    switch_factor = dplyr::case_when(
      any(switch_factor == "biosimilar") ~ "biosimilar",
      TRUE ~ "originator"
    ) |> 
      factor(levels = c("originator", "biosimilar")),
    event = max(event)
  ) |>
  dplyr::select(-start_date, -end_date) |>
  dplyr::ungroup() |>
  dplyr::distinct()

# time-varying + labs (median-imputation) data frame
tfit_lab_df_median_imp <- survival_lab_df_median_imp |>
  dplyr::select(
    tidyselect::all_of(keep_vars), 
    tidyselect::any_of(lab_names),
    start_date, end_date, event, deId
  )

# time-varying + labs (mean-imputation) data frame
tfit_lab_df_mean_imp <- survival_lab_df_mean_imp |>
  dplyr::select(
    tidyselect::all_of(keep_vars), 
    tidyselect::any_of(lab_names),
    start_date, end_date, event, deId
  )

# time-varying + labs (missforest-imputation) data frame
tfit_lab_df_rf_imp <- survival_lab_df_rf_imp |>
  dplyr::select(
    tidyselect::all_of(keep_vars), 
    tidyselect::any_of(lab_names),
    start_date, end_date, event, deId
  )

fit_df_list <- list(
  "no_labs" = tfit_df,
  "labs_median_imp" = tfit_lab_df_median_imp,
  "labs_mean_imp" = tfit_lab_df_mean_imp,
  "labs_rf_imp" = tfit_lab_df_rf_imp
)

for (split_mode in c("validation", "test")) {
  model_ls <- list()
  cindex_ls <- list()
  km_ls <- list()
  for (fit_mode in names(fit_df_list)) {
    keep_models <- c()
    fit_df <- fit_df_list[[fit_mode]]
    
    #### Data Split ####
    if (split_mode == "validation") {
      train_idx <- fit_df$deId %in% train_deids
      test_idx <- fit_df$deId %in% valid_deids
    } else if (split_mode == "test") {
      train_idx <- fit_df$deId %in% c(train_deids, valid_deids)
      test_idx <- fit_df$deId %in% test_deids
    }
    train_df <- fit_df |> 
      dplyr::filter(!!train_idx)
    test_df <- fit_df |> 
      dplyr::filter(!!test_idx)
    
    formula_long <- stringr::str_replace(
      formula_timevar, "\\.", 
      paste(
        setdiff(colnames(fit_df), c("start_date", "end_date", "event", "deId")),
        collapse = " + "
      )
    )
    
    #### Fit Models ####
    
    # fit/predict LTRC with CIF base estimator
    model_name <- sprintf("%s_ltrc_%s", split_mode, fit_mode)
    model_ls[[model_name]] <- fit_wrapper(
      fit_fun = fit_ltrc,
      formula = formula_long,
      train_df = train_df, 
      test_df = test_df,
      times = eval_times,
      ntree = 100, 
      mtry = round(sqrt(ncol(train_df))),
      cores = n_cores,
      model_name = model_name, results_dir = OUT_DIR, use_cache = USE_CACHED
    )
    if ((fit_mode == "labs_rf_imp") && (split_mode == "test")) {
      model_ls[[model_name]]$fit[["importance"]] <- fit_wrapper(
        fit_fun = permutation_importance,
        fit = model_ls[[model_name]]$fit, 
        data = test_df, 
        times = eval_times,
        formula = formula_long, 
        eval_formula = eval_formula,
        cov_data = cov_df,
        pred_fun = function(fit, x, times) {
          LTRCforests::predictProb(
            fit, newdata = x, newdata.id = deId, time.eval = times
          )$survival.probs |>
            t()
        },
        B = B,
        model_name = sprintf("%s_importance", model_name),
        results_dir = OUT_DIR,
        use_cache = USE_CACHED,
        check_size = FALSE
      )
    }
    keep_models <- c(keep_models, model_name)
    
    # fit/predict cox
    model_name <- sprintf("%s_cox_%s", split_mode, fit_mode)
    model_ls[[model_name]] <- fit_wrapper(
      fit_fun = fit_pec,
      model = survival::coxph, 
      formula = formula_long,
      train_df = train_df,
      test_df = test_df,
      times = eval_times,
      x = TRUE,
      model_name = model_name, results_dir = OUT_DIR, use_cache = USE_CACHED
    )
    keep_models <- c(keep_models, model_name)
    
    # fit/predict cox + lasso
    model_name <- sprintf("%s_cox_lasso_%s", split_mode, fit_mode)
    model_ls[[model_name]] <- fit_wrapper(
      fit_fun = fit_glmnet,
      formula = formula_long,
      train_df = train_df,
      test_df = test_df,
      times = eval_times,
      alpha = 1,
      nfolds = 5,
      model_name = model_name, results_dir = OUT_DIR, use_cache = USE_CACHED
    )
    keep_models <- c(keep_models, model_name)
    
    # fit/predict cox + ridge
    model_name <- sprintf("%s_cox_ridge_%s", split_mode, fit_mode)
    model_ls[[model_name]] <- fit_wrapper(
      fit_fun = fit_glmnet,
      formula = formula_long,
      train_df = train_df,
      test_df = test_df,
      times = eval_times,
      alpha = 0,
      nfolds = 5,
      model_name = model_name, results_dir = OUT_DIR, use_cache = USE_CACHED
    )
    keep_models <- c(keep_models, model_name)
    
    #### Evaluation ####
    cindex_ls[[fit_mode]] <- eval_cindex_tv(
      model_ls[keep_models],
      data = test_df,
      times = eval_times,
      formula = eval_formula
    )
    km_ls[[fit_mode]] <- eval_km_curve(
      model_ls = model_ls[keep_models],
      data = fit_df,
      times = eval_times,
      accruals = accruals
    )
  }
  
  saveRDS(
    fit_df_list,
    file = file.path(
      OUT_DIR, 
      sprintf(
        "%s_fits_forecast_%smo%s.rds", 
        split_mode, forecast_name, censor_switcher_name
      )
    )
  )
  saveRDS(
    cindex_ls,
    file = file.path(
      OUT_DIR, 
      sprintf(
        "%s_cindex_forecast_%smo%s.rds", 
        split_mode, forecast_name, censor_switcher_name
      )
    )
  )
  saveRDS(
    km_ls,
    file = file.path(
      OUT_DIR, 
      sprintf(
        "%s_km_forecast_%smo%s.rds", 
        split_mode, forecast_name, censor_switcher_name
      )
    )
  )
}
