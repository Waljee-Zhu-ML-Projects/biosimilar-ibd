rm(list = ls())
library(ggfortify)

# initialize directory variables and load in helper functions
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
seeds <- 1:50
model_names <- c("ltrc", "cox", "cox_ridge", "cox_lasso")
varset_modes <- c("no_labs", "labs_mean_imp", "labs_median_imp", "labs_rf_imp")

# aggregate survival prediction results
for (split_mode in c("validation", "test")) {
  cindex_ls <- list()
  km_ls <- list()
  for (censor_switcher_i in 1:length(censor_switchers)) {
    censor_switcher <- censor_switchers[censor_switcher_i]
    censor_switcher_name <- names(censor_switchers)[censor_switcher_i]
    for (forecast_i in 1:length(forecasts)) {
      forecast <- forecasts[forecast_i]
      forecast_name <- names(forecasts)[forecast_i]
      for (seed in seeds) {
        OUT_DIR <- file.path(
          RESULTS_DIR, sprintf("seed%s", seed), sprintf("forecast%s", forecast_name)
        )
        
        cindex_df <- readRDS(
          file.path(
            OUT_DIR, 
            sprintf(
              "%s_cindex_forecast_%smo%s.rds", 
              split_mode, forecast_name, censor_switcher_name
            )
          )
        ) |>
          purrr::imap(
            function(x, varset_mode) {
              df <- dplyr::bind_cols(time = x$time, x$cindex)
              colnames(df) <- stringr::str_remove(
                colnames(df), sprintf("_%s", varset_mode)
              ) |>
                stringr::str_remove(sprintf("%s_", split_mode))
              return(df)
            }
          ) |>
          dplyr::bind_rows(.id = "varset_mode") |>
          dplyr::mutate(
            split_mode = !!split_mode,
            rep = !!seed,
            forecast = !!forecast_name,
            censor_switcher = !!censor_switcher_name == "_censor_switcher"
          )
        cindex_ls <- c(cindex_ls, list(cindex_df))
        
        km_df <- readRDS(
          file.path(
            OUT_DIR, 
            sprintf(
              "%s_km_forecast_%smo%s.rds", 
              split_mode, forecast_name, censor_switcher_name
            )
          )
        ) |>
          purrr::imap(
            function(x, varset_mode) {
              purrr::imap(
                x,
                function(xx, accrual) {
                  purrr::imap(
                    xx,
                    function(xxx, model_name) {
                      df <- ggplot2::autoplot(xxx)$data |>
                        dplyr::mutate(
                          method = stringr::str_remove(
                            model_name, sprintf("_%s", varset_mode)
                          ) |>
                            stringr::str_remove(sprintf("%s_", split_mode))
                        )
                    }
                  ) |>
                    dplyr::bind_rows()
                }
              ) |>
                dplyr::bind_rows(.id = "accrual")
            }
          ) |>
          dplyr::bind_rows(.id = "varset_mode") |>
          dplyr::mutate(
            split_mode = !!split_mode,
            rep = !!seed,
            forecast = !!forecast_name,
            censor_switcher = !!censor_switcher_name == "_censor_switcher"
          ) |>
          dplyr::filter(
            varset_mode == "labs_rf_imp"
          )
        km_ls <- c(km_ls, list(km_df))
      }
    }
  }
  
  cindex_df <- dplyr::bind_rows(cindex_ls) |>
    dplyr::mutate(
      forecast_name = factor(
        sprintf("Forecast = %smo", forecast),
        levels = sprintf("Forecast = %smo", names(forecasts))
      )
    )
  km_df <- dplyr::bind_rows(km_ls) |>
    dplyr::mutate(
      accrual_name = factor(
        sprintf("Accrual = %syr", accrual),
        levels = sprintf("Accrual = %syr", names(accruals))
      ),
      forecast_name = factor(
        sprintf("Forecast = %smo", forecast),
        levels = sprintf("Forecast = %smo", names(forecasts))
      )
    ) |>
    tibble::as_tibble()
  
  saveRDS(
    cindex_df,
    file = file.path(RESULTS_DIR, sprintf("%s_cindexs.rds", split_mode))
  )
  saveRDS(
    km_df,
    file = file.path(RESULTS_DIR, sprintf("%s_km.rds", split_mode))
  )
}

# aggregate survival importance results
vimp_ls <- list()
split_mode <- "test"
for (censor_switcher_i in 1:length(censor_switchers)) {
  censor_switcher <- censor_switchers[censor_switcher_i]
  censor_switcher_name <- names(censor_switchers)[censor_switcher_i]
  for (forecast_i in 1:length(forecasts)) {
    forecast <- forecasts[forecast_i]
    forecast_name <- names(forecasts)[forecast_i]
    for (seed in seeds) {
      OUT_DIR <- file.path(
        RESULTS_DIR, 
        sprintf("seed%s", seed), 
        sprintf("forecast%s", forecast_name)
      )
      for (model_name in model_names) {
        for (varset_mode in varset_modes) {
          if (stringr::str_detect(model_name, "cox")) {
            fit_out <- readRDS(
              file.path(
                OUT_DIR, 
                sprintf("%s_%s_%s.rds", split_mode, model_name, varset_mode)
              )
            )
            fit <- fit_out$fit
            if (identical(model_name, "cox")) {
              vimp_df <- broom::tidy(fit) |>
                dplyr::rename(variable = term) |>
                dplyr::mutate(
                  HR = exp(estimate),
                  importance = abs(statistic)
                )
            } else {
              best_lambda_idx <- which(fit$lambda == fit[["lambda.min"]])
              best_coefs <- as.matrix(fit$glmnet.fit$beta)[, best_lambda_idx]
              vimp_df <- tibble::tibble(
                variable = names(best_coefs),
                estimate = best_coefs,
                importance = abs(best_coefs),
                HR = exp(best_coefs)
              )
            }
          } else {
            importance_fname <- file.path(
              OUT_DIR, 
              sprintf(
                "%s_%s_%s_importance.rds", split_mode, model_name, varset_mode
              )
            )
            if (file.exists(importance_fname)) {
              fit_out <- readRDS(importance_fname)
              vimp_df <- fit_out$importance
            } else {
              next
            }
          }
          vimp_df <- vimp_df |>
            dplyr::mutate(
              rep = !!seed,
              varset_mode = !!varset_mode,
              forecast = !!forecast_name,
              censor_switcher = !!censor_switcher_name == "_censor_switcher",
              method = !!model_name
            ) |>
            dplyr::arrange(-importance)
          vimp_ls <- c(vimp_ls, list(vimp_df))
        }
      }
    }
  }
}

vimp_df <- dplyr::bind_rows(vimp_ls) |>
  dplyr::mutate(
    forecast_name = factor(
      sprintf("Forecast = %smo", forecast),
      levels = sprintf("Forecast = %smo", names(forecasts))
    )
  )
saveRDS(
  vimp_df,
  file = file.path(RESULTS_DIR, sprintf("%s_importances.rds", split_mode))
)
