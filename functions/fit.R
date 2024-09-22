#' Fit model with caching wrapper
#' 
#' @param fit_fun Function to fit model.
#' @param times Vector of evaluation times.
#' @param ... Additional arguments to pass to fit_fun.
#' @param model_name Name of model.
#' @param results_dir Directory to save model results.
#' @param use_cache Logical indicating whether to use cached results.
#' @param check_size Logical indicating whether to check the size of the
#'   cached results. If sizes (i.e., number of samples) do not match, the
#'   model will be refit.
#'   
#' @returns Resulting model fit (i.e., output of fit_fun()).
fit_wrapper <- function(fit_fun, times, ..., model_name, results_dir, 
                        use_cache = TRUE, check_size = TRUE) {
  model_fname <- file.path(results_dir, sprintf("%s.rds", model_name))
  if (file.exists(model_fname) && use_cache) {
    results <- readRDS(model_fname)
    if (check_size) {
      dots_ls <- rlang::dots_list(...)
      test_df <- dots_ls[["test_df"]]
      test_predictions <- get_predictions_list(list(results))
      if ((ncol(test_predictions[[1]]) == length(times)) &&
          (nrow(test_predictions[[1]]) == length(unique(test_df$deId)))) {
        return(results)
      }
    } else {
      return(results)
    }
  }
  print(sprintf("Fitting %s", model_fname))
  results <- fit_fun(..., times = times)
  saveRDS(results, file = model_fname)
  return(results)
}

#' Fit pec-compatible model
#' 
#' @param model Function to fit model.
#' @param formula Formula for model.
#' @param train_df Training data.
#' @param test_df Test data.
#' @param times Vector of evaluation times.
#' @param ... Additional arguments to pass to model.
#' 
#' @returns List containing model fit and predictions.
fit_pec <- function(model, formula, train_df, test_df, times, ...) {
  fit <- model(formula = as.formula(formula), data = train_df, ...)
  test_preds <- pec::predictSurvProb(fit, test_df, times)
  if (any(is.na(test_preds))) {
    test_preds[is.na(test_preds)] <- 0
  }
  if (length(as.character(as.formula(formula)[[2]])) == 4) {
    test_preds <- aggregate_timevar_predictions(
      preds = test_preds, data = test_df, times = times
    )
    pec_flag <- FALSE
  } else {
    pec_flag <- TRUE
  }
  out <- list(
    fit = fit,
    test_preds = test_preds,
    pec = pec_flag
  )
  return(out)
}


#' Fit time-varying survival forest using LTRCforests
#' 
#' @param formula Formula for model.
#' @param train_df Training data.
#' @param test_df Test data.
#' @param times Vector of evaluation times.
#' @param ... Additional arguments to pass to LTRCforests::ltrccif.
#' 
#' @returns List containing model fit and predictions.
fit_ltrc <- function(formula, train_df, test_df, times, ...) {
  fit <- LTRCforests::ltrccif(
    formula = as.formula(formula), 
    data = as.data.frame(train_df), 
    id = deId, 
    ...
  )
  test_preds <- LTRCforests::predictProb(
    fit, newdata = test_df, newdata.id = deId, time.eval = times
  )$survival.probs |>
    t()
  rownames(test_preds) <- unique(test_df$deId)
  out <- list(
    fit = fit,
    test_preds = test_preds,
    pec = FALSE
  )
  return(out)
}


#' Fit survival forest using randomForestSRC
#' 
#' @param formula Formula for model.
#' @param train_df Training data.
#' @param test_df Test data.
#' @param ... Additional arguments to pass to randomForestSRC::rfsrc.
#' 
#' @returns List containing model fit and predictions.
fit_rfsrc <- function(formula, train_df, test_df, ...) {
  fit <- randomForestSRC::rfsrc(
    formula = as.formula(formula), 
    data = as.data.frame(train_df),
    ...
  )
  test_preds <- predict(fit, as.data.frame(test_df))$survival
  out <- list(
    fit = fit,
    test_preds = test_preds,
    pec = TRUE
  )
  return(out)
}


#' Fit regularized cox model using glmnet
#' 
#' @param formula Formula for model.
#' @param train_df Training data.
#' @param test_df Test data.
#' @param times Vector of evaluation times.
#' @param type.measure Type of measure to use for cv.glmnet.
#' @param best_lambda Best lambda to use for predictions.
#' @param ... Additional arguments to pass to glmnet::cv.glmnet.
#' 
#' @returns List containing model fit and predictions.
fit_glmnet <- function(formula, train_df, test_df, times, 
                       type.measure = "C", best_lambda = "lambda.min", ...) {
  surv_obj_str <- as.character(as.formula(formula)[[2]])
  if (length(surv_obj_str) == 3) {
    y_surv_obj <- Surv(
      train_df[[surv_obj_str[2]]],
      train_df[[surv_obj_str[3]]]
    )
  } else {
    y_surv_obj <- Surv(
      train_df[[surv_obj_str[2]]],
      train_df[[surv_obj_str[3]]],
      train_df[[surv_obj_str[4]]]
    )
  }
  # convert data frame to centered+scaled numeric matrix for glmnet input
  x <- train_df |>
    dplyr::select(-tidyselect::all_of(c(surv_obj_str[-1], "deId")))
  dummy_var_transform <- caret::dummyVars(~ ., data = x, fullRank = TRUE)
  x <- predict(dummy_var_transform, x)
  xtest <- predict(dummy_var_transform, test_df)
  x_mean <- colMeans(x)
  x_sd <- apply(x, 2, sd)
  x <- scale(x, center = x_mean, scale = x_sd)
  xtest <- scale(xtest, center = x_mean, scale = x_sd)
  
  fit <- glmnet::cv.glmnet(
    x = x, y = y_surv_obj, family = "cox", type.measure = type.measure, ...
  )
  test_surv_out <- survival::survfit(
    fit, s = fit[[best_lambda]], x = x, y = y_surv_obj, newx = xtest
  )
  test_preds <- t(summary(test_surv_out, times = times)$surv)

  if (any(is.na(test_preds))) {
    test_preds[is.na(test_preds)] <- 0
  }
  if (length(as.character(as.formula(formula)[[2]])) == 4) {
    test_preds <- aggregate_timevar_predictions(
      preds = test_preds, data = test_df, times = times
    )
  }
  
  out <- list(
    fit = fit,
    test_preds = test_preds,
    pec = FALSE
  )
  return(out)
}


#' Evaluate permutation importance
#' 
#' @param fit Model fit.
#' @param data Data frame used to permute time-varying covariates.
#' @param times Evaluation times.
#' @param formula Formula for model.
#' @param eval_formula Formula for evaluation.
#' @param cov_data Data frame used to permute non-time-varying covariates.
#'   Should have one row per patient.
#' @param pred_fun Prediction function.
#' @param B Number of permutations.
#' 
#' @returns List containing feature name and feature importance, as measured
#'   by the change in c-index.
permutation_importance <- function(fit, data, times, 
                                   formula, eval_formula,
                                   cov_data = NULL, 
                                   pred_fun = predict, B = 100) {
  
  if (is.null(eval_formula)) {
    eval_formula <- Surv(event_censor_date, event) ~ 1
  }
  times <- unname(times)
  
  n_samples <- length(unique(data$deId))
  if (!is.null(cov_data)) {
    cov_data <- cov_data |>
      dplyr::filter(deId %in% unique(data$deId))
  }
  # permutation indices for non-time-varying covariates
  permute_idxs <- purrr::map(
    1:B, 
    ~ rep(
      sample(1:n_samples, n_samples, replace = FALSE), 
      times = rle(data$deId)$lengths
    )
  )
  # permutation indices for time-varying covariates
  tpermute_idxs <- purrr::map(
    1:B, ~ sample(1:nrow(data), nrow(data), replace = FALSE)
  )
  
  # original predictions and cindex
  preds_orig <- pred_fun(fit, data, times)
  eval_data <- data |>
    dplyr::group_by(deId) |>
    dplyr::mutate(event_censor_date = max(end_date)) |>
    dplyr::arrange(-start_date)
  cindex_orig <- purrr::imap_dbl(
    times,
    function(t, idx) {
      time_df <- eval_data |>
        dplyr::filter(start_date <= !!t) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup(deId)
      pec::cindex(
        purrr::map(list(preds_orig), ~ .x[, idx, drop = FALSE]),
        formula = eval_formula,
        data = time_df,
        eval.times = t,
        splitMethod = "none"
      )[[1]][[1]]
    }
  )
  
  var_names <- labels(terms(as.formula(formula)))
  cindex_perm <- purrr::map(
    var_names,
    function(j) {
      if (j %in% colnames(cov_data)) {
        p_idxs <- permute_idxs
        p_data <- cov_data
      } else {
        p_idxs <- tpermute_idxs
        p_data <- data
      }
        
      cindex_ls <- furrr::future_map(
        p_idxs,
        function(p_idx) {
          perm_data <- data
          perm_data[, j] <- p_data[p_idx, j]
          preds <- pred_fun(fit, perm_data, times)
          scores <- purrr::imap_dbl(
            times,
            function(t, idx) {
              time_df <- eval_data |>
                dplyr::filter(start_date <= !!t) |>
                dplyr::slice_head(n = 1) |>
                dplyr::ungroup(deId)
              pec::cindex(
                purrr::map(list(preds), ~ .x[, idx, drop = FALSE]),
                formula = eval_formula,
                data = time_df,
                eval.times = t,
                splitMethod = "none"
              )[[1]][[1]]
            }
          )
        },
        .options = furrr::furrr_options(
          seed = 123, 
          packages = c("prodlim", "survival")
        )
      )
      purrr::reduce(cindex_ls, `+`) / B
    }
  ) |>
    setNames(var_names)
  
  fi_df <- purrr::map(
    cindex_perm,
    function(cindex_p) {
      tibble::tibble(importance = mean(cindex_orig - cindex_p, na.rm = TRUE))
    }
  ) |>
    purrr::list_rbind(names_to = "variable")
  
  out <- list(
    importance = fi_df,
    cindex_orig = cindex_orig,
    cindex_perm = cindex_perm
  )
  return(out)
}
