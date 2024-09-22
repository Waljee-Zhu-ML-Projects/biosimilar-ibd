#' Evaluate C-index for time-varying survival models
#' 
#' @param model_ls List of models
#' @param data Time-varying (validation) data frame to evaluate c-index
#' @param times Evaluation times
#' @param formula Formula for survival model
#' 
#' @returns List of C-index values and evaluation times
eval_cindex_tv <- function(model_ls, data, times, formula) {
  formula <- as.formula(formula)
  times <- unname(times)
  
  pec_model_ls <- get_pec_model_list(model_ls)
  preds_ls <- get_predictions_list(model_ls)
  
  eval_data <- data |>
    dplyr::group_by(deId) |>
    dplyr::mutate(event_censor_date = max(end_date)) |>
    dplyr::arrange(-start_date)
  out_times <- purrr::imap(
    times,
    function(eval_time, idx) {
      time_df <- eval_data |>
        dplyr::filter(start_date <= !!eval_time) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup(deId)
      pec::cindex(
        purrr::map(pec_model_ls, ~ .x[, idx, drop = FALSE]),
        formula = formula,
        data = time_df,
        eval.times = eval_time,
        splitMethod = "none"
      )
    }
  )
  method_names <- names(out_times[[1]][[1]])
  out <- list(
    cindex = purrr::map(
      method_names, 
      function(method) {
        purrr::map_dbl(out_times, ~ .x[[1]][[method]])
      }
    ) |>
      setNames(method_names),
    time = times
  )
  return(out)
}


#' Evaluate Kaplan-Meier curve for time-varying survival models
#' 
#' @param model_ls List of models
#' @param data Data frame with observed event and event times
#' @param times Evaluation times
#' @param forecast Forecast time
#' @param accruals Accrual times
#' 
#' @returns List of Kaplan-Meier curves for different accrual times and models.
eval_km_curve <- function(model_ls, data, times, accruals) {
  preds_ls <- get_predictions_list(model_ls)
  preds_df <- purrr::map(
    preds_ls,
    ~ as.data.frame(.x) |>
      setNames(paste0("Accrual = ", times)) |>
      tibble::rownames_to_column("deId") |>
      dplyr::left_join(
        data |>
          dplyr::group_by(deId) |>
          dplyr::slice_max(end_date, n = 1) |>
          dplyr::select(deId, event, event_censor_date = end_date),
        by = "deId"
      )
  ) |>
    dplyr::bind_rows(.id = "method")
  
  km_out_ls <- purrr::map(
    accruals,
    function(eval_time) {
      purrr::map(
        names(preds_ls),
        function(method_name) {
          plt_df <- preds_df |>
            dplyr::filter(method == !!method_name) |>
            dplyr::mutate(
              surv_prob = .data[[paste0("Accrual = ", eval_time)]],
              risk = dplyr::case_when(
                surv_prob > quantile(surv_prob, 0.5) ~ "low",
                TRUE ~ "high"
              ) |>
                factor(levels = c("high", "low")),
              event_censor_date = event_censor_date / 365.25
            )
          
          km_fit <- survival::survfit(
            Surv(event_censor_date, event) ~ risk, data = plt_df
          )
        }
      ) |>
        setNames(names(preds_ls))
    }
  )
  
  return(km_out_ls)
}
