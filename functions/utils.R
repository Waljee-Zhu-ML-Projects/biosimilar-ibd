#' Return list of models object to pass to pec::pec or pec::cindex
#' 
#' @param model_ls List of models
#' 
#' @returns Returns a list where each element is a model object to be passed
#'   to pec::pec if model$pec is TRUE, otherwise the valid/test predictions 
#'   matrix to be passed into pec::cindex.
get_pec_model_list <- function(model_ls) {
  pec_model_ls <- purrr::map(
    model_ls,
    function(model) {
      if (model$pec) {
        model$fit
      } else {
        model$test_preds
      }
    }
  )
  return(pec_model_ls)
}


#' Return list of predictions from survival models
#' 
#' @param model_ls List of models
#' 
#' @returns List of predictions
get_predictions_list <- function(model_ls) {
  preds_ls <- purrr::map(model_ls, "test_preds")
  return(preds_ls)
}


#' Aggregate survival predictions for time-varying models
#' 
#' @param preds Predictions matrix
#' @param data Data frame with deId, start_date, and end_date columns
#' @param times Evaluation times
#' 
#' @returns Matrix of size n_samples x length(times) with survival predictions
aggregate_timevar_predictions <- function(preds, data, times) {
  preds_ls <- purrr::map(
    unique(data$deId),
    function(id) {
      keep_idxs <- data$deId == id
      id_valid_df <- data |>
        dplyr::filter(!!keep_idxs) |>
        dplyr::arrange(start_date)
      id_preds_mat <- preds[keep_idxs, , drop = FALSE]
      id_preds <- c()
      for (i in 1:nrow(id_valid_df)) {
        start_date <- id_valid_df$start_date[i]
        if (i == nrow(id_valid_df)) {
          end_date <- max(eval_times) + 1
        } else {
          end_date <- id_valid_df$end_date[i]
        }
        in_times <- which((eval_times >= start_date) & (eval_times < end_date))
        if (length(in_times) > 0) {
          id_preds <- c(id_preds, id_preds_mat[i, in_times])
        }
      }
      return(id_preds)
    }
  )
  preds_df <- do.call(rbind, preds_ls)
  rownames(preds_df) <- unique(data$deId)
  return(preds_df)
}
