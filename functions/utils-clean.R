#' Create one-hot vector of TRUE/FALSE, allowing for all FALSE when input = NA
#' 
#' @param i Index of TRUE value (or NA if all FALSE)
#' @param length Length of vector
#' 
#' @return Vector of TRUE/FALSE values
indicator_fun <- function(i, length) {
  if (is.na(i)) {
    indic <- rep(FALSE, length)
  } else {
    indic <- Matrix::sparseVector(TRUE, i, length = length) |>
      as.vector()
  }
  return(indic)
}


#' Impute NAs with missForest 
#' 
#' @param data Data frame to impute
#' @param keep_vars Variables to keep in imputed data
#' @param lab_vars Lab variables to impute
#' 
#' @return Imputed data frame
impute_missforest <- function(data, keep_vars = NULL, lab_vars = NULL) {
  if (is.null(keep_vars)) {
    keep_vars <- c(
      "switch_factor", "COMORBIDITY_cat", "Age_IFX", "Gender", "race_gp",
      "Ethnicity", "MS_cat", "priority", "concomitant_med", "TimeYrs_IBD_IFX",
      "IBD_type_c", "rurality", "count_anyenc", "count_IBDenc", 
      "ProviderType_gp", "n_flares_preifx"
    )
  }
  
  data_imputed <- data |>
    dplyr::select(
      tidyselect::any_of(c(keep_vars, lab_vars)), 
      tidyselect::contains("lab")
    ) |>
    as.data.frame() |>
    missForest::missForest()
  data[colnames(data_imputed$ximp)] <- data_imputed$ximp
  return(data)
}
