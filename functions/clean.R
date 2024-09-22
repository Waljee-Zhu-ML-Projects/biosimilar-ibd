#' Clean patient data
#' 
#' @param patient_data_orig Original (uncleaned) patient data.
#' @param flares_data A tibble containing the (cleaned) flares data. Only used
#'   if `mode` is not "default".
#' @param accrual Length of the accrual period (in days). Only used if `mode`
#'   is not "default".
#' @param forecast Length of the forecast period (in days). Only used if `mode`
#'   is not "default".
#' @param censor_switcher Logical indicating whether to censor patients who
#'   switched within the forecast window. Only used if `mode` is not "default".
#' @param mode Character string indicating the cleaning mode. Must be
#'   one of "default", "classification", or "survival".
#'   
#' @return A tibble containing the cleaned patient data.
clean_patient_data <- function(patient_data_orig, 
                               flares_data = NULL,
                               accrual = NULL, 
                               forecast = NULL,
                               censor_switcher = FALSE,
                               mode = c("default", "classification", "survival")) {
  mode <- match.arg(mode)
  patient_data <- patient_data_orig |>
    dplyr::select(
      deId, deIdSta3n, 
      switch, newindexswitch_date, 
      ibd_max, n_flares_preifx,
      IBDtype_max,
      Gender, Race, Ethnicity, MaritalStatus,
      priority = EnrollmentPriority, 
      concomitant_med, steroids, steroids_no_dx,
      rurality = GISURH,
      count_anyenc, count_IBDenc, 
      ProviderType,
      Date_Birth, Date_Death, 
      COMORBIDITY, 
      TimeYrs_IBD_IFX, Age_IFX, 
      days_to_death_or_end2021, 
      yrs_to_death_or_end2021
    ) |>
    dplyr::mutate(
      race_gp = dplyr::case_when(
        Race == "BLACK OR AFRICAN AMERICAN" ~ "African Am",
        Race %in% c("WHITE", "WHITE NOT OF HISP ORIG") ~ "White",
        Race == "" ~ "Unknown",
        TRUE ~ "Other"
      ) |>
        factor(levels = c("White", "African Am", "Other", "Unknown")),
      Ethnicity = dplyr::case_when(
        !(Ethnicity %in% c("HISPANIC OR LATINO", "NOT HISPANIC OR LATINO")) ~
          "UNKNOWN",
        TRUE ~ Ethnicity
      ) |>
        factor(
          levels = c("HISPANIC OR LATINO", "NOT HISPANIC OR LATINO", "UNKNOWN")
        ),
      MS_cat = dplyr::case_when(
        MaritalStatus == "MARRIED" ~ "MARRIED",
        MaritalStatus %in% c("NEVER MARRIED", "SINGLE") ~ "SINGLE",
        TRUE ~ "OTHER"
      ) |>
        factor(levels = c("MARRIED", "OTHER", "SINGLE")),
      priority = dplyr::case_when(
        as.numeric(stringr::str_extract_all(priority, "[0-9]")) <= 6 ~ "low",
        TRUE ~ "high"
      ) |>
        factor(
          levels = c("low", "high")
        ),
      COMORBIDITY_cat = dplyr::case_when(
        COMORBIDITY == 0 ~ "None",
        COMORBIDITY %in% c(1, 2) ~ "1-2",
        COMORBIDITY >= 3 ~ "3+"
      ) |>
        factor(levels = c("None", "1-2", "3+")),
      IBD_type_c = dplyr::case_when(
        IBDtype_max == "CR" ~ "CR",
        TRUE ~ "UC/ID"
      ) |>
        factor(levels = c("CR", "UC/ID")),
      rurality = dplyr::case_when(
        rurality == "U" ~ "Urban",
        TRUE ~ "Rural"
      ) |>
        factor(levels = c("Urban", "Rural")),
      ProviderType_gp = dplyr::case_when(
        ProviderType == "Allopathic & Osteopathic Physicians" ~ 
          "Allopathic & Osteopathic Physicians",
        TRUE ~ "Physician Assistants & Nursing Prov"
      ) |>
        factor(levels = c("Allopathic & Osteopathic Physicians",
                          "Physician Assistants & Nursing Prov")),
      rawcount_anyenc = count_anyenc,
      rawcount_IBDenc = count_IBDenc,
      count_anyenc = log(count_anyenc + 1),
      count_IBDenc = log(count_IBDenc + 1)
    )
  
  if (!identical(mode, "default")) {
    if (is.null(forecast) || is.null(flares_data)) {
      stop(
        sprintf(
          "forecast and flares_data must be provided when mode == '%s'", mode
        )
      )
    }
    
    orig_patient_data_cols <- colnames(patient_data)
    patient_data <- patient_data |>
      dplyr::left_join(flares_data, by = "deId") |>
      dplyr::mutate(
        censor_date = pmin(Date_Death, days_to_death_or_end2021, na.rm = TRUE),
        event_censor_date = pmin(censor_date, event_date, na.rm = TRUE),
      )
    if (identical(mode, "classification")) {
      # format data to predict flares within forecast period
      rand_times <- runif(nrow(patient_data), 0, forecast)
      patient_data <- patient_data |>
        dplyr::mutate(
          rand_time = !!rand_times,
          orig_event_censor_date = event_censor_date,
          event_censor_date = rand_time + purrr::map2_dbl(
            event, event_censor_date, 
            function(.x, .y) {
              if (isTRUE(.x)) {
                .y
              } else {
                runif(1, forecast + accrual, .y)
              }
            }
          )
        )
    }
    # update switch indicator (in case switch happened after event)
    patient_data <- patient_data |>
      dplyr::mutate(
        switch = dplyr::case_when(
          newindexswitch_date <= (event_censor_date - forecast) ~ 1,
          newindexswitch_date <= event_censor_date ~ NA,
          TRUE ~ 0
        )
      )
    if (identical(mode, "survival")) {
      patient_data <- patient_data |>
        dplyr::select(tidyselect::all_of(orig_patient_data_cols))
    } else if (identical(mode, "classification")) {
      patient_data <- patient_data |>
        dplyr::select(
          tidyselect::all_of(orig_patient_data_cols),
          rand_time, event_censor_date, orig_event_censor_date
        )
    }
    
    patient_data <- patient_data |>
      dplyr::mutate(
        newindexswitch_date = dplyr::case_when(
          switch == 0 ~ NA,
          TRUE ~ newindexswitch_date
        )
      )
    if (censor_switcher) {
      # remove patients who switched within the forecast window
      patient_data <- patient_data |>
        dplyr::filter(!is.na(switch))
    } else {
      # don't remove patients who switched within the forecast window
      patient_data <- patient_data |>
        dplyr::mutate(
          switch = tidyr::replace_na(switch, 0)
        )
    }
  }
  
  patient_data <- patient_data |>
    dplyr::mutate(
      switch_factor = factor(switch)
    )
  
  return(patient_data)
}


#' Clean flares data
#' 
#' @param flares_data_orig Original (uncleaned) flares data.
#' @param accrual Length of the accrual period (in days). Only used if `mode`
#'   is "classification" or "survival".
#' @param forecast Length of the forecast period (in days). Only used if `mode`
#'   is "classification" or "survival".
#' @param mode Character string indicating which flares to keep, whether only 
#'   the "first" flare per patient, "all" flares per patient, the flares used
#'   for the "classification" problem or for the "survival" problem.
#' 
#' @return A tibble containing the cleaned flares data.
clean_flares_data <- function(flares_data_orig, 
                              accrual = NULL, 
                              forecast = NULL,
                              mode = c("first", 
                                       "classification",
                                       "survival", 
                                       "all")) {
  mode <- match.arg(mode)
  flares_data <- flares_data_orig |>
    dplyr::select(
      deId, event_date, event_type, disease_type, first
    )
  if (identical(mode, "first")) {
    flares_data <- flares_data |>
      dplyr::filter(first)
  } else if (!identical(mode, "all")) {
    if (identical(mode, "classification")) {
      if (is.null(accrual) || is.null(forecast)) {
        stop(
          sprintf("accrual and forecast must be provided when mode == '%s'", mode)
        )
      }
      min_event_date <- accrual + forecast
      max_prior_flares_date <- accrual
    } else if (identical(mode, "survival")) {
      if (is.null(forecast)) {
        stop(
          sprintf("forecast must be provided when mode == '%s'", mode)
        )
      }
      min_event_date <- forecast
      max_prior_flares_date <- forecast
    }
    flares_data <- flares_data |>
      dplyr::group_by(deId) |>
      dplyr::mutate(
        # get flares within the accrual period but after IFX index
        prior_flare = event_date <= !!max_prior_flares_date,
        # get first flare after accrual period
        first = indicator_fun(
          which(event_date > !!min_event_date)[1], length = dplyr::n()
        )
      ) |>
      dplyr::summarise(
        # number of flares within the accrual period but after IFX index
        n_prior_flare = sum(prior_flare),
        prior_event_dates = list(event_date[prior_flare]),
        # whether the patient had any flare after the accrual period
        event = any(first),
        event_date = sum(event_date * first)
      ) |>
      dplyr::mutate(
        event_date = ifelse(event, event_date, NA)
      )
  }
  
  if (mode == "all") {
    flares_data <- flares_data
  } else if (mode == "forecast_class") {
    if (is.null(accrual) || is.null(forecast)) {
      stop("accrual and forecast must be provided when mode == 'forecast_class'")
    }
    flares_data <- flares_data |>
      dplyr::group_by(deId) |>
      dplyr::summarise(
        n_prior_flare = sum(event_date <= accrual),
        event = any(
          (event_date > accrual) & (event_date <= (accrual + forecast))
        )
      ) |>
      dplyr::ungroup()
  } else {
    if (mode == "forecast_class_rolling") {
      min_event_date <- accrual + forecast
      max_prior_flares_date <- accrual
    } else if (mode == "forecast_survival") {
      if (is.null(forecast)) {
        stop("forecast must be provided when mode == 'forecast_survival'")
      }
      min_event_date <- forecast
      max_prior_flares_date <- forecast
    }
    
  }
  return(flares_data)
}


#' Clean lab data
#' 
#' @param lab_data_orig Original (uncleaned) lab data.
#' @param patient_data A tibble containing the cleaned patient data.
#' @param accrual Length of the accrual period (in days). Only used if
#'   `summarize` is TRUE and `index_mode` is "event".
#' @param forecast Length of the forecast period (in days). Only used if
#'   `summarize` is TRUE and `index_mode` is "event".
#' @param summarize Logical indicating whether to summarize lab data by 
#'   computing their mean, max, SD, total variation, slope, and most recent lab 
#'   values.
#' @param index_mode Character string indicating when to start accruing labs
#'   to summarize. Must be one of "IFX" or "event". Only used if `summarize` is
#'   TRUE.
#'   - "IFX": start accruing labs at the IFX index date until event date
#'   - "event": start accruing labs within the accrual period prior to 
#'     the (event date - forecast period) date
#' @param impute_mode Character string indicating how to impute missing lab
#'   values. Must be one of "median", "mean", or "rf".
#' 
#' @return A tibble containing the cleaned lab data.
clean_lab_data <- function(lab_data_orig, 
                           patient_data,
                           accrual = NULL, 
                           forecast = NULL,
                           summarize = TRUE,
                           index_mode = c("IFX", "event"),
                           impute_mode = c("median", "mean", "rf")) {
  index_mode <- match.arg(index_mode)
  impute_mode <- match.arg(impute_mode)
  
  # log-transform some lab variables
  no_transform_vars <- c(
    "co2", "mpv", "basophils", "hematocrit", "hgb", "mch", "mchc", "rbc"
  )
  lab_data <- lab_data_orig |>
    dplyr::mutate(
      LabValue = dplyr::case_when(
        LabName %in% no_transform_vars ~ LabValue,
        TRUE ~ log(LabValue + 1)
      )
    ) |>
    dplyr::left_join(
      patient_data |> 
        dplyr::select(
          deId, newindexswitch_date, tidyselect::any_of("event_censor_date")
        ), 
      by = "deId"
    ) |>
    dplyr::arrange(time_order_within_patient)
  
  if (summarize) {
    if (identical(index_mode, "IFX")) {
      lab_data <- lab_data |>
        # dplyr::filter(
        #   LabDate <= event_censor_date
        # ) |>
        dplyr::mutate(
          period = dplyr::case_when(
            LabDate <= 0 ~ "preindex",
            TRUE ~ "accrual"
          )
        )
    } else if (identical(index_mode, "event")) {
      if (is.null(forecast) || is.null(accrual)) {
        stop("forecast and accrual must be provided when index_mode = 'event'")
      }
      lab_data <- lab_data |>
        dplyr::filter(
          LabDate <= (event_censor_date - !!forecast)
        ) |>
        dplyr::mutate(
          period = dplyr::case_when(
            LabDate >= (event_censor_date - !!forecast - !!accrual) ~ "accrual",
            TRUE ~ "preindex"
          )
        )
    }
    
    recent_lab_data <- lab_data |>
      dplyr::group_by(LabName, deId) |>
      dplyr::summarise(
        recent_lab_value = LabValue[dplyr::n()],
        .groups = "drop"
      )
    if (identical(index_mode, "event")) {
      lab_data <- lab_data |>
        dplyr::filter(period == "accrual")
    }
    
    get_ols_slope <- function(x, y, recent_times = NULL) {
      if (length(x) == 1) {
        return(NA)
      } else {
        if (!is.null(recent_times)) {
          y <- rev(y)[1:min(length(y), recent_times)]
          x <- rev(x)[1:min(length(x), recent_times)]
        }
        return(lm.fit(cbind(1, x), y)$coefficients[["x"]])
      }
    }
    lab_data <- lab_data |>
      dplyr::group_by(LabName, deId, period) |>
      dplyr::summarise(
        n_labs = dplyr::n(),
        mean_lab_value = mean(LabValue),
        max_lab_value = max(LabValue),
        sd_lab_value = sd(LabValue),
        tv_lab_value = mean(abs(LabValue[-1] - LabValue[-dplyr::n()])),
        ols_slope_lab_value = get_ols_slope(LabDate, LabValue),
        recent_ols_slope_lab_value = get_ols_slope(LabDate, LabValue, 3),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        id_cols = c(LabName, deId), 
        names_from = period,
        values_from = -c(LabName, deId, period)
      )
    
    if (identical(index_mode, "event")) {
      # impute mean/max labs with most recent if no labs within the accrual period
      lab_data <- dplyr::full_join(
        lab_data, recent_lab_data, by = c("deId", "LabName")
      ) |>
        dplyr::mutate(
          dplyr::across(
            c(mean_lab_value_accrual, max_lab_value_accrual),
            ~ ifelse(is.na(.x), recent_lab_value, .x)
          )
        ) |>
        dplyr::select(-recent_lab_value)
    }
    
    # get all possible (id, lab) combinations
    full_group_ids <- expand.grid(
      deId = unique(patient_data$deId), 
      LabName = unique(lab_data$LabName)
    ) |>
      tibble::as_tibble() |> 
      dplyr::mutate(
        dplyr::across(tidyselect::everything(), as.character)
      )
    # get lab data for all possible (id, lab) combinations
    lab_data <- dplyr::left_join(
      full_group_ids, lab_data, by = c("deId", "LabName")
    ) |>
      # add most recent lab values
      dplyr::left_join(recent_lab_data, by = c("deId", "LabName")) |>
      dplyr::mutate(
        # impute slope and # labs with 0 if NA (i.e., 0 labs or NA slope since only 1 lab value)
        dplyr::across(
          c(tidyselect::contains("slope"), tidyselect::starts_with("n_labs")),
          ~ tidyr::replace_na(.x, 0)
        )
      ) |>
      tidyr::pivot_wider(
        id_cols = deId,
        names_from = LabName,
        values_from = -c(LabName, deId)
      )
    
    # remove lab columns with high rates of missingness
    rm_na_cols <- colnames(lab_data)[
      stringr::str_ends(colnames(lab_data), "FeCal|crp|esr")
    ]
    lab_data <- lab_data |>
      dplyr::select(-tidyselect::all_of(rm_na_cols))
    
    # impute mean, max, sd, tv, recent columns
    if (impute_mode == "median") {
      lab_data <- lab_data |>
        dplyr::mutate(
          dplyr::across(
            c(
              tidyselect::starts_with("mean"),
              tidyselect::starts_with("max"),
              tidyselect::starts_with("sd"),
              tidyselect::starts_with("tv"),
              tidyselect::starts_with("recent")
            ), 
            ~ tidyr::replace_na(.x, median(.x, na.rm = TRUE))
          )
        )
    } else if (impute_mode == "mean") {
      lab_data <- lab_data |>
        dplyr::mutate(
          dplyr::across(
            c(
              tidyselect::starts_with("mean"),
              tidyselect::starts_with("max"),
              tidyselect::starts_with("sd"),
              tidyselect::starts_with("tv"),
              tidyselect::starts_with("recent")
            ), 
            ~ tidyr::replace_na(.x, mean(.x, na.rm = TRUE))
          )
        )
    } else if (impute_mode == "rf") {
      lab_data <- lab_data |>
        dplyr::left_join(patient_data, by = "deId") |>
        impute_missforest() |>
        dplyr::select(
          tidyselect::all_of(colnames(lab_data))
        )
    }
  } else {
    # summarize preindex labs into a single baseline value via their mean
    lab_data <- lab_data |>
      dplyr::mutate(
        period = dplyr::case_when(
          LabDate <= 0 ~ 0,
          TRUE ~ LabDate
        )
      ) |>
      dplyr::group_by(LabName, deId, period) |>
      dplyr::summarise(LabValue = mean(LabValue), .groups = "drop")
    baseline_lab_data <- lab_data |>
      dplyr::filter(period == 0) |>
      tidyr::pivot_wider(
        id_cols = c(deId, period), 
        names_from = LabName,
        values_from = LabValue
      ) |>
      dplyr::right_join(
        tibble::tibble(deId = unique(patient_data$deId)), by = "deId"
      )
    # remove FeCal and crp due to high levels of missingness at baseline
    # nas_by_col <- apply(baseline_lab_data, 2, function(x) mean(is.na(x)))
    rm_na_cols <- c("FeCal", "crp", "esr")
    baseline_lab_data <- baseline_lab_data |>
      dplyr::select(-tidyselect::all_of(rm_na_cols))
    lab_data <- lab_data |>
      dplyr::filter(!(LabName %in% !!rm_na_cols))
    
    # impute baseline labs
    if (impute_mode == "median") {
      baseline_lab_data <- baseline_lab_data |>
        dplyr::mutate(
          dplyr::across(
            -deId, 
            ~ tidyr::replace_na(.x, median(.x, na.rm = TRUE))
          )
        )
    } else if (impute_mode == "mean") {
      baseline_lab_data <- baseline_lab_data |>
        dplyr::mutate(
          dplyr::across(
            -deId, 
            ~ tidyr::replace_na(.x, mean(.x, na.rm = TRUE))
          )
        )
    } else if (impute_mode == "rf") {
      baseline_lab_data <- baseline_lab_data |>
        dplyr::left_join(patient_data, by = "deId") |>
        dplyr::mutate(period = 0) |>
        impute_missforest(
          lab_vars = setdiff(colnames(baseline_lab_data), c("deId", "period"))
        ) |>
        dplyr::select(
          tidyselect::all_of(colnames(baseline_lab_data))
        )
    }
    
    baseline_lab_data <- baseline_lab_data |>
      tidyr::pivot_longer(
        cols = -c(deId, period),
        names_to = "LabName",
        values_to = "LabValue"
      )
    
    lab_data <- dplyr::bind_rows(
      lab_data |> dplyr::filter(period != 0),
      baseline_lab_data
    ) |>
      dplyr::arrange(deId, LabName, period)
  }

  return(lab_data)
}


#' Merge patient, flares, and lab data
#' 
#' @param patient_data A tibble containing the cleaned patient data.
#' @param flares_data A tibble containing the cleaned flares data.
#' @param lab_data A tibble containing the cleaned lab data. If NULL, only
#'   patient and flares data are merged.
#'   
#' @return A tibble containing the merged data.
merge_data <- function(patient_data, flares_data, lab_data = NULL) {
  merged_df <- dplyr::left_join(
    patient_data, flares_data, by = "deId"
  ) |>
    dplyr::mutate(
      event_date = ifelse(
        !is.na(event_date) & (event_date < days_to_death_or_end2021),
        event_date,
        NA
      ),
      switch_factor = factor(switch),
      dplyr::across(
        tidyselect::any_of(c("event_type", "disease_type")),
        ~ ifelse(is.na(event_date), "None", .x)
      ),
      dplyr::across(
        tidyselect::any_of("n_prior_flare"),
        ~ tidyr::replace_na(.x, 0)
      ),
      event = !is.na(event_date),
      censor_date = pmin(Date_Death, days_to_death_or_end2021, na.rm = TRUE),
      event_censor_date = pmin(censor_date, event_date, na.rm = TRUE),
      dplyr::across(where(is.integer), as.numeric)
    )
  
  if (!is.null(lab_data)) {
    merged_df <- dplyr::left_join(merged_df, lab_data, by = "deId")
  }
  
  return(merged_df)
}


#' Merge data for forecasting survival problem
#' 
#' @param patient_data A tibble containing the cleaned patient data.
#' @param flares_data A tibble containing the cleaned flares data.
#' @param lab_data A tibble containing the cleaned lab data.
#' @param forecast Length of the forecast period (in days).
#' 
#' @return A tibble containing the merged data.
merge_survival_data <- function(patient_data, 
                                flares_data, 
                                lab_data = NULL,
                                forecast = NULL) {
  
  merged_df <- merge_data(patient_data, flares_data)
  if (!is.null(forecast)) {
    merged_df <- merged_df |>
      dplyr::mutate(
        event_censor_date = event_censor_date - !!forecast
      ) |>
      dplyr::filter(event_censor_date > 0)
  }
  
  baseline_df <- merged_df |>
    dplyr::mutate(
      treatment = "originator",
      start_date = 0,
      end_date = ifelse(switch, newindexswitch_date, event_censor_date),
      event = dplyr::case_when(
        (switch == 1) & (event_date > newindexswitch_date) ~ 0,
        TRUE ~ event
      )
    )
  updated_df <- merged_df |>
    dplyr::filter(
      # only update switchers
      switch == 1, 
      # who didn't flare (or die) before switching
      event_censor_date > newindexswitch_date
    ) |>
    dplyr::mutate(
      treatment = "biosimilar",
      start_date = newindexswitch_date,
      end_date = event_censor_date
    )
  tmerged_df <- dplyr::bind_rows(baseline_df, updated_df) |>
    dplyr::arrange(deId, start_date) |>
    dplyr::mutate(
      treatment = relevel(as.factor(treatment), ref = "originator")
    )
  
  # merge in lab data
  if (!is.null(lab_data)) {
    lab_df <- lab_data |>
      dplyr::left_join(
        merged_df |> dplyr::select(deId, event_censor_date), 
        by = "deId"
      ) |>
      dplyr::filter(period < event_censor_date) |>
      tidyr::pivot_wider(
        id_cols = c(deId, period), 
        names_from = LabName, 
        values_from = LabValue
      ) |>
      dplyr::filter(deId %in% unique(tmerged_df$deId)) |>
      dplyr::arrange(deId, period)
    
    # merge labs with non-lab data
    tmerged_df <- dplyr::full_join(
      tmerged_df, lab_df, by = c("deId", "start_date" = "period")
    ) |>
      dplyr::arrange(deId, start_date) |>
      dplyr::select(-event) |>
      dplyr::group_by(deId) |>
      # fill in missing labs at time t with that of previous labs
      dplyr::mutate(
        dplyr::across(
          -c(start_date, end_date),
          ~ zoo::na.locf(.x, na.rm = FALSE)
        ),
        end_date = c(start_date[-1], max(end_date, na.rm = TRUE))
      ) |>
      dplyr::ungroup() |>
      dplyr::left_join(
        tmerged_df |>
          dplyr::filter(end_date == event_censor_date) |>
          dplyr::select(deId, end_date, event),
        by = c("deId", "end_date")
      ) |>
      dplyr::mutate(
        event = ifelse(is.na(event), 0, event)
      )
  }
  
  tmerged_df <- tmerged_df |>
    dplyr::rename(
      switch_orig = switch_factor,
      switch_factor = treatment
    )
  
  if (!isTRUE(all(tmerged_df$end_date > tmerged_df$start_date))) {
    stop("Error in merging time-varying covariates data. End dates are not all greater than start dates.")
  }
  
  return(tmerged_df)
}


#' Merge data for forecasting classification problem
#' 
#' @param patient_data A tibble containing the cleaned patient data.
#' @param flares_data A tibble containing the cleaned flares data.
#' @param lab_data A tibble containing the cleaned lab data.
#' @param accrual Length of the accrual period (in days).
#' @param forecast Length of the forecast period (in days).
#' 
#' @return A tibble containing the merged data.
merge_classification_data <- function(patient_data, 
                                      flares_data, 
                                      lab_data = NULL,
                                      accrual, 
                                      forecast) {
  merged_df <- dplyr::left_join(patient_data, flares_data, by = "deId") |>
    dplyr::select(-tidyselect::any_of("prior_event_dates")) |>
    dplyr::mutate(
      dplyr::across(
        c(n_prior_flare, event), 
        ~ tidyr::replace_na(.x, 0)
      ),
      censor_date = pmin(Date_Death, days_to_death_or_end2021, na.rm = TRUE)
    ) |>
    dplyr::filter(censor_date > (accrual + forecast)) |>
    dplyr::select(-newindexswitch_date)
  if (!is.null(lab_data)) {
    merged_df <- dplyr::left_join(merged_df, lab_data, by = "deId")
  }
  return(merged_df)
}
