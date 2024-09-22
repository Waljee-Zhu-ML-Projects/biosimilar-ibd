#' Load in patient data
#' 
#' @param data_dir The directory containing the data files
#' 
#' @return A tibble containing the full patient data
load_patient_data <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- file.path(Sys.getenv("Biosimilar"))
  }
  data_dir1 <- file.path(data_dir, "Biosimilar Data Original Files")
  data_dir2 <- file.path(data_dir, "Biosimilar_DeIdentified_11Jan2024")
  data_dir3 <- file.path(data_dir, "Biosimilar_DeIdentified_12March2024")
  data_dir4 <- file.path(data_dir, "Biosimilar_DeIdentified_12April2024_TT")

  demo <- data.table::fread(
    file.path(data_dir1, "Patient_Char", "PtDemog.csv")
	) |>
    tibble::as_tibble()
  ibdtype <- data.table::fread(
    file.path(data_dir1, "IBDtype", "ibdtype_alldef022022.csv")
  ) |>
    tibble::as_tibble()
  idx_switch <- data.table::fread(
    file.path(data_dir1, "Outcomes", "Index_Switch", "ifx_ml_newswitchdef2023.csv")
  ) |>
    tibble::as_tibble()
  steroids <- data.table::fread(
    file.path(data_dir3, "Concomitent_IBDmeds", "steroid032024.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::group_by(deId) |>
    dplyr::summarise(
      steroids = 1,
      .groups = "keep"
    ) |> 
    dplyr::ungroup()
  steroids2 <- data.table::fread(
    file.path(data_dir3, "Concomitent_IBDmeds", "steroid032024.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(flag_dx_withSteroid == 0) |>
    dplyr::group_by(deId) |>
    dplyr::summarise(
      steroids_no_dx = 1,
      .groups = "keep"
    ) |> 
    dplyr::ungroup()
  meds <- data.table::fread(
    file.path(data_dir3, "Concomitent_IBDmeds", "concmeds032024.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::group_by(deId) |>
    dplyr::summarise(
      concomitant_med = 1,
      .groups = "keep"
    ) |>
    dplyr::ungroup()
  priority <- data.table::fread(
    file.path(data_dir2, "Patient_Charc", "ml_priority.csv")
  ) |>
    tibble::as_tibble()
  rurality <- data.table::fread(
    file.path(data_dir2, "Patient_Charc", "ml_rurality.csv")
  ) |>
    tibble::as_tibble()
  enc <- data.table::fread(
    file.path(data_dir3, "Patient_Charc", "enc_InOut_12mpre_ifxindex032024.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::select(
      deId, 
      count_IBDenc = count_IBDencounter_12mPreIFX,
      count_anyenc = count_encounter_12mPreIFX
    )
  providertype <- data.table::fread(
    file.path(data_dir2, "Provider_Charc", "ml_providertype.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::select(-first_IFXdate_post17)
  censortimes <- data.table::fread(
    file.path(data_dir2, "Censoring", "pt_censoring.csv")
  ) |>
    tibble::as_tibble()
  
  preifx_allergy <- data.table::fread(
    file.path(data_dir4, "Pre1yr_IFX", "ml_allergy_1yrpreifx.csv")
  ) |>
    tibble::as_tibble()
  preifx_er <- data.table::fread(
    file.path(data_dir4, "Pre1yr_IFX", "ml_outcome_ervisits_1yrpreifx.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::select(deId)
  preifx_hospital <- data.table::fread(
    file.path(data_dir4, "Pre1yr_IFX", "ml_outcome_hospital_1yrpreifx.csv")
  ) |>
    tibble::as_tibble() |>
    dplyr::select(deId)
  preifx_flares <- dplyr::bind_rows(preifx_er, preifx_hospital) |>
    dplyr::group_by(deId) |>
    dplyr::summarize(
      n_flares_preifx = dplyr::n()
    )
  preifx_steroids <- data.table::fread(
    file.path(data_dir4, "Pre1yr_IFX", "steroid_1yrpreifxindex.csv")
  ) |>
    tibble::as_tibble()
  preifx_anaphylaxix <- data.table::fread(
    file.path(data_dir4, "Pre1yr_IFX", "ml_anaphylaxix_1yrpreifx.csv")
  ) |>
    tibble::as_tibble()

  patient_df <- idx_switch |>
    dplyr::left_join(preifx_flares, by = "deId") |>
    dplyr::left_join(ibdtype, by = "deId") |>
    dplyr::left_join(demo, by = "deId") |>
    dplyr::left_join(steroids, by = "deId") |>
    dplyr::left_join(steroids2, by = "deId") |>
    dplyr::left_join(meds, by = "deId") |>
    dplyr::left_join(priority, by = "deId") |>
    dplyr::left_join(rurality, by = "deId") |>
    dplyr::left_join(enc, by = "deId") |>
    dplyr::left_join(providertype, by = "deId") |>
    dplyr::left_join(censortimes, by = "deId") |>
    dplyr::mutate(
      n_flares_preifx = tidyr::replace_na(n_flares_preifx, 0),
      Age_IFX = (first_IFXdate_post17 - Date_Birth) / 365.25,
      yrs_ifx_biosim = (newindexswitch_date - first_IFXdate_post17) / 365.25,
      ifx_yr = (first_IFXdate_post17 - min(first_IFXdate_post17)) / 365.25,
      Gender = relevel(as.factor(Gender), ref = "M"),
      IBDtype_max = as.factor(IBDtype_max),
      concomitant_med = ifelse(is.na(concomitant_med), 0, concomitant_med),
      steroids = ifelse(is.na(steroids), 0, steroids),
      steroids_no_dx = ifelse(is.na(steroids_no_dx), 0, steroids_no_dx),
      count_anyenc = ifelse(is.na(count_anyenc), 0, count_anyenc),
      count_IBDenc = ifelse(is.na(count_IBDenc), 0, count_IBDenc),
      yrs_to_death_or_end2021 = days_to_death_or_end2021 / 365.25,
      yrs_to_death_or_end2022 = days_to_death_or_end2022 / 365.25
    )
	return(patient_df)
}


#' Load in flares data
#' 
#' @param data_dir The directory containing the data files
#' 
#' @return A tibble containing the flares data
load_flares_data <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    base_data_dir <- Sys.getenv("Biosimilar")
    data_dir <- file.path(base_data_dir, "Biosimilar Data Original Files")
  } else {
    base_data_dir <- dirname(data_dir)
  }

  hosp <- data.table::fread(
    file.path(data_dir, "Outcomes", "Hospital", "ml_Outcome_hospital.csv")
  ) |>
    tibble::as_tibble() |> 
    dplyr::rename(event_date = AdmitDate) |> 
    dplyr::mutate(event_type = "hosp")
  er <- data.table::fread(
    file.path(data_dir, "Outcomes", "ER_visits", "ml_outcome_ervisits.csv")
  ) |>
    tibble::as_tibble() |> 
    dplyr::rename(event_date = visitdate) |> 
    dplyr::mutate(event_type = "er")
  censor_df <- load_patient_data(base_data_dir) |>
    dplyr::select(deId, days_to_death_or_end2021)

  flares_all_df <- dplyr::bind_rows(er, hosp) |>
    dplyr::mutate(
      disease_type = dplyr::case_when(
        Outpat_ICD10Code_CRind == "CR" ~ "CR",
        TRUE ~ "UC"
      )
    ) |>
    dplyr::select(-Outpat_ICD10Code_CRind, -Outpat_ICD10Code_UCind) |>
    dplyr::arrange(deId, event_date) |>
    dplyr::left_join(censor_df, by = "deId") |>
    # filter out events that occur after the end of the study
    dplyr::filter(event_date <= days_to_death_or_end2021) |>
    dplyr::group_by(deId) |>
    dplyr::mutate(
      # create a flag for the first flare
      first = c(TRUE, rep(FALSE, dplyr::n() - 1)),
      # number of total flares
      n = dplyr::n(),
      n_cr = sum(disease_type == "CR"),
      n_uc = sum(disease_type == "UC"),
      n_cr_and_uc = any(disease_type == "CR") && any(disease_type == "UC")
    ) |>
    dplyr::ungroup()
  
	return(flares_all_df)
}


#' Load in lab data
#' 
#' @param data_dir The directory containing the data files
#' 
#' @return A tibble containing the lab data in long format
load_lab_data <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- file.path(Sys.getenv("Biosimilar"), "Biosimilar Data Original Files")
  }
	lab_files <- list.files(file.path(data_dir, "Outcomes", "labs"), pattern = ".csv")
  lab_names <- stringr::str_remove(lab_files, "\\.csv$")
  names(lab_names) <- lab_names
  lab_data_ls <- purrr::map(
    lab_names, 
    ~ data.table::fread(
      file.path(data_dir, "Outcomes", "labs", sprintf("%s.csv", .x))
    ) |>
      tibble::as_tibble()
  )
  lab_data <- lab_data_ls |>
    dplyr::bind_rows(.id = "LabName")
  return(lab_data)
}
