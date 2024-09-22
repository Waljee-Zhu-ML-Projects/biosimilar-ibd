rm(list = ls())

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
accrual_i <- as.integer(args[2])
forecast_i <- as.integer(args[3])
censor_switcher_i <- as.integer(args[4])

# initialize directory variables and load in helper functions
DATA_DIR <- here::here("data", "classification_data")
if (!dir.exists(DATA_DIR)) {
  dir.create(DATA_DIR, recursive = TRUE)
}
RESULTS_DIR <- here::here("results")
FUNCTIONS_DIR <- here::here("functions")
for (fname in list.files(FUNCTIONS_DIR, pattern = "\\.R$", full.names = TRUE)) {
  source(fname, chdir = TRUE)
}

# helper variables
yr <- 365.25
month <- 30
accruals <- c("1" = 1 * yr, "2" = 2 * yr)
forecasts <- c("3" = 3 * month, "6" = 0.5 * yr, "12" = 1 * yr)
censor_switchers <- c(FALSE, "_censor_switcher" = TRUE)

set.seed(1234 + seed)
censor_switcher <- censor_switchers[censor_switcher_i]
censor_switcher_name <- names(censor_switchers)[censor_switcher_i]
accrual <- accruals[accrual_i]
accrual_name <- names(accruals)[accrual_i]
forecast <- forecasts[forecast_i]
forecast_name <- names(forecasts)[forecast_i]

# remove labs (too much missing data or irrelevant)
rm_labs <- c("mch", "mchc", "mpv", "rbc", "hematocrit", "neutrophils", "FeCal", "crp")

# load in data
flares_df_orig <- load_flares_data()
patient_df_orig <- load_patient_data()
lab_df_orig <- load_lab_data()

# clean data
flares_df <- clean_flares_data(
  flares_df_orig, 
  accrual = accrual,
  forecast = forecast,
  mode = "classification"
)
patient_df <- clean_patient_data(
  patient_df_orig, 
  flares_df,
  accrual = accrual,
  forecast = forecast,
  censor_switcher = censor_switcher,
  mode = "classification"
)

# impute missing labs with median, mean, and RF imputation
for (impute_mode in c("median", "mean", "rf")) {
  lab_df <- clean_lab_data(
    lab_df_orig, 
    patient_df, 
    accrual = accrual,
    forecast = forecast,
    summarize = TRUE,
    index_mode = "event", 
    impute_mode = impute_mode
  )
  merged_df <- merge_classification_data(
    patient_df, flares_df, lab_df, accrual = accrual, forecast = forecast
  ) |>
    dplyr::select(-tidyselect::contains(rm_labs))
  write.csv(
    merged_df,
    file = file.path(
      DATA_DIR, 
      sprintf(
        "classification_accrual_%syr_forecast_%smo%s_%s_imputed_rep%s.csv",
        accrual_name, forecast_name, censor_switcher_name, impute_mode, seed
      )
    ),
    row.names = FALSE
  )
}

