rm(list = ls())

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
forecast_i <- as.integer(args[1])
censor_switcher_i <- as.integer(args[2])

# initialize directory variables and load in helper functions
DATA_DIR <- here::here("data", "survival_data")
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
forecasts <- c("3" = 3 * month, "6" = 0.5 * yr, "12" = 1 * yr)
censor_switchers <- c(FALSE, "_censor_switcher" = TRUE)

set.seed(1234)
censor_switcher <- censor_switchers[censor_switcher_i]
censor_switcher_name <- names(censor_switchers)[censor_switcher_i]
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
  forecast = forecast,
  mode = "survival"
)
patient_df <- clean_patient_data(
  patient_df_orig, 
  flares_df,
  forecast = forecast,
  censor_switcher = censor_switcher,
  mode = "survival"
)

# impute missing labs with median, mean, and RF imputation
lab_df_ls <- list()
for (impute_mode in c("median", "mean", "rf")) {
  lab_df_ls[[impute_mode]] <- clean_lab_data(
    lab_df_orig, 
    patient_df,
    summarize = FALSE,
    impute_mode = impute_mode
  )
}

# merge all data together
survival_df <- merge_survival_data(
  patient_df, 
  flares_df,
  forecast = forecast
)
survival_lab_df_median_imp <- merge_survival_data(
  patient_df, 
  flares_df,
  lab_df_ls[["median"]],
  forecast = forecast
)
survival_lab_df_mean_imp <- merge_survival_data(
  patient_df, 
  flares_df,
  lab_df_ls[["mean"]],
  forecast = forecast
)
survival_lab_df_rf_imp <- merge_survival_data(
  patient_df, 
  flares_df,
  lab_df_ls[["rf"]],
  forecast = forecast
)

save(
  survival_df, 
  survival_lab_df_median_imp, 
  survival_lab_df_mean_imp, 
  survival_lab_df_rf_imp,
  file = file.path(
    DATA_DIR,
    sprintf("survival_forecast_%smo%s.RData", forecast_name, censor_switcher_name)
  )
)
