# Step 1: Generate classification data
forecasts=(1 2 3)
accruals=(1 2)
censorswitchers=(1 2)

for forecast in "${forecasts[@]}"
do
  for accrual in "${accruals[@]}"
  do
    for censorswitcher in "${censorswitchers[@]}"
    do
      sbatch --job-name="classification_data_accrual${accrual}_forecast${forecast}_censor${censorswitcher}" submit_classification_data_job.sh  ${accrual} ${forecast} ${censorswitcher}
    done
  done
done

# Step 2: Fit classification models
forecasts=(3 6 12)
accruals=(1 2)
censorswitchers=("" "_censor_switcher")

for forecast in "${forecasts[@]}"
do
  for accrual in "${accruals[@]}"
  do
    for censorswitcher in "${censorswitchers[@]}"
    do
      sbatch --job-name="classification_model_accrual${accrual}_forecast${forecast}${censorswitcher}_mean" submit_classification_models_job.sh  "classification_accrual_${accrual}yr_forecast_${forecast}mo${censorswitcher}_mean_imputed"
      sbatch --job-name="classification_model_accrual${accrual}_forecast${forecast}${censorswitcher}_median" submit_classification_models_job.sh  "classification_accrual_${accrual}yr_forecast_${forecast}mo${censorswitcher}_median_imputed"
      sbatch --job-name="classification_model_accrual${accrual}_forecast${forecast}${censorswitcher}_rf" submit_classification_models_job.sh  "classification_accrual_${accrual}yr_forecast_${forecast}mo${censorswitcher}_rf_imputed"
    done
  done
done

# Step 3: Aggregate all classification model results
sbatch --job-name="classification_results" submit_python_script.sh 01c_classification_results
