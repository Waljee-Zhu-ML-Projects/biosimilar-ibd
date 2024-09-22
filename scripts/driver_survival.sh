# Step 1: Generate survival data
forecasts=(1 2 3)
censorswitchers=(1 2)

for forecast in "${forecasts[@]}"
do
  for censorswitcher in "${censorswitchers[@]}"
  do
    sbatch --job-name="survival_data_forecast${forecast}_censor${censorswitcher}" --mem-per-cpu=1GB --time=00:30:00 submit_r_script.sh 02a_survival_data ${forecast} ${censorswitcher}
  done
done

# Step 2: Fit survival models
forecasts=(1 2 3)
censorswitchers=(1 2)

for forecast in "${forecasts[@]}"
do
  for censorswitcher in "${censorswitchers[@]}"
  do
    sbatch --job-name="survival_models_forecast${forecast}_censor${censorswitcher}" submit_survival_models_job.sh ${forecast} ${censorswitcher}
  done
done


# Step 3: Aggregate all survival model results
sbatch --job-name="survival_results" submit_r_script.sh 02c_survival_results
