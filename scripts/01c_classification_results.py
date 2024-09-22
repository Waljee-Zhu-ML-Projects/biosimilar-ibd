import pandas as pd
import numpy as np
import os
from os.path import join as oj
import pickle as pkl
import sys

from sklearn.metrics import roc_curve
from dcurves import dca

sys.path.append("..")
from functions.clean import clean_data


if __name__ == "__main__":

  root_dir = oj("..")
  data_dir = oj(root_dir, "data", "classification_data")
  out_dir = oj(root_dir, "results", "classification_models")

  seeds = np.arange(1, 51)
  forecasts = ["3", "6", "12"]
  accruals = ["1", "2"]
  censor_switchers = ["", "_censor_switcher"]
  impute_modes = ["rf", "median", "mean"]
  varset_modes = ["no_labs", "recent_labs", "recent_mean_max_labs", "recent_mean_max_recentolsslope_labs", "recent_accrual_labs"]

  # get prediction results
  for split_mode in ["validation", "test"]:
    test_errs_ls = []
    best_cutoff_ls = []
    calibration_ls = []
    dca_ls = []
    vimps_ls = []
    for censor_switcher in censor_switchers:
      for impute_mode in impute_modes:
        for forecast in forecasts:
          for accrual in accruals:
            for varset_mode in varset_modes:
              for seed in seeds:
                results_fname = oj(
                  out_dir,
                  f"classification_accrual_{accrual}yr_forecast_{forecast}mo_{impute_mode}_imputed{censor_switcher}_rep{seed}",
                  f"{split_mode}_{varset_mode}.pkl"
                )
                results = pkl.load(open(results_fname, "rb"))

                # test errors
                test_errs = results["test_errs"].reset_index(names="method")
                test_errs["censor_switcher"] = censor_switcher
                test_errs["impute_mode"] = impute_mode
                test_errs["forecast"] = forecast
                test_errs["accrual"] = accrual
                test_errs["varset_mode"] = varset_mode
                test_errs["seed"] = seed
                test_errs_ls.append(test_errs)

                # best cutoffs (sensitivity/specificity)
                test_preds_dict = results["test_prob_preds"]
                roc_out_ls = []
                for method, test_preds in test_preds_dict.items():
                  roc_out = roc_curve(results["test_y"], test_preds)
                  best_roc_idx = np.argmin((1 - roc_out[1])**2 + roc_out[0]**2)
                  roc_out_ls.append(
                    pd.DataFrame({
                      "seed": [seed],
                      "method": [method],
                      "sensitivity": [roc_out[1][best_roc_idx]],
                      "specificity": [1 - roc_out[0][best_roc_idx]],
                      "thr": [roc_out[2][best_roc_idx]]
                    })
                  )
                best_cutoff = pd.concat(roc_out_ls)
                best_cutoff["censor_switcher"] = censor_switcher
                best_cutoff["impute_mode"] = impute_mode
                best_cutoff["forecast"] = forecast
                best_cutoff["accrual"] = accrual
                best_cutoff["varset_mode"] = varset_mode
                best_cutoff["seed"] = seed
                best_cutoff_ls.append(best_cutoff)
                
                # calibration
                prob_preds = results["test_prob_preds"]["rfplus"]
                y_test = results["test_y"].ravel()
                q1 = np.quantile(prob_preds, 0.33)
                q2 = np.quantile(prob_preds, 0.66)
                q1_idx = np.where(prob_preds <= q1)[0]
                q2_idx = np.where((prob_preds > q1) & (prob_preds <= q2))[0]
                q3_idx = np.where(prob_preds > q2)[0]
                n_flares_q1 = np.sum(y_test[q1_idx])
                n_flares_q2 = np.sum(y_test[q2_idx])
                n_flares_q3 = np.sum(y_test[q3_idx])
                predicted_flares_q1 = len(q1_idx) * np.mean(prob_preds[q1_idx])
                predicted_flares_q2 = len(q2_idx) * np.mean(prob_preds[q2_idx])
                predicted_flares_q3 = len(q3_idx) * np.mean(prob_preds[q3_idx])
                calibration_df = pd.DataFrame({
                  "seed": [seed, seed, seed],
                  "method": ["rfplus", "rfplus", "rfplus"],
                  "censor_switcher": [censor_switcher, censor_switcher, censor_switcher],
                  "impute_mode": [impute_mode, impute_mode, impute_mode],
                  "forecast": [forecast, forecast, forecast],
                  "accrual": [accrual, accrual, accrual],
                  "varset_mode": [varset_mode, varset_mode, varset_mode],
                  "quantile": [1, 2, 3],
                  "num_flares": [n_flares_q1, n_flares_q2, n_flares_q3],
                  "predicted_flares": [predicted_flares_q1, predicted_flares_q2, predicted_flares_q3]
                })
                calibration_ls.append(calibration_df)
                
                for method in ["logistic", "rfplus"]:
                  dca_df = pd.DataFrame({
                    "y": results["test_y"],
                    "pred": results["test_prob_preds"][method]
                  })

                # decision curve analysis
                for method in ["logistic", "rfplus"]:
                  dca_df = pd.DataFrame({
                    "y": results["test_y"],
                    "pred": results["test_prob_preds"][method]
                  })
                  dca_out = dca(data=dca_df, outcome="y", modelnames=["pred"])
                  dca_out["method"] = method
                  dca_out["censor_switcher"] = censor_switcher
                  dca_out["impute_mode"] = impute_mode
                  dca_out["forecast"] = forecast
                  dca_out["accrual"] = accrual
                  dca_out["varset_mode"] = varset_mode
                  dca_out["seed"] = seed
                  dca_ls.append(dca_out)

                # feature importances
                vimps_dict = results["vimps"]
                for key in vimps_dict:
                  vimps_dict[key]["method"] = key
                  vimps_dict[key]["var"] = results["feature_names"]
                vimps = pd.concat(vimps_dict.values(), ignore_index=True)
                vimps["censor_switcher"] = censor_switcher
                vimps["impute_mode"] = impute_mode
                vimps["forecast"] = forecast
                vimps["accrual"] = accrual
                vimps["varset_mode"] = varset_mode
                vimps["seed"] = seed
                vimps_ls.append(vimps)

    test_errs = pd.concat(test_errs_ls)
    test_errs.to_csv(oj(out_dir, f"{split_mode}_errors.csv"))
    best_cutoff = pd.concat(best_cutoff_ls)
    best_cutoff.to_csv(oj(out_dir, f"{split_mode}_best_cutoff_results.csv"))
    calibration = pd.concat(calibration_ls)
    calibration.to_csv(oj(out_dir, f"{split_mode}_calibration_results.csv"))
    dca_out = pd.concat(dca_ls)
    dca_out.to_csv(oj(out_dir, f"{split_mode}_dca_results.csv"))
    vimps = pd.concat(vimps_ls)
    vimps.to_csv(oj(out_dir, f"{split_mode}_importances.csv"))

  print("Completed!")
  
