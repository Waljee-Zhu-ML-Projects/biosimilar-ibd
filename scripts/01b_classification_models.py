import argparse
import pandas as pd
import numpy as np
import os
from os.path import join as oj
import pickle as pkl
import sys
from functools import partial

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, average_precision_score
from sklearn.inspection import partial_dependence

from sklearn.linear_model import LogisticRegression
from imodels.importance import RandomForestPlusClassifier
from imodels.importance.ppms import LogisticClassifierPPM

sys.path.append("..")
from functions.fit import run_binary_classification_pipeline
from functions.clean import clean_data
from functions.importance import get_coefficient, get_feature_importances, get_mdiplus_importances


if __name__ == "__main__":

  default_root_dir = oj("..")
  default_data_dir = oj(default_root_dir, "data", "classification_data")
  default_out_dir = oj(default_root_dir, "results", "classification_models")

  parser = argparse.ArgumentParser()
  parser.add_argument('--seed', type=int, default=0)
  parser.add_argument('--data_dir', type=str, default=default_data_dir)
  parser.add_argument('--data_file', type=str, default="classification_accrual_1yr_forecast_3mo_mean_imputed_rep1.csv")
  parser.add_argument('--out_dir', type=str, default=default_out_dir)
  args = parser.parse_args()

  CV_PARAM_GRID_ALL = {
    "randf__min_samples_leaf": [1, 3, 5, 10],
    "logistic__C": [1],
    "lasso__C": np.logspace(-3, 3, 100),
    "ridge__C": np.logspace(-3, 3, 100),
    "elnet__C": np.logspace(-3, 3, 100),
    "elnet__l1_ratio": [0.1, 0.25, 0.5, 0.75, 0.9],
    "gb__learning_rate": [0.05, 0.1, 0.15],
    "gb__min_samples_leaf": [1, 5, 10],
    "gb__max_depth": [3, 5],
    "figs__max_rules": [5, 10, 12, 15, 20, 30, 50],
    "rfplus__prediction_model": [LogisticClassifierPPM()],
    "rulefit__max_rules": [5, 10, 30, 50]
  }

  MODELS = {
    "randf": RandomForestClassifier(n_estimators=500, min_samples_leaf=1),
    "gb": GradientBoostingClassifier(),
    "logistic": LogisticRegression(penalty=None),
    "lasso": LogisticRegression(penalty="l1", solver="liblinear"),
    "ridge": LogisticRegression(penalty="l2"),
    "elnet": LogisticRegression(penalty="elasticnet", solver="saga"),
    "rfplus": RandomForestPlusClassifier(
        rf_model=RandomForestClassifier(n_estimators=100, min_samples_leaf=1)
    ),
    "autogluon_medium": {"label": "event", "presets": "medium_quality"}
  }

  FI_MODELS = {
    "randf": get_feature_importances,
    "gb": get_feature_importances,
    "logistic": get_coefficient,
    "lasso": partial(get_coefficient, return_sparsity=True),
    "ridge": get_coefficient,
    "elnet": partial(get_coefficient, return_sparsity=True),
    "rfplus": get_mdiplus_importances,
    "autogluon_medium": None
  }

  METRICS = {
    "accuracy": (accuracy_score, "predict"),
    "f1": (f1_score, "predict"),
    "auroc": (roc_auc_score, "predict_proba"),
    "auprc": (average_precision_score, "predict_proba")
  }

  os.makedirs(oj(args.out_dir, os.path.splitext(args.data_file)[0]), exist_ok=True)

  data_orig = pd.read_csv(oj(args.data_dir, args.data_file))
  base_vars = [
    "switch_factor", "COMORBIDITY_cat", "Age_IFX", "Gender", "race_gp",
    "Ethnicity", "MS_cat", "priority", "concomitant_med", "TimeYrs_IBD_IFX",
    "IBD_type_c", "rurality", "count_anyenc", "count_IBDenc", "ProviderType_gp",
    "n_flares_preifx", "n_prior_flare"
  ]
  lab_vars = [v for v in data_orig.columns if ("lab" in v) and not ("n_labs" in v)]

  y = data_orig["event"]
  X_orig = data_orig[base_vars + lab_vars]

  var_sets = {
    "no_labs": base_vars,
    "recent_labs": base_vars + [v for v in lab_vars if "recent" in v and not ("recent_ols" in v)],
    "recent_mean_max_labs": base_vars + [v for v in lab_vars if ("recent" in v and not ("recent_ols" in v)) or "mean_lab_value_accrual" in v or "max_lab_value_accrual" in v],
    "recent_mean_max_recentolsslope_labs": base_vars + [v for v in lab_vars if ("recent" in v and not ("recent_ols" in v)) or "mean_lab_value_accrual" in v or "max_lab_value_accrual" in v or "recent_ols_slope_lab_value_accrual" in v],
    "recent_accrual_labs": base_vars + [v for v in lab_vars if ("recent" in v and not ("recent_ols" in v)) or "accrual" in v]
  }

  for split_mode in ["validation", "test"]:
    for var_set_name, var_set in var_sets.items():
      X = clean_data(X_orig[var_set])
      out_file = oj(args.out_dir, os.path.splitext(args.data_file)[0], f"{split_mode}_{var_set_name}.pkl")
      if os.path.exists(out_file):
        print(f"{split_mode}_{var_set_name} has already been run. Skipping...")
        continue
      X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, stratify=y, random_state=args.seed)
      X_valid, X_test, y_valid, y_test = train_test_split(X_test, y_test, test_size=0.5, stratify=y_test, random_state=100 + args.seed)
      if split_mode == "test":
        X_train = pd.concat((X_train, X_valid))
        y_train = pd.concat((y_train, y_valid))
      else:
        X_test = X_valid
        y_test = y_valid
      test_errs, test_preds, test_prob_preds, tuned_pipelines, vimps, boot_test_errs = run_binary_classification_pipeline(
        X_train=X_train, y_train=y_train, X_test=X_test, y_test=y_test,
        models_all=MODELS, fi_models_all=FI_MODELS, cv_param_grid_all=CV_PARAM_GRID_ALL,
        metrics=METRICS, importance=True, 
        path=oj(args.out_dir, os.path.splitext(args.data_file)[0], f"autogluon_{split_mode}_{var_set_name}")
      )
      pkl.dump(
        {
          "test_errs": test_errs,
          "test_preds": test_preds,
          "test_prob_preds": test_prob_preds,
          "test_y": y_test,
          "test_switch": X_test["switch_factor"],
          "tuned_pipelines": tuned_pipelines,
          "vimps": vimps,
          "feature_names": list(X_train.columns),
          "boot_test_errs": boot_test_errs
        },
        open(out_file, "wb")
      )

      if split_mode == "test" and var_set_name == "recent_mean_max_labs":
        # compute partial dependence plots
        clf = tuned_pipelines["rfplus"].best_estimator_[0]
        clf.classes_ = [False, True]
        features = [
          "n_prior_flare",
          "count_anyenc",
          "count_IBDenc",
          "max_lab_value_accrual_wbc",
          "recent_lab_value_hgb",
          "mean_lab_value_accrual_hgb",
          "max_lab_value_accrual_platelet",
          "n_flares_preifx",
          "mean_lab_value_accrual_co2"
        ]
        for feature in features:
          pdp_out = partial_dependence(
            clf, X_train, [feature], percentiles=(0, 1)
          )
          pdp_df = pd.DataFrame(
            {
              "grid_values": pdp_out["grid_values"][0],
              "average": pdp_out["average"].ravel()
            }
          )
          pdp_df.to_csv(
            oj(args.out_dir, os.path.splitext(args.data_file)[0], f"{split_mode}_{var_set_name}_pdp_{feature}.csv")
          )

  print("Completed!")
