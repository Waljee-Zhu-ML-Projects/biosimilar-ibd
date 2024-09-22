import pandas as pd
import numpy as np
import copy

from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, roc_auc_score, average_precision_score
from autogluon.tabular import TabularPredictor


def run_binary_classification_pipeline(
  X_train, y_train, X_test, y_test,
  models_all, cv_param_grid_all, fi_models_all,
  metrics=None, keep_models=None, importance=False, B=100,
  path=None
):
  """
  Fit and evaluate binary classification models
  
  Parameters
  ----------
  X_train : pd.DataFrame
    Training data
  y_train : pd.Series
    Training labels
  X_test : pd.DataFrame
    Test data
  y_test : pd.Series
    Test labels
  models_all : dict
    Dictionary of models to fit
  cv_param_grid_all : dict
    Dictionary of CV parameters for each model
  fi_models_all : dict
    Dictionary of feature importance functions for each model
  metrics : dict
    Dictionary of metrics to evaluate
  keep_models : list
    List of models to fit
  importance : bool
    Whether to compute feature importances
  B : int
    Number of bootstrap samples
  path : str
    Scratch path to save autogluon output
    
  Returns
  -------
  test_errs : pd.DataFrame
    Test errors
  test_preds : dict
    Test predictions
  test_prob_preds : dict
    Test probability predictions
  tuned_pipelines : dict
    Tuned model pipelines
  vimps : dict
    Feature importances
  boot_test_errs : dict
    Bootstrap test errors
  """
  if metrics is None:
    metrics = {
      "accuracy": (accuracy_score, "predict"),
      "auroc": (roc_auc_score, "predict_proba"),
      "auprc": (average_precision_score, "predict_proba")
    }
  if keep_models is not None:
    models = {}
    for keep_model in keep_models:
        models[keep_model] = copy.deepcopy(models_all[keep_model])
  else:
    models = models_all
  pipes = {}
  for model_name, model in models.items():
    if "autogluon" in model_name:
      pipe = model
    else:
      pipe = Pipeline(steps=[(model_name, model)])
    pipes[model_name] = pipe

  test_errs = {}
  boot_test_errs = {}
  test_preds = {}
  test_prob_preds = {}
  tuned_pipelines = {}
  vimps = {}
  for pipe_name, pipe in pipes.items():
    print(pipe_name)
    # run CV for pipeline
    if "autogluon" in pipe_name:
      train_df = pd.concat([X_train, y_train], axis=1)
      pipe_search = TabularPredictor(label=pipe["label"], problem_type="binary", path=path)
    else:
      # get relevant CV parameters given the steps of the pipeline
      cv_param_grid = {
        key: cv_param_grid_all[key] for key in cv_param_grid_all.keys() \
          if key.startswith(tuple(pipe.named_steps.keys()))
      }
      if len(cv_param_grid) == 0:
        pipe_search = pipe
      else:
        pipe_search = GridSearchCV(pipe, cv_param_grid)
    if "rfplus" in pipe_name:
      pipe_search.fit(X_train.to_numpy(), y_train.to_numpy().ravel())
    elif "autogluon" in pipe_name:
      pipe_search.fit(train_df, presets=pipe["presets"])
    else:
      pipe_search.fit(X_train, y_train.to_numpy().ravel())
    tuned_pipelines[pipe_name] = copy.deepcopy(pipe_search)

    # make predictions
    preds_out = pipe_search.predict(X_test)
    if isinstance(preds_out, pd.Series):
      preds_out = preds_out.to_numpy()
    if len(np.unique(preds_out)) > 2:
      test_prob_preds[pipe_name] = preds_out
      test_preds[pipe_name] = (preds_out > 0.5).astype(int)
    else:
      test_preds[pipe_name] = preds_out
      if hasattr(pipe_search, "predict_proba"):
        preds_proba_out = pipe_search.predict_proba(X_test)
        if isinstance(preds_proba_out, pd.DataFrame):
          preds_proba_out = preds_proba_out.to_numpy()
        if preds_proba_out.ndim == 2:
          test_prob_preds[pipe_name] = preds_proba_out[:, 1]
        else:
          test_prob_preds[pipe_name] = preds_proba_out

    # evaluate predictions
    err_out = {}
    for metric_name, metric in metrics.items():
      metric_fun = metric[0]
      metric_type = metric[1]
      if metric_type == "predict":
        err_out[metric_name] = metric_fun(y_test, test_preds[pipe_name])
      elif metric_type == "predict_proba":
        if pipe_name in test_prob_preds.keys():
          err_out[metric_name] = metric_fun(y_test, test_prob_preds[pipe_name])
        else:
          err_out[metric_name] = np.nan
    test_errs[pipe_name] = copy.deepcopy(err_out)

    # evaluate predictions via bootstrap to get uncertainty estimate
    boot_test_errs[pipe_name] = {}
    for metric_name, metric in metrics.items():
      metric_fun = metric[0]
      metric_type = metric[1]
      if metric_type == "predict":
        y_est = test_preds[pipe_name]
      elif metric_type == "predict_proba":
        if pipe_name in test_prob_preds.keys():
          y_est = test_prob_preds[pipe_name]
        else:
          continue
      boot_test_errs_ls = []
      for b in range(B):
        boot_idx = np.random.choice(len(y_test), len(y_test))
        boot_test_errs_ls.append(metric_fun(y_test.iloc[boot_idx], y_est[boot_idx]))
      boot_test_errs[pipe_name][metric_name] = copy.deepcopy(boot_test_errs_ls)

    # evaluate feature importances
    if importance:
      fi_model = fi_models_all[pipe_name]
      if fi_model is not None:
        if hasattr(pipe_search, "best_estimator_"):
          best_estimator = pipe_search.best_estimator_[0]
        else:
          best_estimator = pipe_search
        vimp_df = fi_model(
          best_estimator,
          X_train=X_train, y_train=y_train, X_test=X_test, y_test=y_test
        )
        vimp_df = vimp_df.rename(columns={vimp_df.columns[0]: "var"})
        vimps[pipe_name] = copy.deepcopy(vimp_df)

    # if "rfplus" in pipe_name:
    #   tuned_pipelines[pipe_name] = None

  test_errs = pd.DataFrame(test_errs).T

  return test_errs, test_preds, test_prob_preds, tuned_pipelines, vimps, boot_test_errs
