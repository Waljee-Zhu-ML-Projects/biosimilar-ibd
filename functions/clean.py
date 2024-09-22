import pandas as pd
from sklearn.preprocessing import StandardScaler, OneHotEncoder


def clean_data(X):
  # convert NaNs to "none" in COMORBIDITY_cat column in X
  X.loc[:, "COMORBIDITY_cat"] = X.loc[:, "COMORBIDITY_cat"].fillna("0")

  # one hot encode categorical variables and standardize X
  ohe_vars = [
    "COMORBIDITY_cat", "Gender", "race_gp", "Ethnicity", "MS_cat", "priority",
    "IBD_type_c", "rurality", "ProviderType_gp"
  ]
  ohe_vars = [v for v in ohe_vars if v in X.columns]
  ohe = OneHotEncoder(drop="if_binary", sparse_output=False)
  ohe.fit(X[ohe_vars])
  X_cat = ohe.transform(X[ohe_vars])
  X_cat = pd.DataFrame(X_cat, columns=ohe.get_feature_names_out(ohe_vars), index=X.index)
  X_ohe = pd.concat([X.drop(columns=ohe_vars), X_cat], axis=1)

  # scale data
  scaler = StandardScaler()
  X_ohe_scaled = scaler.fit_transform(X_ohe)
  X_ohe_scaled = pd.DataFrame(X_ohe_scaled, columns=scaler.get_feature_names_out(X_ohe.columns), index=X.index)

  return X_ohe_scaled