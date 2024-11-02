# Prediction of Infliximab Loss of Response or Adverse Event in a National Cohort of Patients with IBD: Comparison of Traditional and Machine Learning Models

This repository contains all code to reproduce the results and figures from:

Jason K. Hou, Tiffany M. Tang, Shubhada Sansgiry, Tony Van, Peter Richardson, Codey Pham, Francessca Cunningham, Jessica A. Baker, Ji Zhu, Akbar K. Waljee. "Prediction of Infliximab Loss of Response or Adverse Event in a National Cohort of Patients with IBD: Comparison of Traditional and Machine Learning Models." (2024).

In this work, we aimed to investigate the safety and effectiveness of non-medical biosimilar switch of Infliximab (IFX), a highly effective and cost effective treatment for Crohn’s disease (CD) and ulcerative colitis (UC). To that end, we developed various statistical and machine learning models to predict adverse clinical events among patients on maintenance IFX.

## Installing Dependencies

This project includes both and R and Python code. In general, R was used to perform the survival analysis while Python was used to perform the classification analysis. To facilitate reproducibility, we have created reproducible environments using [`renv`](https://rstudio.github.io/renv/articles/renv.html) and [`conda`](https://docs.conda.io/en/latest/) for R and Python, respectively. These reproducible environments, which include all necessary package dependencies, can be installed as follows:

In R/Rstudio, navigate to the root repository directory and run
```r
renv::restore()
```
to reproduce the R environment for this project.

In your terminal, navigate to the root repository directory and run
```bash
conda create --name <env> --file requirements.txt
```
to reproduce the Python environment for this project.

## Directory Structure

- **[scripts/](./scripts/)**: contains all scripts to reproduce survival and classification analyses
	- To run classification analysis:
		- [01a_classification_data.R](./scripts/01a_classification_data.R): loads and cleans data used to train classification models
		- [01b_classification_models.py](./scripts/01b_classification_models.py): fits classification models to predict IBD-related flares
		- [01c_classification_results.py](./scripts/01c_classification_results.py): aggregates prediction and interpretability results from classification models across multiple seeds, accruals, and forecasting periods
		- [driver_classification.sh](./scripts/driver_classification.sh): main driver script to run classification analysis across multiple data splits as well as different accrual and forecasting periods
	- To run survival analysis:
		- [02a_survival_data.R](./scripts/02a_survival_data.R): loads and cleans data used to train survival models
		- [02b_survival_models.R](./scripts/02b_survival_models.R): fits survival models to estimate risk of experiencing an IBD-related flare
		- [02c_survival_results.R](./scripts/02c_survival_results.R): aggregates prediction and interpretability results from survival models across multiple seeds and forecasting periods
		- [driver_survival.sh](./scripts/driver_survival.sh): main driver script to run survival analysis across multiple data splits as well as different forecasting periods
	- [03_paper_figures.R](./scripts/03_paper_figures.R): script to reproduce all figures presented in the manuscript
- **[functions/](./functions/)**: contains all functions necessary to reproduce survival and classification analyses
	- [load.R](./functions/load.R): functions to load in raw patient, flares, and laboratory data
	- [clean.R](./functions/clean.R): functions to clean and merge patient, flares, and laboratory data for survival and classification analyses
	- [clean.py](./functions/clean.py): functions to clean (mainly, one hot encode and standardize) data for classification analyses in Python
	- [fit.R](./functions/fit.R): functions to fit survival models in R
	- [fit.py](./functions/fit.py): functions to fit classification models in Python
	- [eval.R](./functions/eval.R): functions to evaluate the performance of survival models
	- [importance.py](./functions/importance.py): functions to extract and format feature importances from classifcation models
	- [utils-clean.R](./functions/utils-clean.R): utility functions used for data cleaning
	- [utils.R](./functions/utils.R): miscellaneous utility functions
- **[notebooks/](./notebooks/)**: contains files to reproduce [supplementary R Markdown documentation](https://waljee-zhu-ml-projects.github.io/biosimilar-ibd)
- **[renv/](./renv/)**: the R project library, containing all packages used by this project; automatically created by `renv`
- [.Rprofile](.Rprofile): the R project profile; automatically created by `renv`
- [renv.lock](renv.lock): the R project lockfile, containing metadata about every package to be re-installed on a new machine; automatically created by `renv`
- [requirements.txt](requirements.txt): the exported `conda` environment, needed to reproduce Python environment with all required dependencies
