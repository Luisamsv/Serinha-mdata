# Serinhaem monthly-discharge reproducibility package

This repository package contains the public monthly data and the complete R
workflow used for the manuscript.

## Files

- `T1_monthly_data.txt` — T1 temporal configuration:
  calibration 1966-01 to 2005-12; validation 2006-01 to 2019-04.
- `T2_monthly_data.txt` — T2 temporal configuration:
  calibration 1966-01 to 2011-12; validation 2012-01 to 2019-04.
- `data_dictionary.txt` — variable definitions, units, missing-value convention,
  and internal lag notation.
- `Serinhaem_complete_reproducible_analysis.R` — complete analysis workflow.

## Public data schema

The two tab-delimited data files use only descriptive English field names with
units:

- `date`
- `freshwater_discharge_m3_s`
- `precipitation_mm_month`
- `data_role`

`NA` denotes an unavailable original observation. Missing values are not
imputed.

## Running the analysis

Place the four files in the same directory and run:

```bash
Rscript Serinhaem_complete_reproducible_analysis.R T1_monthly_data.txt T2_monthly_data.txt
```

If no command-line arguments are supplied, the script searches for those two
filenames in the current working directory.

A timestamped `analysis_outputs_*` directory is created. It contains model
selection results, validation metrics, exact common-date T1/T2 comparisons,
predictive-error tests, sensitivity analyses, predictive intervals, figures,
data audits, package versions, the fixed random seed, and full `sessionInfo()`.

## Primary MLP implementation

The primary MLPs are single-hidden-layer feedforward networks fitted with
`nnet`, logistic hidden activation, linear output, BFGS optimization, squared
error, optional L2 weight decay, and a maximum of 1,000 iterations. No
validation-based early stopping, MLP frequency parameter, or MLP differencing
parameter is used. Three fixed-seed random starts are used per rolling-origin
fold and 20 final starts are averaged.

Primary held-out evaluation is sequential one-step-ahead prediction conditioned
on observed antecedent discharge. Precipitation-informed models use observed
historical precipitation at the required lags. Fully recursive discharge and
1-, 3-, 6-, 12-, and 24-month horizon analyses are reported separately as
sensitivity analyses.
