# ============================================================================== 
# COMPLETE REPRODUCIBLE ANALYSIS — MONTHLY DISCHARGE OF THE SERINHAEM RIVER
# ============================================================================== 
# This script reproduces the statistical and modeling analyses used in the
# revised manuscript. Output files use descriptive names and are intentionally
# independent of manuscript table or figure numbering.
#
# Expected inputs: T1_monthly_data.txt and T2_monthly_data.txt.
# Both are tab-delimited files with the strict public schema:
#   date, freshwater_discharge_m3_s, precipitation_mm_month, data_role.
# Units are embedded in the field names. Missing observations are encoded as NA
# and are retained at their original monthly positions; no values are imputed.
#
# Primary workflow:
#   - temporal partitions T1 and T2;
#   - rolling-origin model selection using calibration data only;
#   - linear regression, ARMA, ARIMA, SARIMA, seasonal-naive and MLP models;
#   - one-step, recursive and multi-horizon evaluation;
#   - absolute, normalized, bias, skill and variance metrics;
#   - predictor-matched comparisons and paired predictive-error inference;
#   - post-2012 common-period comparisons;
#   - hydroclimatic change, stationarity, trend and autocorrelation analyses;
#   - flow-dependent error analyses;
#   - MLP lag-window, activation-function and learning-rate sensitivity checks;
#   - moving-block residual bootstrap predictive intervals.
# ============================================================================== 


# ------------------------------------------------------------------------------
# Utility: moving-block resampling for serially dependent paired quantities
# ------------------------------------------------------------------------------
sample_moving_blocks <- function(x, n = length(x), block_length = 12L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < block_length) stop("Insufficient observations for moving-block resampling.")
  starts <- seq_len(length(x) - block_length + 1L)
  n_blocks <- ceiling(n / block_length)
  chosen <- sample(starts, n_blocks, replace = TRUE)
  out <- unlist(lapply(chosen, function(i) x[i:(i + block_length - 1L)]), use.names = FALSE)
  out[seq_len(n)]
}

moving_block_mean_ci <- function(x, B = 2000L, block_length = 12L, seed = 123L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(lower = NA_real_, upper = NA_real_))
  set.seed(seed)
  sims <- replicate(B, mean(sample_moving_blocks(x, length(x), block_length)))
  c(lower = unname(quantile(sims, 0.025, na.rm = TRUE)),
    upper = unname(quantile(sims, 0.975, na.rm = TRUE)))
}


# ------------------------------------------------------------------------------
# 0. Configuration and package setup
# ------------------------------------------------------------------------------
PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)

T1_DATA_FILE <- if (length(args) >= 1L) {
  args[[1]]
} else {
  file.path(PROJECT_DIR, "T1_monthly_data.txt")
}

T2_DATA_FILE <- if (length(args) >= 2L) {
  args[[2]]
} else {
  file.path(PROJECT_DIR, "T2_monthly_data.txt")
}

RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(PROJECT_DIR, paste0("analysis_outputs_", RUN_ID))
SEED <- 123L
MLP_FINAL_RESTARTS <- 20L
MLP_CV_RESTARTS <- 3L

T1_CAL_START <- as.Date("1966-01-01")
T1_CAL_END    <- as.Date("2005-12-01")
T1_VAL_START <- as.Date("2006-01-01")
T2_CAL_END    <- as.Date("2011-12-01")
T2_VAL_START <- as.Date("2012-01-01")
END_DATE    <- as.Date("2019-04-01")

CV_INITIAL_YEARS <- 20L
CV_LAST_MONTHS <- 120L
CV_STEP <- 3L

HORIZONS <- c(1L, 3L, 6L, 12L, 24L)
HORIZON_STEP <- 6L

# Moving-block bootstrap settings used for uncertainty and paired inference.
BOOT_B <- 2000L
BOOT_BLOCK_LENGTH <- 12L

PRIMARY_EVALUATION <- "one_step"

packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "lubridate", "tibble", "ggplot2",
  "forecast", "nnet", "trend", "tseries"
)

missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the missing packages before running:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lubridate)
  library(tibble)
  library(ggplot2)
  library(forecast)
  library(nnet)
  library(trend)
  library(tseries)
})

set.seed(SEED)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 1. Public-data reading, configuration audit and missing-data audit
# ------------------------------------------------------------------------------
PUBLIC_FIELDS <- c(
  "date",
  "freshwater_discharge_m3_s",
  "precipitation_mm_month",
  "data_role"
)

read_configuration_data <- function(path, configuration) {
  if (!file.exists(path)) stop("File not found: ", path)

  x <- readr::read_tsv(
    path,
    na = c("NA", ""),
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  missing_fields <- setdiff(PUBLIC_FIELDS, names(x))
  if (length(missing_fields) > 0L) {
    stop(
      "Missing required public-data field(s) in ", basename(path), ": ",
      paste(missing_fields, collapse = ", ")
    )
  }

  x <- x %>%
    select(all_of(PUBLIC_FIELDS)) %>%
    mutate(
      date = as.Date(date),
      freshwater_discharge_m3_s = as.numeric(freshwater_discharge_m3_s),
      precipitation_mm_month = as.numeric(precipitation_mm_month),
      data_role = tolower(trimws(as.character(data_role)))
    ) %>%
    arrange(date)

  if (anyNA(x$date)) stop("Invalid date values were found in ", basename(path), ".")
  if (anyDuplicated(x$date)) stop("Duplicate monthly dates were found in ", basename(path), ".")

  expected_dates <- seq(min(x$date), max(x$date), by = "month")
  if (!identical(x$date, expected_dates)) {
    stop("The monthly time grid is incomplete or out of order in ", basename(path), ".")
  }

  expected_role <- if (configuration == "T1") {
    ifelse(x$date <= T1_CAL_END, "calibration", "validation")
  } else if (configuration == "T2") {
    ifelse(x$date <= T2_CAL_END, "calibration", "validation")
  } else {
    stop("Unknown configuration: ", configuration)
  }

  if (!identical(x$data_role, expected_role)) {
    stop(
      "The data_role field does not match the declared ", configuration,
      " calibration/validation boundaries in ", basename(path), "."
    )
  }

  x
}

t1_public <- read_configuration_data(T1_DATA_FILE, "T1")
t2_public <- read_configuration_data(T2_DATA_FILE, "T2")

# The T1 and T2 files intentionally contain the same physical monthly observations.
# They differ only in the calibration/validation role assigned to each date.
if (!identical(t1_public$date, t2_public$date)) {
  stop("T1 and T2 public files do not contain the same monthly dates.")
}

same_numeric_with_na <- function(x, y, tolerance = 1e-12) {
  same_na <- identical(is.na(x), is.na(y))
  if (!same_na) return(FALSE)
  ok <- !is.na(x)
  if (!any(ok)) return(TRUE)
  all(abs(x[ok] - y[ok]) <= tolerance)
}

if (!same_numeric_with_na(
  t1_public$freshwater_discharge_m3_s,
  t2_public$freshwater_discharge_m3_s
)) {
  stop("T1 and T2 files contain different freshwater-discharge observations.")
}

if (!same_numeric_with_na(
  t1_public$precipitation_mm_month,
  t2_public$precipitation_mm_month
)) {
  stop("T1 and T2 files contain different precipitation observations.")
}

monthly_data <- t1_public %>%
  select(date, freshwater_discharge_m3_s, precipitation_mm_month) %>%
  filter(date >= T1_CAL_START, date <= END_DATE) %>%
  mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    freshwater_discharge_observed = is.finite(freshwater_discharge_m3_s),
    precipitation_observed = is.finite(precipitation_mm_month),
    freshwater_discharge_original_m3_s = freshwater_discharge_m3_s,
    precipitation_original_mm_month = precipitation_mm_month
  )

# Audit the exact observation counts used in the revised manuscript.
data_audit <- tibble(
  item = c(
    "calendar_months",
    "freshwater_discharge_observed",
    "precipitation_observed",
    "both_variables_observed"
  ),
  value = c(
    nrow(monthly_data),
    sum(monthly_data$freshwater_discharge_observed),
    sum(monthly_data$precipitation_observed),
    sum(monthly_data$freshwater_discharge_observed &
          monthly_data$precipitation_observed)
  )
)
write_csv(data_audit, file.path(OUTPUT_DIR, "public_data_audit.csv"))

expected_audit <- c(640L, 634L, 631L, 628L)
if (!identical(as.integer(data_audit$value), expected_audit)) {
  stop(
    "Public-data audit does not match the revised analysis. Expected: ",
    paste(expected_audit, collapse = ", "),
    "; obtained: ", paste(data_audit$value, collapse = ", "), "."
  )
}

missing_rows <- monthly_data %>%
  filter(!freshwater_discharge_observed | !precipitation_observed)
write_csv(missing_rows, file.path(OUTPUT_DIR, "missing_data_audit.csv"))

period_audit <- tibble(
  Period_label = c(
    "Full series",
    "T1 calibration",
    "T1 validation",
    "T2 calibration",
    "T2 validation",
    "Common post-2012 period"
  ),
  Start = c(
    min(monthly_data$date),
    T1_CAL_START,
    T1_VAL_START,
    T1_CAL_START,
    T2_VAL_START,
    T2_VAL_START
  ),
  End = c(
    END_DATE,
    T1_CAL_END,
    END_DATE,
    T2_CAL_END,
    END_DATE,
    END_DATE
  )
) %>%
  rowwise() %>%
  mutate(
    N_months = sum(monthly_data$date >= Start & monthly_data$date <= End),
    N_discharge_observed = sum(
      monthly_data$date >= Start & monthly_data$date <= End &
        monthly_data$freshwater_discharge_observed
    ),
    N_precipitation_observed = sum(
      monthly_data$date >= Start & monthly_data$date <= End &
        monthly_data$precipitation_observed
    )
  ) %>%
  ungroup()
write_csv(period_audit, file.path(OUTPUT_DIR, "period_audit.csv"))

# Repository-facing data dictionary. Qk and Pk below are internal model notation
# for k-month lags and are not public raw-data column names.
data_dictionary <- tribble(
  ~field, ~definition, ~unit_or_values,
  "date",
  "First day of the month representing the monthly observation.",
  "YYYY-MM-DD",
  "freshwater_discharge_m3_s",
  "Observed monthly river discharge at the Itubera gauging station.",
  "m^3 s^-1",
  "precipitation_mm_month",
  "Observed monthly precipitation at the Camamu rain gauge.",
  "mm month^-1",
  "data_role",
  "Role of the observation in the named temporal configuration.",
  "calibration | validation",
  "NA",
  "Missing-value code. Missing observations are not imputed.",
  "not applicable",
  "Qk",
  "Internal model notation for freshwater discharge lagged by k months.",
  "m^3 s^-1",
  "Pk",
  "Internal model notation for precipitation lagged by k months.",
  "mm month^-1"
)
write_tsv(data_dictionary, file.path(OUTPUT_DIR, "data_dictionary.tsv"))

# ------------------------------------------------------------------------------
# 2. Lagged predictors, performance metrics and helper functions
# ------------------------------------------------------------------------------
# Each modeling row targets freshwater discharge at month t. Qk and Pk denote
# the observed freshwater-discharge and precipitation values at t-k months.
# Missing targets or required lagged predictors remain missing and are excluded
# only from the model configuration that requires them. No imputation is used.
#
# A feedforward MLP has no intrinsic recurrent or seasonal state. Temporal and
# annual information is therefore supplied explicitly through lagged inputs,
# including Q12/P12 where applicable. Continuous 12- and 24-month lag windows
# are evaluated later as sensitivity analyses.
max_lag <- 24L
lagged_data <- monthly_data
for (k in seq_len(max_lag)) {
  lagged_data[[paste0("Q", k)]] <- dplyr::lag(lagged_data$freshwater_discharge_m3_s, k)
  lagged_data[[paste0("P", k)]] <- dplyr::lag(lagged_data$precipitation_mm_month, k)
}

univariate_sets <- list(
  U1 = c("Q1"),
  U2 = c("Q1", "Q2"),
  U3 = c("Q1", "Q2", "Q3"),
  U4 = c("Q1", "Q2", "Q12"),
  U5 = c("Q1", "Q2", "Q3", "Q12")
)

multivariate_sets <- list(
  M1 = c("Q1", "P1"),
  M2 = c("Q1", "Q2", "P1"),
  M3 = c("Q1", "Q2", "P1", "P2"),
  M4 = c("Q1", "Q2", "Q12", "P1", "P2"),
  M5 = c("Q1", "Q2", "Q12", "P1", "P12"),
  M6 = c("Q1", "Q2", "Q3", "P1", "P2"),
  M7 = c("Q1", "Q2", "Q3", "Q12", "P1", "P2"),
  M8 = c("Q1", "Q2", "Q12", "P1", "P2", "P12")
)

mase_scale <- function(obs_cal, seasonal_lag = 12L) {
  current <- obs_cal[(seasonal_lag + 1L):length(obs_cal)]
  previous <- obs_cal[1L:(length(obs_cal) - seasonal_lag)]
  ok <- is.finite(current) & is.finite(previous)
  if (!any(ok)) return(NA_real_)
  mean(abs(current[ok] - previous[ok]))
}

compute_metrics <- function(obs, pred, obs_cal) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  if (length(obs) < 2L) return(tibble(
    N = length(obs), RMSE = NA_real_, MAE = NA_real_, NSE = NA_real_,
    KGE2012 = NA_real_, PBIAS = NA_real_, RSR = NA_real_,
    MASE = NA_real_, R2 = NA_real_, Bias = NA_real_,
    Relative_variance_error = NA_real_
  ))
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  denom_nse <- sum((obs - mean(obs))^2)
  r <- suppressWarnings(cor(obs, pred))
  beta <- mean(pred) / mean(obs)
  cv_obs <- sd(obs) / mean(obs)
  cv_pred <- sd(pred) / mean(pred)
  gamma <- cv_pred / cv_obs
  tibble(
    N = length(obs),
    RMSE = rmse,
    MAE = mae,
    NSE = ifelse(denom_nse > 0,
                 1 - sum((obs - pred)^2) / denom_nse, NA_real_),
    KGE2012 = ifelse(all(is.finite(c(r, beta, gamma))),
                     1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2),
                     NA_real_),
    PBIAS = ifelse(sum(obs) != 0,
                   100 * sum(pred - obs) / sum(obs), NA_real_),
    RSR = rmse / sd(obs),
    MASE = mae / mase_scale(obs_cal, 12L),
    R2 = ifelse(is.finite(r), r^2, NA_real_),
    Bias = mean(pred - obs),
    Relative_variance_error = ifelse(var(obs) > 0,
                                     (var(pred) - var(obs)) / var(obs),
                                     NA_real_)
  )
}

standardize_train_test <- function(training, new_data, vars) {
  mu <- vapply(training[vars], mean, numeric(1), na.rm = TRUE)
  sg <- vapply(training[vars], sd, numeric(1), na.rm = TRUE)
  sg[!is.finite(sg) | sg == 0] <- 1
  tr <- training; nv <- new_data
  tr[vars] <- Map(function(x, m, s) (x - m) / s, training[vars], mu, sg)
  nv[vars] <- Map(function(x, m, s) (x - m) / s, new_data[vars], mu, sg)
  list(train = tr, new = nv, mean = mu, sd = sg)
}

rolling_origins <- function(df_cal) {
  first_origin <- max(CV_INITIAL_YEARS * 12L, nrow(df_cal) - CV_LAST_MONTHS)
  seq(first_origin, nrow(df_cal) - 1L, by = CV_STEP)
}

safe_rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))


# ------------------------------------------------------------------------------
# 3. Rolling-origin cross-validation for regression and primary MLP models
# ------------------------------------------------------------------------------
cv_linear_candidate <- function(df_cal, predictors) {
  df <- df_cal %>% select(date, freshwater_discharge_m3_s, freshwater_discharge_observed, all_of(predictors))
  idx_cv <- rolling_origins(df)
  errors <- map_dbl(idx_cv, function(i) {
    tr <- df[seq_len(i), ] %>% filter(if_all(all_of(c("freshwater_discharge_m3_s", predictors)), is.finite))
    te <- df[i + 1L, , drop = FALSE]
    if (!te$freshwater_discharge_observed || !all(is.finite(unlist(te[c("freshwater_discharge_m3_s", predictors)])))) return(NA_real_)
    if (nrow(tr) <= length(predictors) + 1L) return(NA_real_)
    fit <- lm(reformulate(predictors, response = "freshwater_discharge_m3_s"), data = tr)
    as.numeric(te$freshwater_discharge_m3_s - predict(fit, newdata = te))
  })
  eligible_flag <- vapply(idx_cv + 1L, function(j) {
    df$freshwater_discharge_observed[j] && all(is.finite(unlist(df[j, c("freshwater_discharge_m3_s", predictors)])))
  }, logical(1))
  cv_total <- sum(eligible_flag)
  tibble(
    CV_RMSE = sqrt(mean(errors^2, na.rm = TRUE)),
    CV_MAE = mean(abs(errors), na.rm = TRUE),
    CV_N = sum(is.finite(errors)),
    CV_total = cv_total,
    CV_success = CV_N / cv_total
  ) %>%
    mutate(
      CV_RMSE = ifelse(CV_success < 0.90, Inf, CV_RMSE),
      CV_MAE = ifelse(CV_success < 0.90, Inf, CV_MAE)
    )
}

# Primary MLP implementation:
# - one hidden layer with logistic sigmoid activation and a linear output unit;
# - nnet/BFGS optimization with a squared-error objective;
# - optional L2 weight-decay regularization through the decay parameter;
# - no frequency or differencing parameter in the MLP specification;
# - no validation-based early stopping and no stopping epoch;
# - optimization ends at BFGS convergence or at maxit = 1000;
# - three fixed-seed random starts per rolling-origin fold and 20 final starts.
# Alternative training loss functions were not evaluated in the primary nnet
# framework.
cv_mlp_candidate <- function(df_cal, predictors, hidden, decay, seed) {
  vars <- c("freshwater_discharge_m3_s", predictors)
  df <- df_cal %>% select(date, freshwater_discharge_observed, all_of(vars))
  idx_cv <- rolling_origins(df)
  errors <- map_dbl(idx_cv, function(i) {
    tr <- df[seq_len(i), ] %>% filter(if_all(all_of(vars), is.finite))
    te <- df[i + 1L, , drop = FALSE]
    if (!te$freshwater_discharge_observed || !all(is.finite(unlist(te[vars])))) return(NA_real_)
    if (nrow(tr) <= length(predictors) + 2L) return(NA_real_)
    z <- standardize_train_test(tr, te, vars)
    preds_rep <- map_dbl(seq_len(MLP_CV_RESTARTS), function(rep_id) {
      set.seed(seed + 10000L * rep_id + i)
      fit <- tryCatch(
        nnet(reformulate(predictors, response = "freshwater_discharge_m3_s"), data = z$train,
             size = hidden, linout = TRUE, decay = decay, maxit = 1000,
             MaxNWts = 10000, trace = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      pred_z <- as.numeric(predict(fit, newdata = z$new))
      pred_z * z$sd[["freshwater_discharge_m3_s"]] + z$mean[["freshwater_discharge_m3_s"]]
    })
    if (mean(is.finite(preds_rep)) < 0.90) return(NA_real_)
    te$freshwater_discharge_m3_s - mean(preds_rep, na.rm = TRUE)
  })
  cv_n <- sum(is.finite(errors))
  eligible_flag <- vapply(idx_cv + 1L, function(j) {
    df$freshwater_discharge_observed[j] && all(is.finite(unlist(df[j, vars])))
  }, logical(1))
  cv_total <- sum(eligible_flag)
  tibble(
    CV_RMSE = sqrt(mean(errors^2, na.rm = TRUE)),
    CV_MAE = mean(abs(errors), na.rm = TRUE),
    CV_N = cv_n,
    CV_total = cv_total,
    CV_success = cv_n / cv_total
  )
}

select_linear_models <- function(df_cal, predictor_sets, period) {
  result <- imap_dfr(predictor_sets, function(preds, name) {
    bind_cols(tibble(Period = period, Model = name,
                     Predictors = paste(preds, collapse = "+")),
              cv_linear_candidate(df_cal, preds))
  }) %>%
    filter(is.finite(CV_RMSE), is.finite(CV_MAE), CV_success >= 0.90) %>%
    arrange(CV_RMSE, CV_MAE)
  if (nrow(result) == 0L) {
    stop("No numerically valid linear model was found for ", period, ".")
  }
  result
}

select_mlp_models <- function(df_cal, predictor_sets, period, model_type) {
  monthly_grid <- tidyr::crossing(
    Model = names(predictor_sets),
    Hidden = c(1L, 2L, 3L, 4L, 5L, 6L, 8L, 10L),
    Decay = c(0, 0.0001, 0.001, 0.01)
  )
  result <- pmap_dfr(monthly_grid, function(Model, Hidden, Decay) {
    preds <- predictor_sets[[Model]]
    bind_cols(
      tibble(Period = period, Type = model_type, Model = Model,
             Predictors = paste(preds, collapse = "+"),
             Hidden = Hidden, Decay = Decay,
             Architecture = paste(length(preds), Hidden, 1, sep = "-")),
      cv_mlp_candidate(df_cal, preds, Hidden, Decay, SEED)
    )
  }) %>%
    filter(is.finite(CV_RMSE), is.finite(CV_MAE), CV_success >= 0.90) %>%
    arrange(CV_RMSE, CV_MAE, Hidden, Decay)
  if (nrow(result) == 0L) {
    stop("No numerically valid MLP configuration was found for ", period,
         " (", model_type, ").")
  }
  result
}


arma_grid <- tidyr::crossing(p = 0:4, q = 0:4) %>% filter(p + q > 0)
arima_grid <- tidyr::crossing(p = 0:5, q = 0:4) # d = 1
sarima_grid <- tidyr::crossing(p = 0:3, q = 0:3, P = 0:2, Q = 0:2) %>%
  filter(p + q + P + Q > 0) # d = 1, D = 0, period = 12

# ------------------------------------------------------------------------------
# 4. Rolling-origin cross-validation for ARMA, ARIMA and SARIMA models
# ------------------------------------------------------------------------------
cv_arima_candidate <- function(y, observed_flag, p, d, q, P = 0L, D = 0L, Q = 0L) {
  idx <- rolling_origins(tibble(x = y))
  errors <- map_dbl(idx, function(i) {
    if (!observed_flag[i + 1L] || !is.finite(y[i + 1L])) return(NA_real_)
    yi <- y[seq_len(i)]
    fit <- tryCatch(
      forecast::Arima(ts(yi, frequency = 12), order = c(p, d, q),
                      seasonal = list(order = c(P, D, Q), period = 12),
                      include.mean = (d == 0), method = "ML"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    y[i + 1L] - as.numeric(forecast(fit, h = 1)$mean)
  })
  tibble(CV_RMSE = sqrt(mean(errors^2, na.rm = TRUE)),
         CV_MAE = mean(abs(errors), na.rm = TRUE),
         CV_N = sum(is.finite(errors)),
         CV_total = sum(observed_flag[idx + 1L] & is.finite(y[idx + 1L])),
         CV_success = CV_N / CV_total) %>%
    mutate(
      CV_RMSE = ifelse(CV_success < 0.90, Inf, CV_RMSE),
      CV_MAE = ifelse(CV_success < 0.90, Inf, CV_MAE)
    )
}

select_arima_grid <- function(y, observed_flag, monthly_grid, model_class, period) {
  result <- pmap_dfr(monthly_grid, function(...) {
    pars <- list(...)
    p <- pars$p; q <- pars$q
    d <- ifelse(model_class == "ARMA", 0L, 1L)
    P <- ifelse("P" %in% names(pars), pars$P, 0L)
    Q <- ifelse("Q" %in% names(pars), pars$Q, 0L)
    met <- cv_arima_candidate(y, observed_flag, p, d, q, P, 0L, Q)
    tibble(Period = period, Class = model_class, p = p, d = d, q = q,
           P = P, D = 0L, Q = Q, s = 12L) %>% bind_cols(met)
  }) %>%
    filter(is.finite(CV_RMSE), is.finite(CV_MAE), CV_success >= 0.90) %>%
    arrange(CV_RMSE, CV_MAE, p + q + P + Q)
  if (nrow(result) == 0L) {
    stop("No numerically valid ", model_class, " model was found for ", period, ".")
  }
  result
}


# ------------------------------------------------------------------------------
# 5. Final model fitting and held-out one-step/recursive prediction
# ------------------------------------------------------------------------------
fit_mlp <- function(train, preds, hidden, decay, seed = SEED) {
  vars <- c("freshwater_discharge_m3_s", preds)
  z <- standardize_train_test(train, train, vars)
  fits <- map(seq_len(MLP_FINAL_RESTARTS), function(rep_id) {
    set.seed(seed + rep_id - 1L)
    tryCatch(
      nnet(reformulate(preds, response = "freshwater_discharge_m3_s"), data = z$train,
           size = hidden, linout = TRUE, decay = decay, maxit = 1000,
           MaxNWts = 10000, trace = FALSE),
      error = function(e) NULL
    )
  })
  fits <- compact(fits)
  if (length(fits) < ceiling(0.90 * MLP_FINAL_RESTARTS)) {
    stop("Fewer than 90% of the final MLP random initializations converged.")
  }
  list(fits = fits, mean = z$mean, sd = z$sd, vars = vars, preds = preds,
       n_success = length(fits), n_requested = MLP_FINAL_RESTARTS)
}

predict_mlp <- function(obj, newdata) {
  z <- newdata
  z[obj$vars] <- Map(function(x, m, s) (x - m) / s,
                     newdata[obj$vars], obj$mean, obj$sd)
  standardized_predictions <- vapply(
    obj$fits,
    function(fit) as.numeric(predict(fit, newdata = z)),
    numeric(nrow(newdata))
  )
  if (is.null(dim(standardized_predictions))) standardized_predictions <- matrix(standardized_predictions, nrow = 1L)
  rowMeans(standardized_predictions) * obj$sd[["freshwater_discharge_m3_s"]] + obj$mean[["freshwater_discharge_m3_s"]]
}

predict_mlp_complete_cases <- function(obj, newdata) {
  out <- rep(NA_real_, nrow(newdata))
  ok <- complete.cases(newdata[, obj$preds, drop = FALSE])
  if (any(ok)) out[ok] <- predict_mlp(obj, newdata[ok, , drop = FALSE])
  out
}

fit_arima_safe <- function(y, pars, context = "") {
  methods <- c("ML", "CSS-ML", "CSS")
  errors <- character()
  for (method in methods) {
    fit <- tryCatch(
      suppressWarnings(
        forecast::Arima(
          ts(y, frequency = 12), order = c(pars$p, pars$d, pars$q),
          seasonal = list(order = c(pars$P, pars$D, pars$Q), period = 12),
          include.mean = (pars$d == 0), method = method
        )
      ),
      error = function(e) {
        errors <<- c(errors, paste0(method, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(fit) && all(is.finite(coef(fit)))) {
      attr(fit, "metodo_usado") <- method
      return(fit)
    }
  }
  stop(
    "ARIMA fitting failed in ", context, " for order ",
    sprintf("(%d,%d,%d)(%d,%d,%d)[12]", pars$p, pars$d, pars$q,
            pars$P, pars$D, pars$Q), ". Attempts: ",
    paste(errors, collapse = " | ")
  )
}

predict_one_step_arima <- function(y_cal, y_val, pars, context = "") {
  fit_cal <- fit_arima_safe(y_cal, pars,
                                  context = paste0(context, ", fixed calibration"))
  historico <- y_cal
  out <- numeric(length(y_val))
  for (i in seq_along(y_val)) {
    updated_fit <- tryCatch(
      suppressWarnings(
        forecast::Arima(ts(historico, frequency = 12), model = fit_cal)
      ),
      error = function(e) stop(
        "Failed to apply the fixed-calibration model in ", context,
        ", one-step origin ", i, ": ", conditionMessage(e)
      )
    )
    out[i] <- as.numeric(forecast(updated_fit, h = 1)$mean)
    historico <- c(historico, y_val[i])
  }
  out
}

predict_recursive_arima <- function(y_cal, h, pars, context = "") {
  fit <- fit_arima_safe(y_cal, pars, context = paste0(context, ", recursive"))
  as.numeric(forecast(fit, h = h)$mean)
}

replace_discharge_lags <- function(row_data, discharge_history) {
  qcols <- grep("^Q[0-9]+$", names(row_data), value = TRUE)
  for (v in qcols) {
    k <- as.integer(sub("Q", "", v))
    row_data[[v]] <- if (length(discharge_history) >= k) discharge_history[length(discharge_history) - k + 1L] else NA_real_
  }
  row_data
}

predict_recursive_regressor <- function(model, val, y_cal, model_type = c("lm", "mlp")) {
  model_type <- match.arg(model_type)
  hist <- y_cal
  out <- numeric(nrow(val))
  preds <- if (model_type == "lm") all.vars(delete.response(terms(model))) else model$preds
  for (i in seq_len(nrow(val))) {
    row_data <- replace_discharge_lags(val[i, , drop = FALSE], hist)
    if (!all(is.finite(unlist(row_data[, preds, drop = FALSE])))) {
      out[i] <- NA_real_
    } else {
      out[i] <- if (model_type == "lm") {
        as.numeric(predict(model, newdata = row_data))
      } else {
        predict_mlp(model, row_data)
      }
    }
    hist <- c(hist, out[i])
  }
  out
}

# Primary held-out predictions are sequential one-step-ahead predictions.
# Observed antecedent freshwater discharge is used at each validation month.
# Precipitation-informed models receive the observed historical precipitation
# values at their required lags; precipitation is never recursively predicted.
# Completely recursive discharge trajectories are evaluated separately below.
evaluate_period <- function(period, cal_end, val_start, selections) {
  cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  val <- lagged_data %>% filter(date >= val_start, date <= END_DATE)
  y_cal <- cal$freshwater_discharge_m3_s; y_val <- val$freshwater_discharge_m3_s
  y_obs_avaliacao <- val$freshwater_discharge_original_m3_s
  output_rows <- list()

  pred_sn_1 <- val$Q12
  pred_sn_r <- rep(tail(y_cal, 12L), length.out = nrow(val))
  output_rows[["SeasonalNaive"]] <- tibble(date = val$date, Observed = y_obs_avaliacao,
                                       OneStep = pred_sn_1, Recursive = pred_sn_r,
                                       Model = "Seasonal naive", Family = "Benchmark")

  for (model_class in c("ARMA", "ARIMA", "SARIMA")) {
    pars <- selections$ar[[model_class]][1, ]
    message("Evaluating ", period, " — ", model_class, "...")
    p1 <- predict_one_step_arima(y_cal, y_val, pars,
                                  context = paste(period, model_class))
    pr <- predict_recursive_arima(y_cal, nrow(val), pars,
                                   context = paste(period, model_class))
    name <- sprintf("%s(%d,%d,%d)(%d,%d,%d)[12]", model_class,
                    pars$p, pars$d, pars$q, pars$P, pars$D, pars$Q)
    output_rows[[model_class]] <- tibble(date = val$date, Observed = y_obs_avaliacao,
                               OneStep = p1, Recursive = pr,
                               Model = name, Family = model_class)
  }

  best_lm <- selections$lm[1, ]
  preds_lm <- strsplit(best_lm$Predictors, "\\+")[[1]]
  cal_lm <- cal %>% select(freshwater_discharge_m3_s, all_of(preds_lm)) %>% drop_na()
  fit_lm <- lm(reformulate(preds_lm, response = "freshwater_discharge_m3_s"), data = cal_lm)
  p1_lm <- as.numeric(predict(fit_lm, newdata = val))
  pr_lm <- predict_recursive_regressor(fit_lm, val, cal$freshwater_discharge_m3_s, "lm")
  output_rows[["Linear"]] <- tibble(date = val$date, Observed = y_obs_avaliacao,
                                OneStep = p1_lm, Recursive = pr_lm,
                                Model = paste0("Linear [", best_lm$Model, "]"),
                                Family = "Linear regression")

  bu <- selections$mlp_uni[1, ]
  preds_u <- strsplit(bu$Predictors, "\\+")[[1]]
  fit_u <- fit_mlp(cal %>% select(freshwater_discharge_m3_s, all_of(preds_u)) %>% drop_na(),
                   preds_u, bu$Hidden, bu$Decay)
  p1_u <- predict_mlp_complete_cases(fit_u, val)
  pr_u <- predict_recursive_regressor(fit_u, val, cal$freshwater_discharge_m3_s, "mlp")
  output_rows[["MLP_uni"]] <- tibble(date = val$date, Observed = y_obs_avaliacao,
                                 OneStep = p1_u, Recursive = pr_u,
                                 Model = paste0("MLP univariate [", bu$Model, "; ", bu$Architecture, "]"),
                                 Family = "MLP univariate")

  bm <- selections$mlp_multi[1, ]
  preds_m <- strsplit(bm$Predictors, "\\+")[[1]]
  fit_m <- fit_mlp(cal %>% select(freshwater_discharge_m3_s, all_of(preds_m)) %>% drop_na(),
                   preds_m, bm$Hidden, bm$Decay)
  p1_m <- predict_mlp_complete_cases(fit_m, val)
  pr_m <- predict_recursive_regressor(fit_m, val, cal$freshwater_discharge_m3_s, "mlp")
  output_rows[["MLP_multi"]] <- tibble(date = val$date, Observed = y_obs_avaliacao,
                                   OneStep = p1_m, Recursive = pr_m,
                                   Model = paste0("MLP multivariate [", bm$Model, "; ", bm$Architecture, "]"),
                                   Family = "MLP multivariate")

  bind_rows(output_rows) %>% mutate(Period = period)
}


# ------------------------------------------------------------------------------
# 6. Model selection and held-out evaluation for T1 and T2
# ------------------------------------------------------------------------------
run_model_selection <- function(period, cal_end) {
  cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  message("Selecting models for ", period, "...")
  list(
    lm = select_linear_models(cal, multivariate_sets, period),
    mlp_uni = select_mlp_models(cal, univariate_sets, period, "Univariate"),
    mlp_multi = select_mlp_models(cal, multivariate_sets, period, "Multivariate"),
    ar = list(
      ARMA = select_arima_grid(cal$freshwater_discharge_m3_s, cal$freshwater_discharge_observed, arma_grid, "ARMA", period),
      ARIMA = select_arima_grid(cal$freshwater_discharge_m3_s, cal$freshwater_discharge_observed, arima_grid, "ARIMA", period),
      SARIMA = select_arima_grid(cal$freshwater_discharge_m3_s, cal$freshwater_discharge_observed, sarima_grid, "SARIMA", period)
    )
  )
}

selected_T1 <- run_model_selection("T1", T1_CAL_END)
selected_T2 <- run_model_selection("T2", T2_CAL_END)


write_csv(selected_T1$lm, file.path(OUTPUT_DIR, "cv_linear_T1.csv"))
write_csv(selected_T2$lm, file.path(OUTPUT_DIR, "cv_linear_T2.csv"))
write_csv(selected_T1$mlp_uni, file.path(OUTPUT_DIR, "cv_mlp_univariate_T1.csv"))
write_csv(selected_T2$mlp_uni, file.path(OUTPUT_DIR, "cv_mlp_univariate_T2.csv"))
write_csv(selected_T1$mlp_multi, file.path(OUTPUT_DIR, "cv_mlp_multivariate_T1.csv"))
write_csv(selected_T2$mlp_multi, file.path(OUTPUT_DIR, "cv_mlp_multivariate_T2.csv"))
walk2(selected_T1$ar, names(selected_T1$ar), ~write_csv(.x, file.path(OUTPUT_DIR, paste0("cv_", tolower(.y), "_T1.csv"))))
walk2(selected_T2$ar, names(selected_T2$ar), ~write_csv(.x, file.path(OUTPUT_DIR, paste0("cv_", tolower(.y), "_T2.csv"))))

audit_model_n <- function(period, cal_end, val_start, selections) {
  cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  val <- lagged_data %>% filter(date >= val_start, date <= END_DATE)
  linhas_preditor <- function(family, row_data) {
    preds <- strsplit(row_data$Predictors, "\\+")[[1]]
    tibble(
      Period = period, Family = family, Model = row_data$Model,
      Predictors = row_data$Predictors,
      N_cal_calendar = nrow(cal),
      N_cal_complete = sum(complete.cases(cal[, c("freshwater_discharge_m3_s", preds), drop = FALSE])),
      N_val_calendar = nrow(val),
      N_val_observed_target = sum(val$freshwater_discharge_observed),
      N_val_complete = sum(val$freshwater_discharge_observed &
                             complete.cases(val[, preds, drop = FALSE]))
    )
  }
  ar <- map_dfr(c("ARMA", "ARIMA", "SARIMA"), function(model_class) {
    row_data <- selections$ar[[model_class]][1, ]
    tibble(
      Period = period, Family = model_class,
      Model = sprintf("(%d,%d,%d)(%d,%d,%d)[12]", row_data$p, row_data$d,
                      row_data$q, row_data$P, row_data$D, row_data$Q),
      Predictors = "monthly discharge series",
      N_cal_calendar = nrow(cal), N_cal_complete = sum(cal$freshwater_discharge_observed),
      N_val_calendar = nrow(val),
      N_val_observed_target = sum(val$freshwater_discharge_observed),
      N_val_complete = sum(val$freshwater_discharge_observed)
    )
  })
  bind_rows(
    ar,
    linhas_preditor("Linear regression", selections$lm[1, ]),
    linhas_preditor("MLP univariate", selections$mlp_uni[1, ]),
    linhas_preditor("MLP multivariate", selections$mlp_multi[1, ])
  )
}

selected_model_n_audit <- bind_rows(
  audit_model_n("T1", T1_CAL_END, T1_VAL_START, selected_T1),
  audit_model_n("T2", T2_CAL_END, T2_VAL_START, selected_T2)
)
write_csv(selected_model_n_audit,
          file.path(OUTPUT_DIR, "selected_model_sample_size_audit.csv"))

predictions_T1 <- evaluate_period("T1", T1_CAL_END, T1_VAL_START, selected_T1)
predictions_T2 <- evaluate_period("T2", T2_CAL_END, T2_VAL_START, selected_T2)
predictions <- bind_rows(predictions_T1, predictions_T2)
write_csv(predictions, file.path(OUTPUT_DIR, "all_model_predictions.csv"))


# ------------------------------------------------------------------------------
# 7. Validation performance and common-period summaries
# ------------------------------------------------------------------------------
calibration_observed <- list(
  T1 = monthly_data %>% filter(date <= T1_CAL_END) %>% pull(freshwater_discharge_original_m3_s),
  T2 = monthly_data %>% filter(date <= T2_CAL_END) %>% pull(freshwater_discharge_original_m3_s)
)

performance_long <- predictions %>%
  pivot_longer(c(OneStep, Recursive), names_to = "Evaluation", values_to = "Predicted") %>%
  group_by(Period, Family, Model, Evaluation) %>%
  group_modify(~compute_metrics(.x$Observed, .x$Predicted, calibration_observed[[as.character(.y$Period)]])) %>%
  ungroup() %>%
  arrange(Evaluation, Period, RMSE)
write_csv(performance_long, file.path(OUTPUT_DIR, "validation_performance_metrics.csv"))

common_period_metrics <- predictions %>%
  filter(date >= T2_VAL_START) %>%
  pivot_longer(c(OneStep, Recursive), names_to = "Evaluation", values_to = "Predicted") %>%
  group_by(Period, Family, Model, Evaluation) %>%
  group_modify(~compute_metrics(.x$Observed, .x$Predicted, calibration_observed[[as.character(.y$Period)]])) %>%
  ungroup() %>%
  arrange(Evaluation, Family, Period)
write_csv(common_period_metrics, file.path(OUTPUT_DIR, "common_post2012_performance_metrics.csv"))

dm_compare <- function(df, family_a, family_b, evaluation = "OneStep") {
  a <- df %>% filter(Family == family_a) %>% select(date, Observed, A = all_of(evaluation))
  b <- df %>% filter(Family == family_b) %>% select(date, B = all_of(evaluation))
  z <- inner_join(a, b, by = "date") %>% drop_na()
  e1 <- z$Observed - z$A; e2 <- z$Observed - z$B
  test_result <- tryCatch(forecast::dm.test(e1, e2, alternative = "two.sided", h = 1, power = 1),
                    error = function(e) NULL)
  tibble(Model_A = family_a, Model_B = family_b, Evaluation = evaluation,
         N = nrow(z), DM = ifelse(is.null(test_result), NA, unname(test_result$statistic)),
         P_value = ifelse(is.null(test_result), NA, test_result$p.value))
}

dm_resultados <- bind_rows(lapply(split(predictions, predictions$Period), function(df) {
  bind_rows(
    dm_compare(df, "ARIMA", "MLP univariate"),
    dm_compare(df, "Linear regression", "MLP multivariate")
  ) %>% mutate(Period = unique(df$Period))
}))
write_csv(dm_resultados, file.path(OUTPUT_DIR, "selected_model_dm_tests.csv"))

paired_selected_model_test <- function(df, family_a, family_b, seed_offset = 0L) {
  a <- df %>% filter(Family == family_a) %>%
    select(date, Observed, A = OneStep)
  b <- df %>% filter(Family == family_b) %>%
    select(date, B = OneStep)
  z <- inner_join(a, b, by = "date") %>%
    filter(if_all(c(Observed, A, B), is.finite)) %>%
    arrange(date)
  e_a <- z$Observed - z$A
  e_b <- z$Observed - z$B
  dm <- tryCatch(
    forecast::dm.test(e_a, e_b, alternative = "two.sided", h = 1, power = 1),
    error = function(e) NULL
  )
  delta_abs <- abs(e_b) - abs(e_a)
  ci <- moving_block_mean_ci(delta_abs, B = BOOT_B,
                             block_length = BOOT_BLOCK_LENGTH,
                             seed = SEED + seed_offset)
  tibble(
    Model_A = family_a,
    Model_B = family_b,
    N = nrow(z),
    RMSE_A = sqrt(mean(e_a^2)), RMSE_B = sqrt(mean(e_b^2)),
    MAE_A = mean(abs(e_a)), MAE_B = mean(abs(e_b)),
    Delta_MAE_B_minus_A = mean(delta_abs),
    Delta_MAE_CI_low = ci[["lower"]],
    Delta_MAE_CI_high = ci[["upper"]],
    DM_statistic = ifelse(is.null(dm), NA_real_, unname(dm$statistic)),
    DM_p_value = ifelse(is.null(dm), NA_real_, dm$p.value)
  )
}

selected_h1_tests <- bind_rows(lapply(split(predictions, predictions$Period), function(df) {
  p <- unique(df$Period)
  bind_rows(
    paired_selected_model_test(df, "ARMA", "MLP univariate", 101L),
    paired_selected_model_test(df, "ARIMA", "MLP univariate", 102L),
    paired_selected_model_test(df, "SARIMA", "MLP univariate", 103L),
    paired_selected_model_test(df, "Linear regression", "MLP multivariate", 104L)
  ) %>% mutate(Period = p)
})) %>%
  mutate(DM_p_adjusted_BH = p.adjust(DM_p_value, method = "BH"))

write_csv(selected_h1_tests,
          file.path(OUTPUT_DIR, "selected_model_paired_error_tests.csv"))


# ------------------------------------------------------------------------------
# 8. Structural-change, trend, stationarity and rolling-autocorrelation analyses
# ------------------------------------------------------------------------------
annual_series <- monthly_data %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    freshwater_discharge_m3_s = mean(freshwater_discharge_original_m3_s, na.rm = TRUE),
    N_calendar_months = n(),
    N_observed = sum(freshwater_discharge_observed),
    .groups = "drop"
  ) %>%
  filter(N_calendar_months == 12L, N_observed >= 10L, is.finite(freshwater_discharge_m3_s))

extract_test <- function(name, obj, estimate = NA_real_) {
  tibble(Test = name, Statistic = unname(obj$statistic)[1], Estimate = estimate,
         P_value = obj$p.value, Method = obj$method)
}

sen_obj <- trend::sens.slope(annual_series$freshwater_discharge_m3_s)
structural_trend_tests <- bind_rows(
  extract_test("Pettitt", trend::pettitt.test(annual_series$freshwater_discharge_m3_s)),
  extract_test("Buishand U", trend::bu.test(annual_series$freshwater_discharge_m3_s)),
  extract_test("SNHT", trend::snh.test(annual_series$freshwater_discharge_m3_s)),
  extract_test("Mann-Kendall", trend::mk.test(annual_series$freshwater_discharge_m3_s)),
  extract_test("Sen slope", sen_obj, estimate = unname(sen_obj$estimates)[1])
)
write_csv(structural_trend_tests, file.path(OUTPUT_DIR, "structural_change_and_trend_tests.csv"))

rolling_window <- 120L
rolling_acf <- map_dfr(seq(rolling_window, nrow(monthly_data)), function(i) {
  x <- monthly_data$freshwater_discharge_original_m3_s[(i - rolling_window + 1L):i]
  ac_lag <- function(z, k) {
    a <- z[(k + 1L):length(z)]; b <- z[1L:(length(z) - k)]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 3L) NA_real_ else cor(a[ok], b[ok])
  }
  tibble(Window_end_date = monthly_data$date[i],
         ACF_lag1 = ac_lag(x, 1L), ACF_lag12 = ac_lag(x, 12L),
         Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE),
         N_observed = sum(is.finite(x)))
})
write_csv(rolling_acf, file.path(OUTPUT_DIR, "rolling_10year_acf.csv"))

annual_ts <- ts(annual_series$freshwater_discharge_m3_s, start = min(annual_series$year), frequency = 1)
adf_obj <- tseries::adf.test(annual_ts, alternative = "stationary")
kpss_obj <- tseries::kpss.test(annual_ts, null = "Level")
stationarity_results <- tibble(
  Test = c("ADF", "KPSS"),
  Null_hypothesis = c("Unit root / non-stationary", "Level stationary"),
  Statistic = c(unname(adf_obj$statistic), unname(kpss_obj$statistic)),
  P_value = c(adf_obj$p.value, kpss_obj$p.value),
  Suggested_d = c(
    forecast::ndiffs(annual_ts, test = "adf"),
    forecast::ndiffs(annual_ts, test = "kpss")
  )
)
write_csv(stationarity_results, file.path(OUTPUT_DIR, "stationarity_tests.csv"))


# ------------------------------------------------------------------------------
# 9. Fold-level diagnostics for selected MLP configurations
# ------------------------------------------------------------------------------
cv_mlp_detailed <- function(df_cal, selection, period, family) {
  preds <- strsplit(selection$Predictors, "\\+")[[1]]
  vars <- c("freshwater_discharge_m3_s", preds)
  df <- df_cal %>% select(date, freshwater_discharge_observed, all_of(vars))
  idx_cv <- rolling_origins(df)

  map_dfr(seq_along(idx_cv), function(j) {
    i <- idx_cv[j]
    tr_grade <- df[seq_len(i), ]
    tr <- tr_grade %>% filter(if_all(all_of(vars), is.finite))
    te <- df[i + 1L, , drop = FALSE]
    if (!te$freshwater_discharge_observed || !all(is.finite(unlist(te[vars])))) {
      return(tibble(Fold = j, Origin = tr_grade$date[nrow(tr_grade)], Target = te$date,
                    Observed = NA_real_, Predicted = NA_real_,
                    Error = NA_real_, Absolute_error = NA_real_))
    }
    z <- standardize_train_test(tr, te, vars)
    pred_reps <- map_dbl(seq_len(MLP_CV_RESTARTS), function(rep_id) {
      set.seed(SEED + 10000L * rep_id + i)
      fit <- tryCatch(
        nnet(reformulate(preds, response = "freshwater_discharge_m3_s"), data = z$train,
             size = selection$Hidden, linout = TRUE, decay = selection$Decay,
             maxit = 1000, MaxNWts = 10000, trace = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      pred_z <- as.numeric(predict(fit, newdata = z$new))
      pred_z * z$sd[["freshwater_discharge_m3_s"]] + z$mean[["freshwater_discharge_m3_s"]]
    })
    if (mean(is.finite(pred_reps)) < 0.90) {
      pred <- NA_real_
    } else {
      pred <- mean(pred_reps, na.rm = TRUE)
    }
    error_value <- te$freshwater_discharge_m3_s - pred
    tibble(Fold = j, Origin = tr_grade$date[nrow(tr_grade)], Target = te$date,
           Observed = te$freshwater_discharge_m3_s, Predicted = pred, Error = error_value,
           Absolute_error = abs(error_value))
  }) %>%
    mutate(Period = period, Family = family,
           Model = selection$Model, Architecture = selection$Architecture)
}

cv_mlp_folds <- bind_rows(
  cv_mlp_detailed(lagged_data %>% filter(date <= T1_CAL_END),
                   selected_T1$mlp_uni[1, ], "T1", "MLP univariate"),
  cv_mlp_detailed(lagged_data %>% filter(date <= T1_CAL_END),
                   selected_T1$mlp_multi[1, ], "T1", "MLP multivariate"),
  cv_mlp_detailed(lagged_data %>% filter(date <= T2_CAL_END),
                   selected_T2$mlp_uni[1, ], "T2", "MLP univariate"),
  cv_mlp_detailed(lagged_data %>% filter(date <= T2_CAL_END),
                   selected_T2$mlp_multi[1, ], "T2", "MLP multivariate")
)
write_csv(cv_mlp_folds, file.path(OUTPUT_DIR, "selected_mlp_cv_folds.csv"))


# ------------------------------------------------------------------------------
# 10. Multi-horizon sensitivity analysis
# ------------------------------------------------------------------------------
prepare_horizon_models <- function(cal_end, selections, period) {
  cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)

  fits_ar <- map(c("ARMA", "ARIMA", "SARIMA"), function(model_class) {
    pars <- selections$ar[[model_class]][1, ]
    list(
      Family = model_class,
      pars = pars,
      fit = fit_arima_safe(cal$freshwater_discharge_m3_s, pars,
                                 context = paste(period, model_class, "multi-horizon"))
    )
  })
  names(fits_ar) <- c("ARMA", "ARIMA", "SARIMA")

  best_lm <- selections$lm[1, ]
  preds_lm <- strsplit(best_lm$Predictors, "\\+")[[1]]
  fit_lm_obj <- lm(reformulate(preds_lm, response = "freshwater_discharge_m3_s"),
                   data = cal %>% select(freshwater_discharge_m3_s, all_of(preds_lm)) %>% drop_na())

  bu <- selections$mlp_uni[1, ]
  preds_u <- strsplit(bu$Predictors, "\\+")[[1]]
  fit_u_obj <- fit_mlp(cal %>% select(freshwater_discharge_m3_s, all_of(preds_u)) %>% drop_na(),
                       preds_u, bu$Hidden, bu$Decay,
                       seed = SEED + ifelse(period == "T1", 1000L, 2000L))

  bm <- selections$mlp_multi[1, ]
  preds_m <- strsplit(bm$Predictors, "\\+")[[1]]
  fit_m_obj <- fit_mlp(cal %>% select(freshwater_discharge_m3_s, all_of(preds_m)) %>% drop_na(),
                       preds_m, bm$Hidden, bm$Decay,
                       seed = SEED + ifelse(period == "T1", 3000L, 4000L))

  list(cal = cal, ar = fits_ar, lm = fit_lm_obj,
       mlp_uni = fit_u_obj, mlp_multi = fit_m_obj)
}

evaluate_horizons_period <- function(period, val_start, objects) {
  first_origin_index <- max(which(lagged_data$date < val_start))
  origins <- seq(first_origin_index, nrow(lagged_data) - 1L,
                 by = HORIZON_STEP)

  map_dfr(seq_along(origins), function(j) {
    origin_index <- origins[j]
    h_max <- min(max(HORIZONS), nrow(lagged_data) - origin_index)
    hs <- HORIZONS[HORIZONS <= h_max]
    if (length(hs) == 0L) return(tibble())

    future_data <- lagged_data[(origin_index + 1L):(origin_index + h_max), , drop = FALSE]
    discharge_history <- lagged_data$freshwater_discharge_m3_s[seq_len(origin_index)]
    origin_date <- lagged_data$date[origin_index]

    prediction_list <- list()

    prediction_list[["Benchmark"]] <- tibble(
      Family = "Benchmark",
      Predicted = rep(tail(discharge_history, 12L), length.out = h_max)
    )

    for (model_class in c("ARMA", "ARIMA", "SARIMA")) {
      updated_fit <- tryCatch(
        suppressWarnings(forecast::Arima(ts(discharge_history, frequency = 12),
                                         model = objects$ar[[model_class]]$fit)),
        error = function(e) NULL
      )
      pred <- if (is.null(updated_fit)) rep(NA_real_, h_max) else
        as.numeric(forecast(updated_fit, h = h_max)$mean)
      prediction_list[[model_class]] <- tibble(Family = model_class, Predicted = pred)
    }

    prediction_list[["Linear regression"]] <- tibble(
      Family = "Linear regression",
      Predicted = predict_recursive_regressor(objects$lm, future_data,
                                                discharge_history, "lm")
    )
    prediction_list[["MLP univariate"]] <- tibble(
      Family = "MLP univariate",
      Predicted = predict_recursive_regressor(objects$mlp_uni, future_data,
                                                discharge_history, "mlp")
    )
    prediction_list[["MLP multivariate"]] <- tibble(
      Family = "MLP multivariate",
      Predicted = predict_recursive_regressor(objects$mlp_multi, future_data,
                                                discharge_history, "mlp")
    )

    bind_rows(prediction_list) %>%
      group_by(Family) %>%
      mutate(Horizon = seq_len(n()),
             date = future_data$date,
             Observed = ifelse(future_data$freshwater_discharge_observed,
                               future_data$freshwater_discharge_original_m3_s, NA_real_)) %>%
      ungroup() %>%
      filter(Horizon %in% hs) %>%
      mutate(Period = period, Origin_id = j, Origin = origin_date)
  })
}

horizon_models_T1 <- prepare_horizon_models(T1_CAL_END, selected_T1, "T1")
horizon_models_T2 <- prepare_horizon_models(T2_CAL_END, selected_T2, "T2")

horizon_predictions <- bind_rows(
  evaluate_horizons_period("T1", T1_VAL_START, horizon_models_T1),
  evaluate_horizons_period("T2", T2_VAL_START, horizon_models_T2)
)
write_csv(horizon_predictions,
          file.path(OUTPUT_DIR, "multihorizon_predictions.csv"))

horizon_metrics <- horizon_predictions %>%
  filter(is.finite(Observed), is.finite(Predicted)) %>%
  group_by(Period, Family, Horizon) %>%
  summarise(
    N = n(),
    RMSE = sqrt(mean((Observed - Predicted)^2)),
    MAE = mean(abs(Observed - Predicted)),
    Bias = mean(Predicted - Observed),
    PBIAS = 100 * sum(Predicted - Observed) / sum(Observed),
    NSE = ifelse(sum((Observed - mean(Observed))^2) > 0,
                 1 - sum((Observed - Predicted)^2) /
                   sum((Observed - mean(Observed))^2), NA_real_),
    .groups = "drop"
  )
write_csv(horizon_metrics,
          file.path(OUTPUT_DIR, "multihorizon_performance.csv"))


family_order <- c("Benchmark", "ARMA", "ARIMA", "SARIMA",
                    "Linear regression", "MLP univariate", "MLP multivariate")

tema_artigo <- theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = "grey45"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(colour = "black")
  )

fig3_data <- performance_long %>%
  filter(Evaluation == "OneStep") %>%
  mutate(Family = factor(Family, levels = family_order)) %>%
  select(Period, Family, RMSE, MAE, RSR, MASE) %>%
  pivot_longer(c(RMSE, MAE, RSR, MASE),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("RMSE", "MAE", "RSR", "MASE")))

fig3 <- ggplot(fig3_data, aes(Family, Value, colour = Period, shape = Period)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(T1 = "#2166AC", T2 = "#B2182B")) +
  labs(x = NULL, y = NULL, colour = "Configuration", shape = "Configuration") +
  tema_artigo +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "model_performance_metrics.png"),
       fig3, width = 9.2, height = 7.2, dpi = 600)

fig4 <- predictions %>%
  filter(Family %in% c("MLP univariate", "MLP multivariate")) %>%
  transmute(Period, Family, Observed, Predicted = OneStep) %>%
  filter(is.finite(Observed), is.finite(Predicted)) %>%
  ggplot(aes(Observed, Predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = "grey35", linewidth = 0.6) +
  geom_point(alpha = 0.70, size = 1.5, colour = "#2C7FB8") +
  facet_grid(Family ~ Period) +
  coord_equal() +
  labs(x = expression(Observed~discharge~(m^3~s^{-1})),
       y = expression(Predicted~discharge~(m^3~s^{-1}))) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "mlp_observed_predicted.png"),
       fig4, width = 8.5, height = 7.0, dpi = 600)

cv_medians <- cv_mlp_folds %>%
  group_by(Period, Family) %>%
  summarise(Median = median(Absolute_error, na.rm = TRUE), .groups = "drop")

fig5 <- cv_mlp_folds %>%
  filter(is.finite(Absolute_error)) %>%
  ggplot(aes(Target, Absolute_error)) +
  geom_line(colour = "grey55", linewidth = 0.35) +
  geom_point(colour = "#2C7FB8", size = 1.4) +
  geom_hline(data = cv_medians, aes(yintercept = Median),
             linetype = 2, colour = "#B2182B", linewidth = 0.55) +
  facet_grid(Family ~ Period, scales = "free_x") +
  labs(x = "Target month", y = expression(Absolute~error~(m^3~s^{-1}))) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "rolling_origin_mlp.png"),
       fig5, width = 9.2, height = 6.8, dpi = 600)

fig6 <- predictions %>%
  transmute(Period, Family, Observed, Predicted = OneStep) %>%
  filter(is.finite(Observed), is.finite(Predicted)) %>%
  mutate(Family = factor(Family, levels = family_order)) %>%
  ggplot(aes(Observed, Predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = "grey35", linewidth = 0.45) +
  geom_point(alpha = 0.62, size = 0.9, colour = "#2C7FB8") +
  facet_grid(Family ~ Period) +
  coord_equal() +
  labs(x = expression(Observed~discharge~(m^3~s^{-1})),
       y = expression(Predicted~discharge~(m^3~s^{-1}))) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "all_models_observed_predicted.png"),
       fig6, width = 8.8, height = 13.5, dpi = 600)
fig6 <- predictions %>%
  transmute(
    Period,
    Family,
    Observed,
    Predicted = OneStep
  ) %>%
  filter(
    is.finite(Observed),
    is.finite(Predicted)
  ) %>%
  mutate(
    Family = factor(Family, levels = family_order)
  ) %>%
  ggplot(aes(Observed, Predicted)) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    colour = "grey35",
    linewidth = 0.45
  ) +
  
  geom_point(
    alpha = 0.62,
    size = 0.9,
    colour = "#2C7FB8"
  ) +
  
  facet_grid(Period ~ Family) +
  
  coord_equal() +
  
  labs(
    x = expression(Observed~discharge~(m^3~s^{-1})),
    y = expression(Predicted~discharge~(m^3~s^{-1}))
  ) +
  
  tema_artigo +
  
  theme(
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 9),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 10),
    panel.spacing = unit(0.35, "lines")
  )


ggsave(
  file.path(
    OUTPUT_DIR,
    "all_models_observed_predicted_landscape.png"
  ),
  fig6,
  width = 14,
  height = 5.5,
  units = "in",
  dpi = 600
)


fig7 <- predictions %>%
  transmute(Period, Family, date, Observed, Predicted = OneStep) %>%
  mutate(Family = factor(Family, levels = family_order)) %>%
  pivot_longer(c(Observed, Predicted), names_to = "Series",
               values_to = "freshwater_discharge_m3_s") %>%
  ggplot(aes(date, freshwater_discharge_m3_s, colour = Series, linewidth = Series)) +
  geom_line(na.rm = TRUE) +
  facet_grid(Family ~ Period, scales = "fixed") +
  scale_colour_manual(values = c(Observed = "black", Predicted = "#2C7FB8")) +
  scale_linewidth_manual(values = c(Observed = 0.62, Predicted = 0.45)) +
  labs(x = NULL, y = expression(freshwater_discharge_m3_s~(m^3~s^{-1})),
       colour = NULL, linewidth = NULL) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "all_models_time_series_one_step.png"),
       fig7, width = 10.2, height = 13.5, dpi = 600)

figS_recursive <- predictions %>%
  transmute(Period, Family, date, Observed, Predicted = Recursive) %>%
  mutate(Family = factor(Family, levels = family_order)) %>%
  pivot_longer(c(Observed, Predicted), names_to = "Series",
               values_to = "freshwater_discharge_m3_s") %>%
  ggplot(aes(date, freshwater_discharge_m3_s, colour = Series, linewidth = Series)) +
  geom_line(na.rm = TRUE) +
  facet_grid(Family ~ Period, scales = "free_y") +
  scale_colour_manual(values = c(Observed = "black", Predicted = "#D95F02")) +
  scale_linewidth_manual(values = c(Observed = 0.62, Predicted = 0.45)) +
  labs(x = NULL, y = expression(freshwater_discharge_m3_s~(m^3~s^{-1})),
       colour = NULL, linewidth = NULL) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "recursive_time_series.png"),
       figS_recursive, width = 10.2, height = 13.5, dpi = 600)

figS_horizon <- horizon_metrics %>%
  select(Period, Family, Horizon, RMSE, MAE) %>%
  pivot_longer(c(RMSE, MAE), names_to = "Metric", values_to = "Value") %>%
  mutate(Family = factor(Family, levels = family_order)) %>%
  ggplot(aes(Horizon, Value, colour = Family, group = Family)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.7) +
  facet_grid(Metric ~ Period, scales = "free_y") +
  scale_x_continuous(breaks = HORIZONS) +
  labs(x = "Forecast horizon (months)", y = NULL, colour = "Model") +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "multihorizon_performance_plot.png"),
       figS_horizon, width = 9.5, height = 7.0, dpi = 600)

q10 <- quantile(monthly_data$freshwater_discharge_original_m3_s, 0.10, na.rm = TRUE)
q90 <- quantile(monthly_data$freshwater_discharge_original_m3_s, 0.90, na.rm = TRUE)
figS_series <- monthly_data %>%
  ggplot(aes(date, freshwater_discharge_original_m3_s)) +
  geom_line(colour = "black", linewidth = 0.45, na.rm = TRUE) +
  geom_hline(yintercept = c(q10, q90), linetype = 3,
             colour = c("#B2182B", "#2166AC")) +
  geom_vline(xintercept = as.numeric(c(T1_VAL_START, T2_VAL_START)),
             linetype = c(2, 1), colour = c("grey40", "#B2182B")) +
  labs(x = NULL, y = expression(freshwater_discharge_m3_s~(m^3~s^{-1}))) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "observed_series_partitions.png"),
       figS_series, width = 10, height = 4.8, dpi = 600)

climatology <- monthly_data %>%
  filter(freshwater_discharge_observed) %>%
  group_by(month) %>%
  summarise(Mean = mean(freshwater_discharge_original_m3_s), SD = sd(freshwater_discharge_original_m3_s), .groups = "drop")
figS_clim <- ggplot(climatology, aes(month, Mean)) +
  geom_ribbon(aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
              fill = "grey75", alpha = 0.6) +
  geom_line(colour = "#2C7FB8", linewidth = 0.8) +
  geom_point(colour = "#2C7FB8", size = 1.8) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(x = NULL, y = expression(Mean~monthly~discharge~(m^3~s^{-1}))) +
  tema_artigo
ggsave(file.path(OUTPUT_DIR, "monthly_climatology.png"),
       figS_clim, width = 7.5, height = 4.5, dpi = 600)

fig_acf <- rolling_acf %>%
  select(Window_end_date, ACF_lag1, ACF_lag12) %>%
  pivot_longer(-Window_end_date, names_to = "Lag", values_to = "ACF") %>%
  ggplot(aes(Window_end_date, ACF, colour = Lag)) + geom_line(linewidth = 0.7) +
  geom_vline(xintercept = as.numeric(T2_VAL_START), linetype = 2) +
  theme_bw() + labs(x = NULL, y = "10-year moving-window ACF", colour = NULL)
ggsave(file.path(OUTPUT_DIR, "rolling_10year_acf.png"), fig_acf,
       width = 8, height = 4.5, dpi = 600)


selected_predictor_summary <- bind_rows(
  selected_T1$lm[1, ] %>% mutate(Family = "Linear regression"),
  selected_T2$lm[1, ] %>% mutate(Family = "Linear regression"),
  selected_T1$mlp_uni[1, ] %>% mutate(Family = "MLP univariate"),
  selected_T2$mlp_uni[1, ] %>% mutate(Family = "MLP univariate"),
  selected_T1$mlp_multi[1, ] %>% mutate(Family = "MLP multivariate"),
  selected_T2$mlp_multi[1, ] %>% mutate(Family = "MLP multivariate")
)
write_csv(selected_predictor_summary, file.path(OUTPUT_DIR, "selected_predictor_models.csv"))

selected_ar_summary <- bind_rows(
  map_dfr(selected_T1$ar, ~.x[1, ]),
  map_dfr(selected_T2$ar, ~.x[1, ])
)
write_csv(selected_ar_summary, file.path(OUTPUT_DIR, "selected_stochastic_models.csv"))

capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "session_info.txt"))
writeLines(c(
  paste("Seed:", SEED),
  paste("Final MLP initializations:", MLP_FINAL_RESTARTS),
  paste("T1 data file:", normalizePath(T1_DATA_FILE)),
  paste("T2 data file:", normalizePath(T2_DATA_FILE)),
  paste("date range:", min(monthly_data$date), "to", max(monthly_data$date)),
  paste("Primary evaluation:", PRIMARY_EVALUATION),
  paste("CV initial years:", CV_INITIAL_YEARS),
  paste("CV last months:", CV_LAST_MONTHS),
  paste("CV step:", CV_STEP),
  paste("Bootstrap replicates:", BOOT_B),
  paste("Bootstrap block length:", BOOT_BLOCK_LENGTH),
  "Primary MLP hidden activation: logistic sigmoid",
  "Primary MLP output activation: linear",
  "Primary MLP optimization: BFGS",
  "Primary MLP objective: squared error with optional L2 weight decay",
  "Primary MLP maximum iterations: 1000",
  "Primary MLP early stopping: not used",
  "Primary MLP frequency/differencing parameters: not applicable",
  "Missing-data treatment: no imputation; monthly grid retained; model-specific complete cases after lag construction",
  "Primary one-step evaluation: observed antecedent discharge; observed historical precipitation at required lags",
  "Recursive sensitivity: model-generated discharge returned to predictor history; precipitation is not recursively predicted",
  "ARMA/ARIMA/SARIMA: NA retained at original dates and handled by ML state-space estimation"
), file.path(OUTPUT_DIR, "run_configuration.txt"))

message("Primary analysis completed; running reviewer-requested sensitivity analyses...")


# ------------------------------------------------------------------------------
# 11. Predictor-matched linear-versus-MLP comparisons
# ------------------------------------------------------------------------------
H1_PREDICTORS_Q  <- c("Q1", "Q2", "Q3")
H1_PREDICTORS_QP <- c("Q1", "Q2", "Q3", "P1", "P2")

fit_predict_matched <- function(period, cal_end, val_start,
                                     information_domain = c("Discharge only", "Discharge + precipitation")) {
  information_domain <- match.arg(information_domain)
  preds <- if (information_domain == "Discharge only") H1_PREDICTORS_Q else H1_PREDICTORS_QP
  cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  val <- lagged_data %>% filter(date >= val_start, date <= END_DATE)

  cv_lm <- select_linear_models(cal, list(Matched = preds), period)
  cv_nn <- select_mlp_models(cal, list(Matched = preds), period,
                          ifelse(information_domain == "Discharge only", "Matched Q", "Matched Q+P"))
  best_nn <- cv_nn[1, ]

  calibration_complete <- cal %>% select(freshwater_discharge_m3_s, all_of(preds)) %>% drop_na()
  fit_l <- lm(reformulate(preds, response = "freshwater_discharge_m3_s"), data = calibration_complete)
  fit_n <- fit_mlp(calibration_complete, preds, best_nn$Hidden, best_nn$Decay,
                   seed = SEED + ifelse(period == "T2", 2000L, 0L))

  linear_prediction <- rep(NA_real_, nrow(val))
  ok_l <- complete.cases(val[, preds, drop = FALSE])
  linear_prediction[ok_l] <- as.numeric(predict(fit_l, newdata = val[ok_l, , drop = FALSE]))
  mlp_prediction <- predict_mlp_complete_cases(fit_n, val)

  prediction_list <- bind_rows(
    tibble(date = val$date, Observed = val$freshwater_discharge_original_m3_s, Predicted = linear_prediction,
           Model_class = "Linear"),
    tibble(date = val$date, Observed = val$freshwater_discharge_original_m3_s, Predicted = mlp_prediction,
           Model_class = "MLP")
  ) %>%
    mutate(Period = period, Information_domain = information_domain,
           Predictors = paste(preds, collapse = "+"))

  list(predictions = prediction_list,
       cv_linear = cv_lm %>% mutate(Information_domain = information_domain),
       cv_mlp = cv_nn %>% mutate(Information_domain = information_domain))
}

matched_scenarios <- list(
  T1_Q = fit_predict_matched("T1", T1_CAL_END, T1_VAL_START, "Discharge only"),
  T1_QP = fit_predict_matched("T1", T1_CAL_END, T1_VAL_START, "Discharge + precipitation"),
  T2_Q = fit_predict_matched("T2", T2_CAL_END, T2_VAL_START, "Discharge only"),
  T2_QP = fit_predict_matched("T2", T2_CAL_END, T2_VAL_START, "Discharge + precipitation")
)

matched_h1_predictions <- map_dfr(matched_scenarios, "predictions")
matched_h1_cv_linear <- map_dfr(matched_scenarios, "cv_linear")
matched_h1_cv_mlp <- map_dfr(matched_scenarios, "cv_mlp")
write_csv(matched_h1_predictions, file.path(OUTPUT_DIR, "matched_predictor_predictions.csv"))
write_csv(matched_h1_cv_linear, file.path(OUTPUT_DIR, "matched_predictor_cv_linear.csv"))
write_csv(matched_h1_cv_mlp, file.path(OUTPUT_DIR, "matched_predictor_cv_mlp.csv"))

calculate_h1_contrast <- function(df) {
  z <- df %>%
    select(date, Observed, Model_class, Predicted) %>%
    pivot_wider(names_from = Model_class, values_from = Predicted) %>%
    filter(if_all(c(Observed, Linear, MLP), is.finite)) %>%
    arrange(date)
  e_l <- z$Observed - z$Linear
  e_n <- z$Observed - z$MLP
  dm <- tryCatch(
    forecast::dm.test(e_l, e_n, alternative = "two.sided", h = 1, power = 1),
    error = function(e) NULL
  )
  delta_abs <- abs(e_n) - abs(e_l)
  ci <- moving_block_mean_ci(delta_abs, B = BOOT_B,
                             block_length = BOOT_BLOCK_LENGTH,
                             seed = SEED)
  tibble(
    N = nrow(z),
    RMSE_linear = sqrt(mean(e_l^2)),
    RMSE_MLP = sqrt(mean(e_n^2)),
    MAE_linear = mean(abs(e_l)),
    MAE_MLP = mean(abs(e_n)),
    Delta_RMSE_MLP_minus_linear = sqrt(mean(e_n^2)) - sqrt(mean(e_l^2)),
    Delta_MAE_MLP_minus_linear = mean(delta_abs),
    Delta_MAE_CI_low = ci[["lower"]],
    Delta_MAE_CI_high = ci[["upper"]],
    DM_statistic = ifelse(is.null(dm), NA_real_, unname(dm$statistic)),
    DM_p_value = ifelse(is.null(dm), NA_real_, dm$p.value)
  )
}

matched_h1_tests <- matched_h1_predictions %>%
  group_by(Period, Information_domain, Predictors) %>%
  group_modify(~calculate_h1_contrast(.x)) %>%
  ungroup() %>%
  mutate(DM_p_adjusted_BH = p.adjust(DM_p_value, method = "BH"))
write_csv(matched_h1_tests, file.path(OUTPUT_DIR, "matched_predictor_error_tests.csv"))

precipitation_increment <- matched_h1_predictions %>%
  select(Period, date, Observed, Model_class, Information_domain, Predicted) %>%
  pivot_wider(names_from = Information_domain, values_from = Predicted) %>%
  filter(if_all(c(Observed, `Discharge only`, `Discharge + precipitation`), is.finite)) %>%
  group_by(Period, Model_class) %>%
  summarise(
    N = n(),
    MAE_Q = mean(abs(Observed - `Discharge only`)),
    MAE_QP = mean(abs(Observed - `Discharge + precipitation`)),
    Delta_MAE_QP_minus_Q = MAE_QP - MAE_Q,
    RMSE_Q = sqrt(mean((Observed - `Discharge only`)^2)),
    RMSE_QP = sqrt(mean((Observed - `Discharge + precipitation`)^2)),
    Delta_RMSE_QP_minus_Q = RMSE_QP - RMSE_Q,
    .groups = "drop"
  )
write_csv(precipitation_increment,
          file.path(OUTPUT_DIR, "incremental_precipitation_value.csv"))


# ------------------------------------------------------------------------------
# 12. Sensitivity to 12- and 24-month lag windows
# ------------------------------------------------------------------------------
lag_sensitivity_sets <- list(
  Q_1_12 = paste0("Q", 1:12),
  Q_1_24 = paste0("Q", 1:24),
  QP_1_12 = c(paste0("Q", 1:12), paste0("P", 1:12))
)

select_mlp_lag_sensitivity <- function(df_cal, period) {
  monthly_grid <- crossing(Model = names(lag_sensitivity_sets),
                    Hidden = c(2L, 4L, 6L, 8L, 10L),
                    Decay = c(0.001, 0.01))
  pmap_dfr(monthly_grid, function(Model, Hidden, Decay) {
    preds <- lag_sensitivity_sets[[Model]]
    bind_cols(
      tibble(Period = period, Model = Model,
             Predictors = paste(preds, collapse = "+"),
             Hidden = Hidden, Decay = Decay,
             Architecture = paste(length(preds), Hidden, 1, sep = "-")),
      cv_mlp_candidate(df_cal, preds, Hidden, Decay, SEED + 5000L)
    )
  }) %>%
    filter(is.finite(CV_RMSE), CV_success >= 0.90) %>%
    group_by(Model) %>%
    arrange(CV_RMSE, CV_MAE, Hidden, Decay, .by_group = TRUE) %>%
    mutate(Rank_within_input = row_number()) %>%
    ungroup()
}

mlp_lag_sensitivity <- bind_rows(
  select_mlp_lag_sensitivity(lagged_data %>% filter(date <= T1_CAL_END), "T1"),
  select_mlp_lag_sensitivity(lagged_data %>% filter(date <= T2_CAL_END), "T2")
)
write_csv(mlp_lag_sensitivity,
          file.path(OUTPUT_DIR, "mlp_lag_window_sensitivity.csv"))
write_csv(mlp_lag_sensitivity %>% filter(Rank_within_input == 1L),
          file.path(OUTPUT_DIR, "mlp_lag_window_best_models.csv"))


# ------------------------------------------------------------------------------
# 13. precipitation_mm_month-discharge climatology, correlations and cross-correlation
# ------------------------------------------------------------------------------
qp_climatology <- monthly_data %>%
  group_by(month) %>%
  summarise(Q_mean = mean(freshwater_discharge_original_m3_s, na.rm = TRUE), Q_sd = sd(freshwater_discharge_original_m3_s, na.rm = TRUE),
            P_mean = mean(precipitation_original_mm_month, na.rm = TRUE), P_sd = sd(precipitation_original_mm_month, na.rm = TRUE),
            N_Q = sum(is.finite(freshwater_discharge_original_m3_s)), N_P = sum(is.finite(precipitation_original_mm_month)), .groups = "drop")
write_csv(qp_climatology, file.path(OUTPUT_DIR, "discharge_precipitation_climatology.csv"))

precipitation_correlations <- map_dfr(c(0L, 1L, 2L, 3L, 12L), function(k) {
  p <- if (k == 0L) lagged_data$precipitation_original_mm_month else lagged_data[[paste0("P", k)]]
  ok <- is.finite(lagged_data$freshwater_discharge_original_m3_s) & is.finite(p)
  test_result <- cor.test(lagged_data$freshwater_discharge_original_m3_s[ok], p[ok], method = "spearman", exact = FALSE)
  tibble(Precipitation_lag = k, N = sum(ok), Spearman_rho = unname(test_result$estimate),
         P_value = test_result$p.value)
}) %>% mutate(P_adjusted_BH = p.adjust(P_value, method = "BH"))
write_csv(precipitation_correlations,
          file.path(OUTPUT_DIR, "discharge_precipitation_correlations.csv"))

qp_ccf <- lagged_data %>% filter(is.finite(freshwater_discharge_original_m3_s), is.finite(precipitation_original_mm_month))
ccf_obj <- ccf(qp_ccf$precipitation_original_mm_month, qp_ccf$freshwater_discharge_original_m3_s, lag.max = 24, plot = FALSE, na.action = na.pass)
ccf_results <- tibble(Lag = as.integer(ccf_obj$lag), CCF = as.numeric(ccf_obj$acf))
write_csv(ccf_results, file.path(OUTPUT_DIR, "precipitation_discharge_ccf.csv"))


# ------------------------------------------------------------------------------
# 14. Flow-dependent validation error structure
# ------------------------------------------------------------------------------
flow_class_errors <- predictions %>%
  transmute(Period, Family, Model, date, Observed, Predicted = OneStep) %>%
  filter(is.finite(Observed), is.finite(Predicted)) %>%
  group_by(Period) %>%
  mutate(Flow_class = ntile(Observed, 3L),
         Flow_class = factor(Flow_class, levels = 1:3,
                             labels = c("Low", "Intermediate", "High"))) %>%
  ungroup() %>%
  group_by(Period, Family, Model, Flow_class) %>%
  summarise(N = n(), RMSE = sqrt(mean((Observed - Predicted)^2)),
            MAE = mean(abs(Observed - Predicted)), Bias = mean(Predicted - Observed),
            Underprediction_fraction = mean(Predicted < Observed), .groups = "drop")
write_csv(flow_class_errors, file.path(OUTPUT_DIR, "flow_class_validation_performance.csv"))

one_step_residuals <- predictions %>%
  transmute(Period, Family, Model, date, Observed, Predicted = OneStep,
            Error = Predicted - Observed, Absolute_error = abs(Predicted - Observed)) %>%
  filter(is.finite(Observed), is.finite(Predicted))
write_csv(one_step_residuals, file.path(OUTPUT_DIR, "one_step_residuals.csv"))


# ------------------------------------------------------------------------------
# 15. T1-versus-T2 comparison on exact common post-2012 prediction dates
# ------------------------------------------------------------------------------
# Every family is paired on the intersection of dates with a finite observed
# target and finite T1 and T2 one-step predictions. Thus T1 and T2 metrics within
# a family are calculated from exactly the same observed months.

h2_paired_predictions <- predictions %>%
  filter(date >= T2_VAL_START) %>%
  select(Period, Family, date, Observed, Predicted = OneStep) %>%
  pivot_wider(names_from = Period, values_from = Predicted) %>%
  filter(if_all(c(Observed, T1, T2), is.finite)) %>%
  arrange(Family, date)

h2_common_date_comparison <- h2_paired_predictions %>%
  group_by(Family) %>%
  group_modify(function(.x, .g) {
    obs <- .x$Observed
    p1 <- .x$T1
    p2 <- .x$T2

    rmse1 <- sqrt(mean((obs - p1)^2))
    rmse2 <- sqrt(mean((obs - p2)^2))
    mae1 <- mean(abs(obs - p1))
    mae2 <- mean(abs(obs - p2))
    observed_sd <- sd(obs)

    scale_t1 <- mase_scale(calibration_observed[["T1"]], 12L)
    scale_t2 <- mase_scale(calibration_observed[["T2"]], 12L)

    tibble(
      N = length(obs),
      First_date = min(.x$date),
      Last_date = max(.x$date),
      RMSE_T1 = rmse1,
      RMSE_T2 = rmse2,
      Delta_RMSE_T2_minus_T1 = rmse2 - rmse1,
      MAE_T1 = mae1,
      MAE_T2 = mae2,
      Delta_MAE_T2_minus_T1 = mae2 - mae1,
      RSR_T1 = rmse1 / observed_sd,
      RSR_T2 = rmse2 / observed_sd,
      MASE_T1 = mae1 / scale_t1,
      MASE_T2 = mae2 / scale_t2,
      Relative_variance_error_T1 =
        (var(p1) - var(obs)) / var(obs),
      Relative_variance_error_T2 =
        (var(p2) - var(obs)) / var(obs)
    )
  }) %>%
  ungroup()

write_csv(
  h2_paired_predictions,
  file.path(OUTPUT_DIR, "common_post2012_exact_paired_predictions.csv")
)
write_csv(
  h2_common_date_comparison,
  file.path(OUTPUT_DIR, "common_post2012_exact_date_metrics.csv")
)

h2_paired_error_tests <- h2_paired_predictions %>%
  group_by(Family) %>%
  group_modify(function(.x, .g) {
    e_t1 <- .x$Observed - .x$T1
    e_t2 <- .x$Observed - .x$T2

    dm <- tryCatch(
      forecast::dm.test(
        e_t1,
        e_t2,
        alternative = "two.sided",
        h = 1,
        power = 1
      ),
      error = function(e) NULL
    )

    delta_abs <- abs(e_t2) - abs(e_t1)
    ci <- moving_block_mean_ci(
      delta_abs,
      B = BOOT_B,
      block_length = BOOT_BLOCK_LENGTH,
      seed = SEED + 17L
    )

    tibble(
      N = nrow(.x),
      Delta_MAE_T2_minus_T1 = mean(delta_abs),
      Delta_MAE_CI_low = ci[["lower"]],
      Delta_MAE_CI_high = ci[["upper"]],
      DM_statistic = ifelse(is.null(dm), NA_real_, unname(dm$statistic)),
      DM_p_value = ifelse(is.null(dm), NA_real_, dm$p.value)
    )
  }) %>%
  ungroup() %>%
  mutate(DM_p_adjusted_BH = p.adjust(DM_p_value, method = "BH"))

write_csv(
  h2_paired_error_tests,
  file.path(OUTPUT_DIR, "common_post2012_paired_error_tests.csv")
)


# ------------------------------------------------------------------------------
# 16. Direct pre-2012 versus post-2012 hydroclimatic comparisons
# ------------------------------------------------------------------------------
h3 <- monthly_data %>%
  filter(is.finite(freshwater_discharge_original_m3_s)) %>%
  mutate(Hydroperiod = factor(ifelse(date < T2_VAL_START, "Pre-2012", "Post-2012"),
                              levels = c("Pre-2012", "Post-2012")))

h3_descriptive <- h3 %>% group_by(Hydroperiod) %>%
  summarise(N = n(), Mean = mean(freshwater_discharge_original_m3_s), Median = median(freshwater_discharge_original_m3_s),
            SD = sd(freshwater_discharge_original_m3_s), Variance = var(freshwater_discharge_original_m3_s),
            CV = SD / Mean, Q05 = quantile(freshwater_discharge_original_m3_s, 0.05),
            Q95 = quantile(freshwater_discharge_original_m3_s, 0.95), .groups = "drop")

mean_test <- t.test(freshwater_discharge_original_m3_s ~ Hydroperiod, data = h3, var.equal = FALSE)
location_test <- wilcox.test(freshwater_discharge_original_m3_s ~ Hydroperiod, data = h3, exact = FALSE)
median_deviations <- h3 %>% group_by(Hydroperiod) %>%
  mutate(Abs_dev_median = abs(freshwater_discharge_original_m3_s - median(freshwater_discharge_original_m3_s))) %>% ungroup()
bf_fit <- lm(Abs_dev_median ~ Hydroperiod, data = median_deviations)
bf_tab <- anova(bf_fit)
seasonal_fit <- lm(freshwater_discharge_original_m3_s ~ factor(month) * Hydroperiod, data = h3)
seasonal_table <- drop1(seasonal_fit, test = "F")
interaction_row <- grep(":", rownames(seasonal_table), fixed = TRUE)
if (length(interaction_row) != 1L) {
  stop("The month-by-period interaction term could not be identified in the seasonal comparison.")
}

h3_direct_tests <- tibble(
  Test = c("Welch mean comparison", "Wilcoxon location comparison",
           "Brown-Forsythe variance comparison", "month-by-period interaction"),
  Statistic = c(unname(mean_test$statistic), unname(location_test$statistic),
                unname(bf_tab$`F value`[1]),
                unname(seasonal_table$`F value`[interaction_row])),
  P_value = c(mean_test$p.value, location_test$p.value,
              bf_tab$`Pr(>F)`[1], seasonal_table$`Pr(>F)`[interaction_row])
) %>% mutate(P_adjusted_BH = p.adjust(P_value, method = "BH"))

acf_at_lag <- function(x, lag) as.numeric(acf(x, lag.max = lag, plot = FALSE,
                                               na.action = na.pass)$acf[lag + 1L])
h3_acf <- h3 %>% group_by(Hydroperiod) %>%
  summarise(ACF_lag1 = acf_at_lag(freshwater_discharge_original_m3_s, 1L),
            ACF_lag12 = acf_at_lag(freshwater_discharge_original_m3_s, 12L), .groups = "drop")

bootstrap_acf_period <- function(x, period, B = BOOT_B, block_length = BOOT_BLOCK_LENGTH) {
  x <- x[is.finite(x)]
  sims <- replicate(B, {
    xb <- sample_moving_blocks(x, length(x), block_length)
    c(lag1 = acf_at_lag(xb, 1L), lag12 = acf_at_lag(xb, 12L))
  })
  tibble(
    Hydroperiod = period,
    Lag = c(1L, 12L),
    Estimate = c(acf_at_lag(x, 1L), acf_at_lag(x, 12L)),
    Lower95 = apply(sims, 1L, quantile, 0.025, na.rm = TRUE),
    Upper95 = apply(sims, 1L, quantile, 0.975, na.rm = TRUE),
    B = B,
    Block_length = block_length
  )
}

h3_acf_bootstrap <- bind_rows(
  bootstrap_acf_period(h3$freshwater_discharge_original_m3_s[h3$Hydroperiod == "Pre-2012"], "Pre-2012"),
  bootstrap_acf_period(h3$freshwater_discharge_original_m3_s[h3$Hydroperiod == "Post-2012"], "Post-2012")
)

h3_climatology <- h3 %>%
  group_by(Hydroperiod, month) %>%
  summarise(
    N = n(), Mean = mean(freshwater_discharge_original_m3_s), SD = sd(freshwater_discharge_original_m3_s),
    SE = SD / sqrt(N), Lower95 = Mean - qt(0.975, pmax(N - 1L, 1L)) * SE,
    Upper95 = Mean + qt(0.975, pmax(N - 1L, 1L)) * SE,
    .groups = "drop"
  )

write_csv(h3_descriptive, file.path(OUTPUT_DIR, "pre_post2012_descriptive_statistics.csv"))
write_csv(h3_direct_tests, file.path(OUTPUT_DIR, "pre_post2012_direct_tests.csv"))
write_csv(h3_acf, file.path(OUTPUT_DIR, "pre_post2012_acf.csv"))
write_csv(h3_acf_bootstrap,
          file.path(OUTPUT_DIR, "pre_post2012_acf_bootstrap.csv"))
write_csv(h3_climatology,
          file.path(OUTPUT_DIR, "pre_post2012_monthly_climatology.csv"))


cv_benchmark_summary <- function(cal_end, period) {
  df <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  idx <- rolling_origins(df)
  target_row <- idx + 1L
  ok <- df$freshwater_discharge_observed[target_row] & is.finite(df$freshwater_discharge_m3_s[target_row]) & is.finite(df$Q12[target_row])
  error_value <- df$freshwater_discharge_m3_s[target_row][ok] - df$Q12[target_row][ok]
  tibble(Period = period, Family = "Benchmark", Model = "Seasonal naive",
         Predictors = "Q12", CV_RMSE = sqrt(mean(error_value^2)),
         CV_MAE = mean(abs(error_value)), CV_N = length(error_value),
         CV_total = sum(ok), CV_success = 1)
}

extract_selected_cv <- function(selections, period) {
  ar <- map_dfr(c("ARMA", "ARIMA", "SARIMA"), function(family) {
    x <- selections$ar[[family]][1, ]
    tibble(
      Period = period, Family = family,
      Model = sprintf("%s(%d,%d,%d)(%d,%d,%d)[12]", family,
                      x$p, x$d, x$q, x$P, x$D, x$Q),
      Predictors = "monthly discharge series",
      CV_RMSE = x$CV_RMSE, CV_MAE = x$CV_MAE, CV_N = x$CV_N,
      CV_total = x$CV_total, CV_success = x$CV_success
    )
  })
  reg <- bind_rows(
    selections$lm[1, ] %>% transmute(
      Period = period, Family = "Linear regression",
      Model = paste0("Linear [", Model, "]"), Predictors,
      CV_RMSE, CV_MAE, CV_N, CV_total, CV_success),
    selections$mlp_uni[1, ] %>% transmute(
      Period = period, Family = "MLP univariate",
      Model = paste0("MLP univariate [", Model, "; ", Architecture, "]"),
      Predictors, CV_RMSE, CV_MAE, CV_N, CV_total, CV_success),
    selections$mlp_multi[1, ] %>% transmute(
      Period = period, Family = "MLP multivariate",
      Model = paste0("MLP multivariate [", Model, "; ", Architecture, "]"),
      Predictors, CV_RMSE, CV_MAE, CV_N, CV_total, CV_success)
  )
  bind_rows(ar, reg)
}

selected_cv_all_families <- bind_rows(
  cv_benchmark_summary(T1_CAL_END, "T1"),
  extract_selected_cv(selected_T1, "T1"),
  cv_benchmark_summary(T2_CAL_END, "T2"),
  extract_selected_cv(selected_T2, "T2")
) %>% arrange(Period, CV_RMSE, CV_MAE)
write_csv(selected_cv_all_families,
          file.path(OUTPUT_DIR, "rolling_origin_selected_models_summary.csv"))

cv_validation_comparison <- selected_cv_all_families %>%
  select(Period, Family, Model, CV_N, CV_RMSE, CV_MAE) %>%
  inner_join(
    performance_long %>% filter(Evaluation == "OneStep") %>%
      select(Period, Family, Model, Validation_N = N,
             Validation_RMSE = RMSE, Validation_MAE = MAE),
    by = c("Period", "Family", "Model")
  ) %>%
  mutate(RMSE_validation_minus_CV = Validation_RMSE - CV_RMSE,
         MAE_validation_minus_CV = Validation_MAE - CV_MAE)
write_csv(cv_validation_comparison,
          file.path(OUTPUT_DIR, "rolling_cv_vs_validation_performance.csv"))

duplicate_predictions <- predictions %>% count(Period, Family, Model, date) %>% filter(n != 1L)
if (nrow(duplicate_predictions) > 0L) stop("Duplicate predictions were detected. See duplicate_prediction_audit.csv")
if (any(!is.finite(performance_long$RMSE[performance_long$N >= 2L]))) {
  stop("A non-finite metric was found for a model with at least two valid observed-predicted pairs.")
}
write_csv(duplicate_predictions, file.path(OUTPUT_DIR, "duplicate_prediction_audit.csv"))

final_audit <- predictions %>%
  group_by(Period, Family, Model) %>%
  summarise(N_calendar = n(), N_observed = sum(is.finite(Observed)),
            N_one_step = sum(is.finite(OneStep)),
            N_valid_pairs = sum(is.finite(Observed) & is.finite(OneStep)),
            First_date = min(date), Last_date = max(date), .groups = "drop")
write_csv(final_audit, file.path(OUTPUT_DIR, "final_prediction_audit.csv"))

fig_h1 <- matched_h1_predictions %>%
  group_by(Period, Information_domain, Model_class) %>%
  group_modify(~compute_metrics(.x$Observed, .x$Predicted,
                         calibration_observed[[as.character(.y$Period)]])) %>%
  ungroup() %>%
  select(Period, Information_domain, Model_class, RMSE, MAE) %>%
  pivot_longer(c(RMSE, MAE), names_to = "Metric", values_to = "Error") %>%
  ggplot(aes(Model_class, Error, colour = Model_class)) +
  geom_point(size = 2.5) +
  facet_grid(Metric ~ Period + Information_domain, scales = "free_y") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = expression(Error~(m^3~s^-1)))
ggsave(file.path(OUTPUT_DIR, "matched_predictor_comparison.png"), fig_h1,
       width = 12, height = 6, dpi = 600)

fig_precip_climatology <- qp_climatology %>%
  select(month, freshwater_discharge_m3_s = Q_mean, precipitation_mm_month = P_mean) %>%
  pivot_longer(-month, names_to = "Variable", values_to = "Mean") %>%
  ggplot(aes(month, Mean, colour = Variable)) + geom_line(linewidth = 0.8) +
  geom_point() + facet_wrap(~Variable, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = 1:12) + theme_bw() +
  labs(x = "month", y = "Monthly mean", colour = NULL)
ggsave(file.path(OUTPUT_DIR, "precipitation_discharge_climatology_plot.png"),
       fig_precip_climatology, width = 7.5, height = 6, dpi = 600)

fig_flow_error <- one_step_residuals %>%
  ggplot(aes(Observed, Absolute_error)) + geom_point(alpha = 0.45, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, colour = "#B2182B") +
  facet_grid(Family ~ Period, scales = "free_y") + theme_bw() +
  labs(x = expression(Observed~discharge~(m^3~s^-1)),
       y = expression(Absolute~error~(m^3~s^-1)))
ggsave(file.path(OUTPUT_DIR, "error_vs_observed_flow.png"), fig_flow_error,
       width = 8, height = 12, dpi = 600)

fig_cv_all <- selected_cv_all_families %>%
  select(Period, Family, CV_RMSE, CV_MAE) %>%
  pivot_longer(c(CV_RMSE, CV_MAE), names_to = "Metric", values_to = "Error") %>%
  mutate(Metric = recode(Metric, CV_RMSE = "RMSE", CV_MAE = "MAE")) %>%
  ggplot(aes(Family, Error, colour = Period, shape = Period)) +
  geom_point(size = 2.7, position = position_dodge(width = 0.45)) +
  facet_wrap(~Metric, scales = "free_y") + theme_bw() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(x = NULL, y = expression(Rolling-origin~error~(m^3~s^-1)),
       colour = "Configuration", shape = "Configuration")
ggsave(file.path(OUTPUT_DIR, "rolling_origin_all_families.png"),
       fig_cv_all, width = 10, height = 5.5, dpi = 600)

h2_plot_data <- h2_common_date_comparison %>%
  select(
    Family,
    RMSE_T1, RMSE_T2,
    MAE_T1, MAE_T2,
    RSR_T1, RSR_T2,
    MASE_T1, MASE_T2,
    Relative_variance_error_T1,
    Relative_variance_error_T2
  ) %>%
  pivot_longer(
    -Family,
    names_to = "Metric_period",
    values_to = "Value"
  ) %>%
  mutate(
    Period = ifelse(grepl("_T1$", Metric_period), "T1", "T2"),
    Metric = sub("_T[12]$", "", Metric_period),
    Metric = recode(
      Metric,
      RMSE = "RMSE",
      MAE = "MAE",
      RSR = "RSR",
      MASE = "MASE",
      Relative_variance_error = "Relative variance error"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "RMSE",
        "MAE",
        "RSR",
        "MASE",
        "Relative variance error"
      )
    ),
    Period = factor(Period, levels = c("T1", "T2"))
  )

fig_h2 <- ggplot(
  h2_plot_data,
  aes(Period, Value, group = Family, colour = Family)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    x = "Temporal configuration",
    y = NULL,
    colour = "Model family"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(OUTPUT_DIR, "common_post2012_comparison.png"),
  fig_h2,
  width = 10,
  height = 7,
  dpi = 600
)

fig_h3 <- h3_climatology %>%
  ggplot(aes(month, Mean, colour = Hydroperiod, group = Hydroperiod)) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.12,
                position = position_dodge(width = 0.12), alpha = 0.55) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) + theme_bw() +
  labs(x = NULL, y = expression(Mean~monthly~discharge~(m^3~s^-1)),
       colour = "Hydroclimatic period")
ggsave(file.path(OUTPUT_DIR, "pre_post2012_monthly_climatology_plot.png"),
       fig_h3, width = 9, height = 5.5, dpi = 600)

fig_qp_scatter <- lagged_data %>%
  filter(is.finite(freshwater_discharge_original_m3_s), is.finite(precipitation_original_mm_month)) %>%
  ggplot(aes(precipitation_original_mm_month, freshwater_discharge_original_m3_s)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", se = TRUE, colour = "#2166AC") + theme_bw() +
  labs(x = "Monthly precipitation (mm)",
       y = expression(Monthly~discharge~(m^3~s^-1)))
ggsave(file.path(OUTPUT_DIR, "precipitation_discharge_scatter.png"),
       fig_qp_scatter, width = 7, height = 5.5, dpi = 600)

ccf_limit <- 1.96 / sqrt(nrow(qp_ccf))
fig_ccf <- ccf_results %>%
  ggplot(aes(Lag, CCF)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_hline(yintercept = c(-ccf_limit, ccf_limit), linetype = 2,
             colour = "#B2182B") +
  geom_segment(aes(xend = Lag, y = 0, yend = CCF), linewidth = 0.6,
               colour = "#2166AC") + theme_bw() +
  scale_x_continuous(breaks = seq(-24, 24, 4)) +
  labs(x = "Lag (months; negative lags indicate precipitation leading discharge)",
       y = "CCF")
ggsave(file.path(OUTPUT_DIR, "precipitation_discharge_ccf_plot.png"),
       fig_ccf, width = 8, height = 5, dpi = 600)

capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "session_info_final.txt"))
writeLines(c(
  paste("Execution ID:", RUN_ID), paste("Seed:", SEED),
  paste("MLP CV initializations per fold:", MLP_CV_RESTARTS),
  paste("MLP final initializations:", MLP_FINAL_RESTARTS),
  paste("Bootstrap replicates:", BOOT_B), paste("Block length:", BOOT_BLOCK_LENGTH),
  "No missing value was imputed.",
  "Primary performance assessment: held-out one-step predictions.",
  "Fully recursive and multi-horizon predictions are supplementary.",
  "Multi-horizon precipitation-informed predictions are conditional on the observed historical precipitation inputs available in the withheld record.",
  "H1 primary contrast: Linear versus MLP with exactly matched predictors.",
  "T1 versus T2 comparison: common post-2012 dates and one-step horizon.",
  "MLPs were fitted with repeated random initializations and weight decay; no separate early-stopping validation subset was used."
), file.path(OUTPUT_DIR, "run_notes.txt"))

message("Core analysis completed. Outputs written to: ", normalizePath(OUTPUT_DIR))

objects_to_save <- list(
  monthly_data = monthly_data,
  lagged_data = lagged_data,
  predictions_T1 = predictions_T1,
  predictions_T2 = predictions_T2,
  predictions = predictions,
  selected_T1 = selected_T1,
  selected_T2 = selected_T2,
  selected_predictor_summary = selected_predictor_summary
)

saveRDS(
  objects_to_save,
  file.path(OUTPUT_DIR, "analysis_objects.rds")
)


# ------------------------------------------------------------------------------
# 17. Calibration-only moving-block residual bootstrap for MLP uncertainty
# ------------------------------------------------------------------------------
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

BOOT_B <- 2000L
BOOT_BLOCK_LENGTH <- 12L
BOOT_SEED <- 20260826L

BOOT_INITIAL_WINDOW <- 120L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


bootstrap_model_specs <- bind_rows(
  
  selected_T1$mlp_uni |>
    slice(1) |>
    transmute(
      Period = "T1",
      Family = "MLP univariate",
      Calibration_start = T1_CAL_START,
      Calibration_end = T1_CAL_END,
      Predictors,
      Hidden = as.integer(Hidden),
      Decay = as.numeric(Decay),
      Model,
      Architecture
    ),
  
  selected_T1$mlp_multi |>
    slice(1) |>
    transmute(
      Period = "T1",
      Family = "MLP multivariate",
      Calibration_start = T1_CAL_START,
      Calibration_end = T1_CAL_END,
      Predictors,
      Hidden = as.integer(Hidden),
      Decay = as.numeric(Decay),
      Model,
      Architecture
    ),
  
  selected_T2$mlp_uni |>
    slice(1) |>
    transmute(
      Period = "T2",
      Family = "MLP univariate",
      Calibration_start = T1_CAL_START,
      Calibration_end = T2_CAL_END,
      Predictors,
      Hidden = as.integer(Hidden),
      Decay = as.numeric(Decay),
      Model,
      Architecture
    ),
  
  selected_T2$mlp_multi |>
    slice(1) |>
    transmute(
      Period = "T2",
      Family = "MLP multivariate",
      Calibration_start = T1_CAL_START,
      Calibration_end = T2_CAL_END,
      Predictors,
      Hidden = as.integer(Hidden),
      Decay = as.numeric(Decay),
      Model,
      Architecture
    )
)

print(bootstrap_model_specs, n = Inf)


generate_calibration_residuals_mlp <- function(specification) {
  
  period <- specification$Period[[1]]
  family <- specification$Family[[1]]
  
  calibration_start <- specification$Calibration_start[[1]]
  calibration_end <- specification$Calibration_end[[1]]
  
  predictors <- strsplit(
    specification$Predictors[[1]],
    split = "\\+"
  )[[1]]
  
  predictors <- trimws(predictors)
  
  hidden <- specification$Hidden[[1]]
  decay <- specification$Decay[[1]]
  model_name <- specification$Model[[1]]
  architecture <- specification$Architecture[[1]]
  
  calibration_data <- lagged_data |>
    filter(
      date >= calibration_start,
      date <= calibration_end
    ) |>
    arrange(date)
  
  if (nrow(calibration_data) <= BOOT_INITIAL_WINDOW) {
    stop(
      "Calibration period is shorter than BOOT_INITIAL_WINDOW for ",
      period, " - ", family
    )
  }
  
  target_indices <- seq.int(
    from = BOOT_INITIAL_WINDOW + 1L,
    to = nrow(calibration_data),
    by = 1L
  )
  
  seed_offset <- case_when(
    period == "T1" & family == "MLP univariate"   ~ 10000L,
    period == "T1" & family == "MLP multivariate" ~ 20000L,
    period == "T2" & family == "MLP univariate"   ~ 30000L,
    TRUE                                             ~ 40000L
  )
  
  map_dfr(target_indices, function(i) {
    
    target_row <- calibration_data[i, , drop = FALSE]
    
    training <- calibration_data[seq_len(i - 1L), , drop = FALSE] |>
      select(freshwater_discharge_m3_s, all_of(predictors)) |>
      drop_na()
    
    observed_valid <- (
      isTRUE(target_row$freshwater_discharge_observed[[1]]) &&
        is.finite(target_row$freshwater_discharge_original_m3_s[[1]])
    )
    
    predictor_values <- unlist(
      target_row[, predictors, drop = FALSE],
      use.names = FALSE
    )
    
    predictors_valid <- (
      length(predictor_values) == length(predictors) &&
        all(is.finite(predictor_values))
    )
    
    if (
      nrow(training) < 60L ||
      !observed_valid ||
      !predictors_valid
    ) {
      return(
        tibble(
          Period = period,
          Family = family,
          Model = model_name,
          Architecture = architecture,
          date = target_row$date[[1]],
          Observed = target_row$freshwater_discharge_original_m3_s[[1]],
          Predicted = NA_real_,
          Residual = NA_real_,
          Success = FALSE
        )
      )
    }
    
    fitted_model <- tryCatch(
      fit_mlp(
        train = training,
        preds = predictors,
        hidden = hidden,
        decay = decay,
        seed = BOOT_SEED + seed_offset + i
      ),
      error = function(e) NULL
    )
    
    if (is.null(fitted_model)) {
      predicted_value <- NA_real_
    } else {
      predicted_value <- tryCatch(
        as.numeric(
          predict_mlp(
            fitted_model,
            target_row[, predictors, drop = FALSE]
          )
        ),
        error = function(e) NA_real_
      )
    }
    
    success_flag <- is.finite(predicted_value)
    
    tibble(
      Period = period,
      Family = family,
      Model = model_name,
      Architecture = architecture,
      date = target_row$date[[1]],
      Observed = target_row$freshwater_discharge_original_m3_s[[1]],
      Predicted = predicted_value,
      Residual = ifelse(
        success_flag,
        target_row$freshwater_discharge_original_m3_s[[1]] - predicted_value,
        NA_real_
      ),
      Success = success_flag
    )
  })
}

set.seed(BOOT_SEED)

mlp_calibration_residuals <- map_dfr(
  seq_len(nrow(bootstrap_model_specs)),
  function(i) {
    message(
      "Generating calibration residuals: ",
      bootstrap_model_specs$Period[i],
      " - ",
      bootstrap_model_specs$Family[i]
    )
    
    generate_calibration_residuals_mlp(
      bootstrap_model_specs[i, , drop = FALSE]
    )
  }
)


bootstrap_residual_audit <- mlp_calibration_residuals |>
  group_by(Period, Family, Model, Architecture) |>
  summarise(
    N_targets = n(),
    N_success = sum(Success),
    Success_rate = mean(Success),
    Mean_residual = mean(Residual, na.rm = TRUE),
    SD_residual = sd(Residual, na.rm = TRUE),
    RMSE_internal = sqrt(mean(Residual^2, na.rm = TRUE)),
    MAE_internal = mean(abs(Residual), na.rm = TRUE),
    First_target = min(date[Success], na.rm = TRUE),
    Last_target = max(date[Success], na.rm = TRUE),
    .groups = "drop"
  )

print(bootstrap_residual_audit, n = Inf)

if (any(bootstrap_residual_audit$Success_rate < 0.80)) {
  warning(
    "At least one MLP configuration had less than 80% successful ",
    "internal one-step predictions."
  )
}

leakage_check <- mlp_calibration_residuals |>
  left_join(
    bootstrap_model_specs |>
      select(Period, Family, Calibration_end),
    by = c("Period", "Family")
  ) |>
  summarise(
    N_residuals_after_calibration = sum(
      date > Calibration_end,
      na.rm = TRUE
    )
  )

stopifnot(
  leakage_check$N_residuals_after_calibration == 0L
)

write_csv(
  mlp_calibration_residuals,
  file.path(
    OUTPUT_DIR,
    "mlp_calibration_rolling_origin_residuals.csv"
  )
)

write_csv(
  bootstrap_residual_audit,
  file.path(
    OUTPUT_DIR,
    "mlp_bootstrap_residual_audit.csv"
  )
)


extract_contiguous_blocks <- function(
    residual_data,
    block_length = BOOT_BLOCK_LENGTH
) {
  
  z <- residual_data |>
    filter(
      Success,
      is.finite(Residual)
    ) |>
    arrange(date) |>
    mutate(
      month_index = year(date) * 12L + month(date)
    )
  
  if (nrow(z) < block_length) {
    stop("Insufficient residuals to construct moving blocks.")
  }
  
  z <- z |>
    mutate(
      Residual_centered = Residual - mean(Residual, na.rm = TRUE),
      sequence_id = cumsum(
        c(
          TRUE,
          diff(month_index) != 1L
        )
      )
    )
  
  sequences <- split(
    z$Residual_centered,
    z$sequence_id
  )
  
  sequences <- sequences[
    lengths(sequences) >= block_length
  ]
  
  if (length(sequences) == 0L) {
    stop(
      "No sequence contains ",
      block_length,
      " consecutive residuals."
    )
  }
  
  block_list <- map(
    sequences,
    function(v) {
      
      starts <- seq_len(
        length(v) - block_length + 1L
      )
      
      t(
        vapply(
          starts,
          function(j) {
            v[j:(j + block_length - 1L)]
          },
          numeric(block_length)
        )
      )
    }
  )
  
  do.call(rbind, block_list)
}

sample_block_series <- function(
    block_matrix,
    n
) {
  
  block_length <- ncol(block_matrix)
  
  n_blocks <- ceiling(
    n / block_length
  )
  
  sampled_blocks <- sample(
    seq_len(nrow(block_matrix)),
    size = n_blocks,
    replace = TRUE
  )
  
  series_values <- as.vector(
    t(
      block_matrix[
        sampled_blocks,
        ,
        drop = FALSE
      ]
    )
  )
  
  series_values[seq_len(n)]
}

bootstrap_block_audit <- mlp_calibration_residuals |>
  group_by(Period, Family) |>
  group_modify(
    function(.x, .g) {
      blocks <- extract_contiguous_blocks(
        .x,
        block_length = BOOT_BLOCK_LENGTH
      )
      
      tibble(
        Block_length = BOOT_BLOCK_LENGTH,
        N_available_blocks = nrow(blocks),
        N_residuals = sum(
          .x$Success & is.finite(.x$Residual)
        )
      )
    }
  ) |>
  ungroup()

print(bootstrap_block_audit, n = Inf)

write_csv(
  bootstrap_block_audit,
  file.path(
    OUTPUT_DIR,
    "mlp_bootstrap_block_audit.csv"
  )
)


mlp_validation_predictions <- predictions |>
  filter(
    Family %in% c(
      "MLP univariate",
      "MLP multivariate"
    )
  ) |>
  transmute(
    Period,
    Family,
    Model,
    date,
    Observed,
    Predicted = OneStep
  ) |>
  arrange(Period, Family, date)


set.seed(BOOT_SEED)

mlp_predictive_intervals <- mlp_validation_predictions |>
  group_by(Period, Family, Model) |>
  group_modify(
    function(.x, .g) {
      
      current_period <- as.character(.g$Period[[1]])
      current_family <- as.character(.g$Family[[1]])
      
      residual_pool <- mlp_calibration_residuals |>
        filter(
          Period == current_period,
          Family == current_family
        )
      
      block_matrix <- extract_contiguous_blocks(
        residual_pool,
        block_length = BOOT_BLOCK_LENGTH
      )
      
      valid_indices <- which(
        is.finite(.x$Predicted)
      )
      
      lower <- rep(NA_real_, nrow(.x))
      median_boot <- rep(NA_real_, nrow(.x))
      upper <- rep(NA_real_, nrow(.x))
      
      if (length(valid_indices) > 0L) {
        
        valid_predictions <- .x$Predicted[valid_indices]
        n_valid <- length(valid_indices)
        
        simulations <- replicate(
          BOOT_B,
          {
            innovations <- sample_block_series(
              block_matrix,
              n = n_valid
            )
            
            pmax(
              0,
              valid_predictions + innovations
            )
          }
        )
        
        if (is.null(dim(simulations))) {
          simulations <- matrix(
            simulations,
            ncol = 1L
          )
        }
        
        lower[valid_indices] <- apply(
          simulations,
          1,
          quantile,
          probs = 0.025,
          names = FALSE,
          type = 8
        )
        
        median_boot[valid_indices] <- apply(
          simulations,
          1,
          quantile,
          probs = 0.50,
          names = FALSE,
          type = 8
        )
        
        upper[valid_indices] <- apply(
          simulations,
          1,
          quantile,
          probs = 0.975,
          names = FALSE,
          type = 8
        )
      }
      
      .x |>
        mutate(
          Lower95 = lower,
          Bootstrap_median = median_boot,
          Upper95 = upper,
          Bootstrap_B = BOOT_B,
          Block_length = BOOT_BLOCK_LENGTH,
          Interval_type =
            "Calibration rolling-origin moving-block residual predictive interval"
        )
    }
  ) |>
  ungroup()

write_csv(
  mlp_predictive_intervals,
  file.path(
    OUTPUT_DIR,
    "mlp_predictive_intervals.csv"
  )
)


mlp_interval_coverage <- mlp_predictive_intervals |>
  filter(
    is.finite(Observed),
    is.finite(Predicted),
    is.finite(Lower95),
    is.finite(Upper95)
  ) |>
  group_by(Period, Family, Model) |>
  summarise(
    N = n(),
    Coverage95 = mean(
      Observed >= Lower95 &
        Observed <= Upper95
    ),
    Mean_width = mean(Upper95 - Lower95),
    Median_width = median(Upper95 - Lower95),
    Below_interval = mean(Observed < Lower95),
    Above_interval = mean(Observed > Upper95),
    .groups = "drop"
  )

print(mlp_interval_coverage, n = Inf)

write_csv(
  mlp_interval_coverage,
  file.path(
    OUTPUT_DIR,
    "mlp_predictive_interval_coverage.csv"
  )
)


plot_mlp_uncertainty <- function(
    interval_data,
    period
) {
  
  plot_data <- interval_data |>
    filter(Period == period) |>
    mutate(
      Family = factor(
        Family,
        levels = c(
          "MLP univariate",
          "MLP multivariate"
        )
      )
    )
  
  ggplot(
    plot_data,
    aes(x = date)
  ) +
    geom_ribbon(
      aes(
        ymin = Lower95,
        ymax = Upper95
      ),
      fill = "grey70",
      alpha = 0.55,
      na.rm = TRUE
    ) +
    geom_line(
      aes(
        y = Observed,
        colour = "Observed"
      ),
      linewidth = 0.85,
      na.rm = TRUE
    ) +
    geom_line(
      aes(
        y = Predicted,
        colour = "Predicted"
      ),
      linewidth = 0.80,
      na.rm = TRUE
    ) +
    facet_wrap(
      vars(Family),
      ncol = 1,
      scales = "fixed"
    ) +
    scale_colour_manual(
      values = c(
        "Observed" = "black",
        "Predicted" = "#2C7FB8"
      ),
      name = NULL
    ) +
    scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y"
    ) +
    labs(
      x = NULL,
      y = expression(
        freshwater_discharge_m3_s~(m^3~s^{-1})
      )
    ) +
    coord_cartesian(
      ylim = c(0, NA)
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      panel.grid.minor = element_blank()
    )
}


bootstrap_plot_T2 <- plot_mlp_uncertainty(
  mlp_predictive_intervals,
  period = "T2"
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "mlp_predictive_intervals_T2.png"
  ),
  plot = bootstrap_plot_T2,
  width = 12,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "mlp_predictive_intervals_T2.pdf"
  ),
  plot = bootstrap_plot_T2,
  width = 12,
  height = 9,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

print(bootstrap_plot_T2)


bootstrap_plot_T1 <- plot_mlp_uncertainty(
  mlp_predictive_intervals,
  period = "T1"
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "mlp_predictive_intervals_T1.png"
  ),
  plot = bootstrap_plot_T1,
  width = 12,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "mlp_predictive_intervals_T1.pdf"
  ),
  plot = bootstrap_plot_T1,
  width = 12,
  height = 9,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

print(bootstrap_plot_T1)


stopifnot(
  all(
    mlp_predictive_intervals$Lower95[
      is.finite(mlp_predictive_intervals$Lower95)
    ] >= 0
  )
)

stopifnot(
  all(
    mlp_predictive_intervals$Lower95[
      is.finite(mlp_predictive_intervals$Lower95)
    ] <=
      mlp_predictive_intervals$Upper95[
        is.finite(mlp_predictive_intervals$Upper95)
      ]
  )
)

message("Bootstrap analysis and uncertainty figures completed.")


# ------------------------------------------------------------------------------
# 18. Flow-dependent MLP errors in pseudo-out-of-sample calibration and validation
# ------------------------------------------------------------------------------
assign_flow_tertiles <- function(df) {
  class_key <- df %>%
    filter(is.finite(Observed)) %>%
    distinct(Period, date, Observed) %>%
    group_by(Period) %>%
    arrange(Observed, date, .by_group = TRUE) %>%
    mutate(Flow_class = factor(ntile(row_number(), 3L), levels = 1:3,
                               labels = c("Low", "Intermediate", "High"))) %>%
    ungroup() %>%
    select(Period, date, Flow_class)
  left_join(df, class_key, by = c("Period", "date"))
}

summarize_flow_classes <- function(df, stage_label) {
  assign_flow_tertiles(df) %>%
    mutate(Error = Predicted - Observed) %>%
    filter(is.finite(Error), !is.na(Flow_class)) %>%
    group_by(Period, Family, Model, Architecture, Flow_class) %>%
    summarise(
      N = n(),
      RMSE = sqrt(mean(Error^2)),
      MAE = mean(abs(Error)),
      Bias = mean(Error),
      Underprediction_fraction = mean(Error < 0),
      .groups = "drop"
    ) %>%
    mutate(Stage = stage_label)
}

calibration_flow_input <- mlp_calibration_residuals %>%
  filter(Success, is.finite(Observed), is.finite(Predicted)) %>%
  select(Period, Family, Model, Architecture, date, Observed, Predicted)

architecture_key <- bootstrap_model_specs %>%
  select(Period, Family, Architecture) %>% distinct()

validation_flow_input <- predictions %>%
  filter(Family %in% c("MLP univariate", "MLP multivariate")) %>%
  transmute(Period, Family, Model, date, Observed, Predicted = OneStep) %>%
  left_join(architecture_key, by = c("Period", "Family")) %>%
  filter(is.finite(Observed), is.finite(Predicted))

mlp_flow_class_performance <- bind_rows(
  summarize_flow_classes(calibration_flow_input, "Calibration pseudo-out-of-sample"),
  summarize_flow_classes(validation_flow_input, "Independent validation")
) %>%
  arrange(Period, Family, Stage, Flow_class)

write_csv(mlp_flow_class_performance,
          file.path(OUTPUT_DIR, "mlp_flow_class_calibration_validation.csv"))


# ------------------------------------------------------------------------------
# 19. Exploratory MLP activation-function sensitivity analysis
# ------------------------------------------------------------------------------
ACTIVATIONS <- c("logistic", "tanh", "relu")
ACTIVATION_HIDDEN <- c(1L, 2L, 3L, 4L, 5L, 6L, 8L, 10L)
ACTIVATION_DECAY <- c(0, 0.0001, 0.001, 0.01)
ACTIVATION_RESTARTS <- 3L

activation_forward <- function(par, X, hidden, activation) {
  p <- ncol(X)
  n_w1 <- (p + 1L) * hidden
  W1 <- matrix(par[seq_len(n_w1)], nrow = p + 1L, ncol = hidden)
  W2 <- par[(n_w1 + 1L):length(par)]
  Hlin <- cbind(1, X) %*% W1
  H <- switch(
    activation,
    logistic = 1 / (1 + exp(-pmax(pmin(Hlin, 35), -35))),
    tanh = tanh(Hlin),
    relu = pmax(0, Hlin),
    stop("Unknown activation function.")
  )
  as.numeric(cbind(1, H) %*% W2)
}

fit_custom_mlp_bfgs <- function(train, predictors, hidden, decay, activation, seed) {
  mu_x <- vapply(train[predictors], mean, numeric(1))
  sd_x <- vapply(train[predictors], sd, numeric(1))
  sd_x[!is.finite(sd_x) | sd_x == 0] <- 1
  mu_y <- mean(train$freshwater_discharge_m3_s)
  sd_y <- sd(train$freshwater_discharge_m3_s)
  if (!is.finite(sd_y) || sd_y == 0) return(NULL)
  X <- sweep(as.matrix(train[, predictors, drop = FALSE]), 2, mu_x, "-")
  X <- sweep(X, 2, sd_x, "/")
  y <- (train$freshwater_discharge_m3_s - mu_y) / sd_y
  p <- ncol(X)
  n_par <- (p + 1L) * hidden + hidden + 1L
  objective <- function(par) {
    pred <- activation_forward(par, X, hidden, activation)
    mse <- mean((y - pred)^2)
    n_w1 <- (p + 1L) * hidden
    W1 <- matrix(par[seq_len(n_w1)], nrow = p + 1L, ncol = hidden)
    W2 <- par[(n_w1 + 1L):length(par)]
    penalty <- sum(W1[-1, , drop = FALSE]^2) + sum(W2[-1]^2)
    mse + decay * penalty
  }
  best <- NULL
  for (r in seq_len(ACTIVATION_RESTARTS)) {
    set.seed(seed + r - 1L)
    init <- rnorm(n_par, sd = 0.15)
    fit <- tryCatch(
      optim(init, objective, method = "BFGS", control = list(maxit = 1000, reltol = 1e-8)),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && (is.null(best) || fit$value < best$value)) best <- fit
  }
  if (is.null(best)) return(NULL)
  list(par = best$par, predictors = predictors, hidden = hidden, activation = activation,
       mu_x = mu_x, sd_x = sd_x, mu_y = mu_y, sd_y = sd_y)
}

predict_custom_mlp_bfgs <- function(object, newdata) {
  X <- sweep(as.matrix(newdata[, object$predictors, drop = FALSE]), 2, object$mu_x, "-")
  X <- sweep(X, 2, object$sd_x, "/")
  z <- activation_forward(object$par, X, object$hidden, object$activation)
  z * object$sd_y + object$mu_y
}

cv_custom_activation_candidate <- function(df_cal, predictors, hidden, decay, activation,
                                           seed = SEED + 7000L) {
  origins <- rolling_origins(df_cal)
  fold_results <- map_dfr(seq_along(origins), function(j) {
    origin <- origins[j]
    train <- df_cal[seq_len(origin), , drop = FALSE]
    target <- df_cal[origin + 1L, , drop = FALSE]
    train_ok <- train %>% select(freshwater_discharge_m3_s, all_of(predictors)) %>% drop_na()
    target_ok <- complete.cases(target[, c("freshwater_discharge_m3_s", predictors), drop = FALSE])
    if (nrow(train_ok) < 30L || !target_ok) {
      return(tibble(Fold = j, date = target$date, Observed = target$freshwater_discharge_original_m3_s,
                    Predicted = NA_real_, Success = FALSE))
    }
    fit <- fit_custom_mlp_bfgs(train_ok, predictors, hidden, decay, activation,
                               seed = seed + 100L * j)
    if (is.null(fit)) {
      return(tibble(Fold = j, date = target$date, Observed = target$freshwater_discharge_original_m3_s,
                    Predicted = NA_real_, Success = FALSE))
    }
    pr <- tryCatch(predict_custom_mlp_bfgs(fit, target), error = function(e) NA_real_)
    tibble(Fold = j, date = target$date, Observed = target$freshwater_discharge_original_m3_s,
           Predicted = pr, Success = is.finite(pr))
  })
  ok <- fold_results$Success & is.finite(fold_results$Observed) & is.finite(fold_results$Predicted)
  tibble(
    CV_RMSE = if (any(ok)) sqrt(mean((fold_results$Observed[ok] - fold_results$Predicted[ok])^2)) else NA_real_,
    CV_MAE = if (any(ok)) mean(abs(fold_results$Observed[ok] - fold_results$Predicted[ok])) else NA_real_,
    CV_N = sum(ok),
    CV_total = nrow(fold_results),
    CV_success = mean(ok)
  )
}

activation_scenarios <- tibble(
  Period = c("T1", "T1", "T2", "T2"),
  Information_domain = c("Discharge only", "Discharge + precipitation",
                         "Discharge only", "Discharge + precipitation"),
  Predictors = c("Q1+Q2+Q3", "Q1+Q2+Q3+P1+P2",
                 "Q1+Q2+Q3", "Q1+Q2+Q3+P1+P2")
)

activation_full_grid <- pmap_dfr(activation_scenarios, function(Period, Information_domain, Predictors) {
  cal_end <- if (Period == "T1") T1_CAL_END else T2_CAL_END
  df_cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  predictors <- strsplit(Predictors, "\\+")[[1]]
  crossing(Hidden = ACTIVATION_HIDDEN, Decay = ACTIVATION_DECAY,
           Activation = ACTIVATIONS) %>%
    pmap_dfr(function(Hidden, Decay, Activation) {
      bind_cols(
        tibble(Period = Period, Information_domain = Information_domain,
               Predictors = Predictors, Hidden = Hidden, Decay = Decay,
               Activation = Activation,
               Architecture = paste(length(predictors), Hidden, 1L, sep = "-")),
        cv_custom_activation_candidate(df_cal, predictors, Hidden, Decay, Activation,
                                       seed = SEED + ifelse(Period == "T2", 2000L, 0L))
      )
    })
}) %>%
  filter(is.finite(CV_RMSE), CV_success >= 0.90) %>%
  arrange(Period, Information_domain, CV_RMSE, CV_MAE)

best_model_by_activation <- activation_full_grid %>%
  group_by(Period, Information_domain, Activation) %>%
  arrange(CV_RMSE, CV_MAE, Hidden, Decay, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup()

activation_comparison <- best_model_by_activation %>%
  group_by(Period, Information_domain) %>%
  mutate(Best_CV_RMSE = min(CV_RMSE),
         Delta_CV_RMSE = CV_RMSE - Best_CV_RMSE,
         Percent_difference = 100 * Delta_CV_RMSE / Best_CV_RMSE) %>%
  ungroup()

best_activation_overall <- best_model_by_activation %>%
  group_by(Period, Information_domain) %>%
  arrange(CV_RMSE, CV_MAE, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup()

write_csv(activation_full_grid, file.path(OUTPUT_DIR, "mlp_activation_full_grid.csv"))
write_csv(best_model_by_activation, file.path(OUTPUT_DIR, "mlp_best_model_by_activation.csv"))
write_csv(activation_comparison, file.path(OUTPUT_DIR, "mlp_activation_comparison.csv"))
write_csv(best_activation_overall, file.path(OUTPUT_DIR, "mlp_best_activation_overall.csv"))

# Compare the custom logistic implementation with nnet using the same predictors,
# hidden size, decay and rolling-origin folds. This comparison is diagnostic only.
compare_logistic_engines <- function(period, information_domain, predictors_text, hidden, decay) {
  cal_end <- if (period == "T1") T1_CAL_END else T2_CAL_END
  df_cal <- lagged_data %>% filter(date >= T1_CAL_START, date <= cal_end)
  predictors <- strsplit(predictors_text, "\\+")[[1]]
  origins <- rolling_origins(df_cal)
  folds <- map_dfr(seq_along(origins), function(j) {
    origin <- origins[j]
    train <- df_cal[seq_len(origin), , drop = FALSE]
    target <- df_cal[origin + 1L, , drop = FALSE]
    train_ok <- train %>% select(freshwater_discharge_m3_s, all_of(predictors)) %>% drop_na()
    if (nrow(train_ok) < 30L || !all(complete.cases(target[, c("freshwater_discharge_m3_s", predictors), drop = FALSE]))) {
      return(tibble(Fold = j, date = target$date, Observed = target$freshwater_discharge_original_m3_s,
                    Pred_nnet = NA_real_, Pred_custom = NA_real_))
    }
    nnet_fit <- tryCatch(fit_mlp(train_ok, predictors, hidden, decay,
                                 seed = SEED + 100L * j), error = function(e) NULL)
    custom_fit <- fit_custom_mlp_bfgs(train_ok, predictors, hidden, decay, "logistic",
                                      seed = SEED + 100L * j)
    p_nnet <- if (is.null(nnet_fit)) NA_real_ else
      tryCatch(predict_mlp(nnet_fit, target), error = function(e) NA_real_)
    p_custom <- if (is.null(custom_fit)) NA_real_ else
      tryCatch(predict_custom_mlp_bfgs(custom_fit, target), error = function(e) NA_real_)
    tibble(Fold = j, date = target$date, Observed = target$freshwater_discharge_original_m3_s,
           Pred_nnet = p_nnet, Pred_custom = p_custom)
  })
  ok <- with(folds, is.finite(Observed) & is.finite(Pred_nnet) & is.finite(Pred_custom))
  summary <- tibble(
    Period = period, Information_domain = information_domain, N = sum(ok),
    RMSE_nnet = sqrt(mean((folds$Observed[ok] - folds$Pred_nnet[ok])^2)),
    RMSE_custom = sqrt(mean((folds$Observed[ok] - folds$Pred_custom[ok])^2)),
    MAE_nnet = mean(abs(folds$Observed[ok] - folds$Pred_nnet[ok])),
    MAE_custom = mean(abs(folds$Observed[ok] - folds$Pred_custom[ok])),
    Prediction_correlation = cor(folds$Pred_nnet[ok], folds$Pred_custom[ok]),
    Mean_absolute_prediction_difference = mean(abs(folds$Pred_nnet[ok] - folds$Pred_custom[ok])),
    RMSE_percent_difference = 100 * abs(RMSE_custom - RMSE_nnet) / RMSE_nnet
  )
  list(folds = folds %>% mutate(Period = period, Information_domain = information_domain,
                                Hidden = hidden, Decay = decay), summary = summary)
}

best_logistic <- best_model_by_activation %>% filter(Activation == "logistic")
logistic_engine_results <- pmap(best_logistic,
  function(Period, Information_domain, Predictors, Hidden, Decay, Activation,
           Architecture, CV_RMSE, CV_MAE, CV_N, CV_total, CV_success, ...) {
    compare_logistic_engines(Period, Information_domain, Predictors, Hidden, Decay)
  })
logistic_engine_folds <- map_dfr(logistic_engine_results, "folds")
logistic_engine_summary <- map_dfr(logistic_engine_results, "summary")
write_csv(logistic_engine_folds, file.path(OUTPUT_DIR, "logistic_nnet_vs_custom_folds.csv"))
write_csv(logistic_engine_summary, file.path(OUTPUT_DIR, "logistic_nnet_vs_custom_summary.csv"))


# ------------------------------------------------------------------------------
# 20. Exploratory learning-rate sensitivity analysis (torch, full-batch SGD)
# ------------------------------------------------------------------------------
RUN_LEARNING_RATE_SENSITIVITY <- TRUE

if (RUN_LEARNING_RATE_SENSITIVITY && !requireNamespace("torch", quietly = TRUE)) {
  warning(
    "Package 'torch' is not installed. The learning-rate sensitivity section will be skipped; ",
    "all other analyses will still run. Install 'torch' to reproduce this exploratory check."
  )
  RUN_LEARNING_RATE_SENSITIVITY <- FALSE
}

if (RUN_LEARNING_RATE_SENSITIVITY && !torch::torch_is_installed()) {
  warning(
    "The torch backend is not installed. The learning-rate sensitivity section will be skipped. ",
    "Run torch::install_torch(), restart R, and rerun the script to reproduce this exploratory check."
  )
  RUN_LEARNING_RATE_SENSITIVITY <- FALSE
}

if (RUN_LEARNING_RATE_SENSITIVITY) {
  suppressPackageStartupMessages(library(torch))

  if (!exists("lagged_data")) {
    stop(
      "Object 'lagged_data' is not available. Run the primary workflow first."
    )
  }


  if (!exists("T1_CAL_START")) {
    T1_CAL_START <- as.Date("1966-01-01")
  }

  if (!exists("T1_CAL_END")) {
    T1_CAL_END <- as.Date("2005-12-01")
  }

  if (!exists("T2_CAL_END")) {
    T2_CAL_END <- as.Date("2011-12-01")
  }

  if (!exists("OUTPUT_DIR")) {
    stop("Object 'OUTPUT_DIR' is not available.")
  }


  if (!exists("CV_INITIAL_YEARS")) {
    CV_INITIAL_YEARS <- 20L
  }

  if (!exists("CV_LAST_MONTHS")) {
    CV_LAST_MONTHS <- 120L
  }

  if (!exists("CV_STEP")) {
    CV_STEP <- 3L
  }


  learning_rate_origins <- function(df_cal) {
  
    first_origin <- max(
      CV_INITIAL_YEARS * 12L,
      nrow(df_cal) - CV_LAST_MONTHS
    )
  
    seq(
      first_origin,
      nrow(df_cal) - 1L,
      by = CV_STEP
    )
  }


  LEARNING_RATE_SEED <- 20260826L

  LEARNING_RATE_RESTARTS <- 3L


  LEARNING_RATE_MAX_EPOCHS <- 2000L

  LEARNING_RATE_MIN_EPOCHS <- 300L

  LEARNING_RATE_CHECK_EVERY <- 50L

  LEARNING_RATE_TOLERANCE <- 1e-6


  LEARNING_RATE_GRID <- c(
    0.0001,
    0.0005,
    0.001,
    0.005,
    0.01,
    0.05,
    0.10
  )


  LEARNING_RATE_DIR <- file.path(
    OUTPUT_DIR,
    "learning_rate_sensitivity"
  )

  dir.create(
    LEARNING_RATE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )


  learning_rate_config <- tibble(
    Period = c(
      "T1",
      "T1",
      "T2",
      "T2"
    ),
  
    Information_domain = c(
      "Discharge only",
      "Discharge + precipitation",
      "Discharge only",
      "Discharge + precipitation"
    ),
  
    Predictors = c(
      "Q1+Q2+Q3",
      "Q1+Q2+Q3+P1+P2",
      "Q1+Q2+Q3",
      "Q1+Q2+Q3+P1+P2"
    ),
  
    Hidden = c(
      8L,
      5L,
      6L,
      5L
    ),
  
    Weight_decay = c(
      0.0001,
      0.01,
      0.01,
      0
    )
  )

  print(learning_rate_config, n = Inf)


  fit_sgd_once <- function(
      training,
      test_result,
      predictors,
      hidden,
      weight_decay,
      learning_rate,
      seed
  ) {
  
  
    mu_x <- vapply(
      training[predictors],
      mean,
      numeric(1)
    )
  
    sd_x <- vapply(
      training[predictors],
      sd,
      numeric(1)
    )
  
    sd_x[
      !is.finite(sd_x) |
        sd_x == 0
    ] <- 1
  
  
    mu_y <- mean(
      training$freshwater_discharge_m3_s
    )
  
    sd_y <- sd(
      training$freshwater_discharge_m3_s
    )
  
    if (
      !is.finite(sd_y) ||
      sd_y == 0
    ) {
      return(
        list(
          Success = FALSE,
          Predicted = NA_real_,
          Epochs = NA_integer_,
          Converged = FALSE,
          Final_loss = NA_real_
        )
      )
    }
  
  
    X_train <- as.matrix(
      training[, predictors, drop = FALSE]
    )
  
    X_test <- as.matrix(
      test_result[, predictors, drop = FALSE]
    )
  
  
    X_train <- sweep(
      X_train,
      2,
      mu_x,
      "-"
    )
  
    X_train <- sweep(
      X_train,
      2,
      sd_x,
      "/"
    )
  
  
    X_test <- sweep(
      X_test,
      2,
      mu_x,
      "-"
    )
  
    X_test <- sweep(
      X_test,
      2,
      sd_x,
      "/"
    )
  
  
    y_train <- (
      training$freshwater_discharge_m3_s - mu_y
    ) / sd_y
  
  
  
    x_train_t <- torch_tensor(
      X_train,
      dtype = torch_float()
    )
  
    y_train_t <- torch_tensor(
      matrix(
        y_train,
        ncol = 1L
      ),
      dtype = torch_float()
    )
  
    x_test_t <- torch_tensor(
      X_test,
      dtype = torch_float()
    )
  
  
  
    set.seed(seed)
  
    torch_manual_seed(
      as.integer(seed)
    )
  
  
  
    model <- nn_sequential(
      nn_linear(
        length(predictors),
        hidden
      ),
    
      nn_sigmoid(),
    
      nn_linear(
        hidden,
        1L
      )
    )
  
  
    optimizer <- optim_sgd(
      model$parameters,
      lr = learning_rate,
      momentum = 0,
      weight_decay = weight_decay
    )
  
  
    model$train()
  
  
    previous_loss_check <- NA_real_
  
    converged_flag <- FALSE
  
    epochs_used <- LEARNING_RATE_MAX_EPOCHS
  
    final_loss <- NA_real_
  
  
  
    for (epoch in seq_len(LEARNING_RATE_MAX_EPOCHS)) {
    
      optimizer$zero_grad()
    
      pred_train <- model(
        x_train_t
      )
    
      loss <- nnf_mse_loss(
        pred_train,
        y_train_t
      )
    
      loss_value <- as.numeric(
        as_array(loss)
      )
    
    
      if (!is.finite(loss_value)) {
      
        return(
          list(
            Success = FALSE,
            Predicted = NA_real_,
            Epochs = epoch,
            Converged = FALSE,
            Final_loss = NA_real_
          )
        )
      }
    
    
      loss$backward()
    
      optimizer$step()
    
    
    
      if (
        epoch >= LEARNING_RATE_MIN_EPOCHS &&
        epoch %% LEARNING_RATE_CHECK_EVERY == 0L
      ) {
      
        if (is.finite(previous_loss_check)) {
        
          relative_change <- abs(
            previous_loss_check -
              loss_value
          ) /
            max(
              abs(previous_loss_check),
              1e-8
            )
        
        
          if (
            is.finite(relative_change) &&
            relative_change < LEARNING_RATE_TOLERANCE
          ) {
          
            converged_flag <- TRUE
          
            epochs_used <- epoch
          
            final_loss <- loss_value
          
            break
          }
        }
      
      
        previous_loss_check <- loss_value
      }
    
    
      if (epoch == LEARNING_RATE_MAX_EPOCHS) {
      
        final_loss <- loss_value
      
        epochs_used <- epoch
      }
    }
  
  
  
    model$eval()
  
  
    pred_z <- with_no_grad({
    
      as.numeric(
        as_array(
          model(
            x_test_t
          )
        )
      )
    
    })
  
  
    prediction_original_scale <- (
      pred_z * sd_y
    ) + mu_y
  
  
    success_flag <- (
      length(prediction_original_scale) == 1L &&
        is.finite(prediction_original_scale)
    )
  
  
    list(
      Success = success_flag,
      Predicted = ifelse(
        success_flag,
        prediction_original_scale,
        NA_real_
      ),
      Epochs = epochs_used,
      Converged = converged_flag,
      Final_loss = final_loss
    )
  }


  evaluate_learning_rate <- function(
      period,
      information_domain,
      predictors,
      hidden,
      weight_decay,
      learning_rate,
      lr_id,
      scenario_id
  ) {
  
    cal_end <- if (
      period == "T1"
    ) {
      T1_CAL_END
    } else {
      T2_CAL_END
    }
  
  
    df <- lagged_data |>
      filter(
        date >= T1_CAL_START,
        date <= cal_end
      ) |>
      select(
        date,
        freshwater_discharge_m3_s,
        freshwater_discharge_observed,
        all_of(predictors)
      )
  
  
    idx_cv <- learning_rate_origins(
      df
    )
  
  
    map_dfr(
      seq_along(idx_cv),
      function(fold_id) {
      
        i <- idx_cv[
          fold_id
        ]
      
      
        training <- df[
          seq_len(i),
          ,
          drop = FALSE
        ] |>
          filter(
            is.finite(freshwater_discharge_m3_s),
            if_all(
              all_of(predictors),
              is.finite
            )
          )
      
      
        test_result <- df[
          i + 1L,
          ,
          drop = FALSE
        ]
      
      
        eligible_flag <- (
          isTRUE(
            test_result$freshwater_discharge_observed[[1]]
          ) &&
            is.finite(
              test_result$freshwater_discharge_m3_s[[1]]
            ) &&
            all(
              is.finite(
                unlist(
                  test_result[
                    ,
                    predictors,
                    drop = FALSE
                  ]
                )
              )
            )
        )
      
      
        if (
          !eligible_flag ||
          nrow(training) <=
          length(predictors) + 2L
        ) {
        
          return(
            tibble(
              Period = period,
              Information_domain = information_domain,
              Learning_rate = learning_rate,
              Hidden = hidden,
              Weight_decay = weight_decay,
            
              Fold = fold_id,
              Origin = df$date[i],
              Target = test_result$date[[1]],
            
              Eligible = eligible_flag,
            
              Observed = ifelse(
                eligible_flag,
                test_result$freshwater_discharge_m3_s[[1]],
                NA_real_
              ),
            
              Predicted = NA_real_,
              Error = NA_real_,
              Absolute_error = NA_real_,
            
              Fold_success = FALSE,
              Mean_epochs = NA_real_,
              Fraction_plateau_converged = NA_real_,
              Fraction_reached_max_epochs = NA_real_
            )
          )
        }
      
      
      
        reps <- map(
          seq_len(LEARNING_RATE_RESTARTS),
          function(rep_id) {
          
            current_seed <-
              LEARNING_RATE_SEED +
              scenario_id * 100000L +
              rep_id * 1000L +
              i
          
          
            tryCatch(
              fit_sgd_once(
                training = training,
                test_result = test_result,
                predictors = predictors,
                hidden = hidden,
                weight_decay = weight_decay,
                learning_rate = learning_rate,
                seed = current_seed
              ),
              error = function(e) {
              
                list(
                  Success = FALSE,
                  Predicted = NA_real_,
                  Epochs = NA_integer_,
                  Converged = FALSE,
                  Final_loss = NA_real_
                )
              }
            )
          }
        )
      
      
        pred_reps <- map_dbl(
          reps,
          "Predicted"
        )
      
      
        success_reps <- is.finite(
          pred_reps
        )
      
      
        fold_success <- (
          mean(success_reps) >= 0.90
        )
      
      
        predicted_value <- if (
          fold_success
        ) {
        
          mean(
            pred_reps[
              success_reps
            ]
          )
        
        } else {
        
          NA_real_
        }
      
      
        error_value <- if (
          is.finite(predicted_value)
        ) {
          test_result$freshwater_discharge_m3_s[[1]] -
            predicted_value
        } else {
          NA_real_
        }
      
      
        epochs_reps <- map_dbl(
          reps,
          function(x) {
            as.numeric(
              x$Epochs
            )
          }
        )
      
      
        conv_reps <- map_lgl(
          reps,
          function(x) {
            isTRUE(
              x$Converged
            )
          }
        )
      
      
        tibble(
          Period = period,
          Information_domain = information_domain,
        
          Learning_rate = learning_rate,
          Hidden = hidden,
          Weight_decay = weight_decay,
        
          Fold = fold_id,
          Origin = df$date[i],
          Target = test_result$date[[1]],
        
          Eligible = TRUE,
        
          Observed = test_result$freshwater_discharge_m3_s[[1]],
          Predicted = predicted_value,
        
          Error = error_value,
        
          Absolute_error = ifelse(
            is.finite(error_value),
            abs(error_value),
            NA_real_
          ),
        
          Fold_success = is.finite(
            predicted_value
          ),
        
          Mean_epochs = ifelse(
            any(is.finite(epochs_reps)),
            mean(
              epochs_reps,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
          Fraction_plateau_converged =
            mean(
              conv_reps
            ),

          Fraction_reached_max_epochs =
            mean(
              epochs_reps >= LEARNING_RATE_MAX_EPOCHS,
              na.rm = TRUE
            )
        )
      }
    )
  }


  learning_rate_grid_cases <- crossing(
    Scenario_id = seq_len(
      nrow(learning_rate_config)
    ),
  
    LR_id = seq_along(
      LEARNING_RATE_GRID
    )
  ) |>
    mutate(
      Learning_rate =
        LEARNING_RATE_GRID[
          LR_id
        ]
    )


  learning_rate_fold_results <- pmap_dfr(
    learning_rate_grid_cases,
    function(
      Scenario_id,
      LR_id,
      Learning_rate
    ) {
    
      cfg <- learning_rate_config[
        Scenario_id,
        ,
        drop = FALSE
      ]
    
    
      predictors <- strsplit(
        cfg$Predictors[[1]],
        "\\+"
      )[[1]]
    
    
      message(
        "Running: ",
        cfg$Period[[1]],
        " | ",
        cfg$Information_domain[[1]],
        " | lr = ",
        Learning_rate
      )
    
    
      evaluate_learning_rate(
        period =
          cfg$Period[[1]],
      
        information_domain =
          cfg$Information_domain[[1]],
      
        predictors =
          predictors,
      
        hidden =
          cfg$Hidden[[1]],
      
        weight_decay =
          cfg$Weight_decay[[1]],
      
        learning_rate =
          Learning_rate,
      
        lr_id =
          LR_id,
      
        scenario_id =
          Scenario_id
      )
    }
  )


  learning_rate_summary <- learning_rate_fold_results |>
    group_by(
      Period,
      Information_domain,
      Learning_rate,
      Hidden,
      Weight_decay
    ) |>
    summarise(
    
      CV_total = sum(
        Eligible
      ),
    
      CV_N = sum(
        Fold_success
      ),
    
      CV_success =
        CV_N /
        CV_total,
    
      CV_RMSE = ifelse(
        CV_N > 0,
        sqrt(
          mean(
            Error[
              Fold_success
            ]^2,
            na.rm = TRUE
          )
        ),
        NA_real_
      ),
    
      CV_MAE = ifelse(
        CV_N > 0,
        mean(
          Absolute_error[
            Fold_success
          ],
          na.rm = TRUE
        ),
        NA_real_
      ),
    
      Mean_epochs = mean(
        Mean_epochs[
          Fold_success
        ],
        na.rm = TRUE
      ),
    
      Mean_plateau_convergence =
        mean(
          Fraction_plateau_converged[
            Fold_success
          ],
          na.rm = TRUE
        ),

      Mean_fraction_reached_max_epochs =
        mean(
          Fraction_reached_max_epochs[
            Fold_success
          ],
          na.rm = TRUE
        ),
    
      .groups = "drop"
    ) |>
    mutate(
      Eligible_model = (
        CV_success >= 0.90
      )
    ) |>
    arrange(
      Period,
      Information_domain,
      CV_RMSE
    )


  print(
    learning_rate_summary,
    n = Inf
  )


  best_learning_rate <- learning_rate_summary |>
    filter(
      Eligible_model,
      is.finite(CV_RMSE)
    ) |>
    group_by(
      Period,
      Information_domain
    ) |>
    arrange(
      CV_RMSE,
      CV_MAE
    ) |>
    slice(1L) |>
    ungroup()


  print(
    best_learning_rate,
    n = Inf
  )


  learning_rate_reference <- learning_rate_summary |>
    filter(
      Learning_rate == 0.01
    ) |>
    select(
      Period,
      Information_domain,
      RMSE_LR_0.01 = CV_RMSE
    )


  learning_rate_comparison <- learning_rate_summary |>
    left_join(
      learning_rate_reference,
      by = c(
        "Period",
        "Information_domain"
      )
    ) |>
    mutate(
      Delta_RMSE_vs_0.01_pct =
        100 *
        (
          CV_RMSE -
            RMSE_LR_0.01
        ) /
        RMSE_LR_0.01
    )


  BOOT_LR_B <- 2000L
  BOOT_LR_BLOCK_FOLDS <- 4L
  RMSE_EQUIVALENCE_PCT <- 2

  sample_circular_block_indices <- function(n, block_length) {
    n_blocks <- ceiling(n / block_length)
    starts <- sample.int(n, n_blocks, replace = TRUE)

    idx <- unlist(
      lapply(
        starts,
        function(s) {
          ((s - 1L + seq_len(block_length) - 1L) %% n) + 1L
        }
      ),
      use.names = FALSE
    )

    idx[seq_len(n)]
  }

  compare_learning_rate_paired <- function(df_scenario) {
    reference_results <- df_scenario |>
      filter(
        Learning_rate == 0.01,
        Fold_success,
        is.finite(Error)
      ) |>
      select(Fold, Target, Error_ref = Error)

    learning_rates <- sort(unique(df_scenario$Learning_rate))

    map_dfr(
      learning_rates,
      function(current_lr) {
        candidate_results <- df_scenario |>
          filter(
            Learning_rate == current_lr,
            Fold_success,
            is.finite(Error)
          ) |>
          select(Fold, Target, Error_candidate = Error)

        paired_results <- inner_join(
          candidate_results,
          reference_results,
          by = c("Fold", "Target")
        ) |>
          arrange(Target)

        n <- nrow(paired_results)

        if (n < 8L) {
          return(
            tibble(
              Learning_rate = current_lr,
              Reference_learning_rate = 0.01,
              N_pairs = n,
              RMSE_candidate = NA_real_,
              RMSE_reference = NA_real_,
              Delta_RMSE = NA_real_,
              Delta_RMSE_pct = NA_real_,
              Lower95_Delta_RMSE = NA_real_,
              Upper95_Delta_RMSE = NA_real_,
              MAE_candidate = NA_real_,
              MAE_reference = NA_real_,
              Delta_MAE = NA_real_,
              Lower95_Delta_MAE = NA_real_,
              Upper95_Delta_MAE = NA_real_,
              Interpretation = "Insufficient paired folds"
            )
          )
        }

        rmse_candidate <- sqrt(mean(paired_results$Error_candidate^2))
        rmse_reference <- sqrt(mean(paired_results$Error_ref^2))
        delta_rmse <- rmse_candidate - rmse_reference
        delta_rmse_pct <- 100 * delta_rmse / rmse_reference

        mae_candidate <- mean(abs(paired_results$Error_candidate))
        mae_reference <- mean(abs(paired_results$Error_ref))
        delta_mae <- mae_candidate - mae_reference

        if (current_lr == 0.01) {
          boot_delta_rmse <- rep(0, BOOT_LR_B)
          boot_delta_mae <- rep(0, BOOT_LR_B)
        } else {
          boot_stats <- replicate(
            BOOT_LR_B,
            {
              idx <- sample_circular_block_indices(
                n,
                BOOT_LR_BLOCK_FOLDS
              )

              e_candidate <- paired_results$Error_candidate[idx]
              e_ref <- paired_results$Error_ref[idx]

              c(
                Delta_RMSE =
                  sqrt(mean(e_candidate^2)) -
                  sqrt(mean(e_ref^2)),
                Delta_MAE =
                  mean(abs(e_candidate)) -
                  mean(abs(e_ref))
              )
            }
          )

          boot_delta_rmse <- boot_stats["Delta_RMSE", ]
          boot_delta_mae <- boot_stats["Delta_MAE", ]
        }

        ci_rmse <- quantile(
          boot_delta_rmse,
          probs = c(0.025, 0.975),
          names = FALSE,
          type = 8
        )

        ci_mae <- quantile(
          boot_delta_mae,
          probs = c(0.025, 0.975),
          names = FALSE,
          type = 8
        )

        margem_equivalencia <-
          RMSE_EQUIVALENCE_PCT / 100 * rmse_reference

        interpretacao <- case_when(
          current_lr == 0.01 ~ "Reference",
          ci_rmse[2] < 0 ~ "Lower RMSE than reference",
          ci_rmse[1] > 0 ~ "Higher RMSE than reference",
          ci_rmse[1] >= -margem_equivalencia &
            ci_rmse[2] <= margem_equivalencia ~
            "Practically equivalent within 2% RMSE",
          TRUE ~ "Inconclusive difference"
        )

        tibble(
          Learning_rate = current_lr,
          Reference_learning_rate = 0.01,
          N_pairs = n,
          RMSE_candidate = rmse_candidate,
          RMSE_reference = rmse_reference,
          Delta_RMSE = delta_rmse,
          Delta_RMSE_pct = delta_rmse_pct,
          Lower95_Delta_RMSE = ci_rmse[1],
          Upper95_Delta_RMSE = ci_rmse[2],
          MAE_candidate = mae_candidate,
          MAE_reference = mae_reference,
          Delta_MAE = delta_mae,
          Lower95_Delta_MAE = ci_mae[1],
          Upper95_Delta_MAE = ci_mae[2],
          Interpretation = interpretacao
        )
      }
    )
  }

  set.seed(LEARNING_RATE_SEED + 900000L)

  learning_rate_inference <- learning_rate_fold_results |>
    group_by(Period, Information_domain) |>
    group_modify(~ compare_learning_rate_paired(.x)) |>
    ungroup()

  print(learning_rate_inference, n = Inf)


  write_csv(
    learning_rate_fold_results,
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_fold_results.csv"
    )
  )


  write_csv(
    learning_rate_summary,
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_summary.csv"
    )
  )


  write_csv(
    best_learning_rate,
    file.path(
      LEARNING_RATE_DIR,
      "best_learning_rate.csv"
    )
  )


  write_csv(
    learning_rate_comparison,
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_comparison.csv"
    )
  )

  write_csv(
    learning_rate_inference,
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_paired_inference.csv"
    )
  )


  fig_learning_rate <- learning_rate_summary |>
    filter(
      is.finite(CV_RMSE)
    ) |>
    ggplot(
      aes(
        x = Learning_rate,
        y = CV_RMSE
      )
    ) +
    geom_line(
      linewidth = 0.6
    ) +
    geom_point(
      size = 2.3
    ) +
    scale_x_log10() +
    facet_grid(
      Information_domain ~ Period,
      scales = "free_y"
    ) +
    labs(
      x = "Learning rate",
      y = expression(
        Rolling-origin~CV~RMSE~(m^3~s^{-1})
      )
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      strip.text = element_text(
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    )


  ggsave(
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_sensitivity.png"
    ),
    fig_learning_rate,
    width = 8.5,
    height = 6.5,
    dpi = 600
  )

  fig_learning_rate_delta <- learning_rate_inference |>
    filter(
      Learning_rate != 0.01,
      is.finite(Delta_RMSE)
    ) |>
    ggplot(
      aes(
        x = Learning_rate,
        y = Delta_RMSE
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_errorbar(
      aes(
        ymin = Lower95_Delta_RMSE,
        ymax = Upper95_Delta_RMSE
      ),
      width = 0
    ) +
    geom_point(size = 2.3) +
    scale_x_log10() +
    facet_grid(
      Information_domain ~ Period,
      scales = "free_y"
    ) +
    labs(
      x = "Learning rate",
      y = expression(
        Delta~RMSE~versus~LR==0.01~(m^3~s^{-1})
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_paired_delta.png"
    ),
    fig_learning_rate_delta,
    width = 8.5,
    height = 6.5,
    dpi = 600,
    bg = "white"
  )


  saveRDS(
    list(
      learning_rate_config = learning_rate_config,
      learning_rate_fold_results =
        learning_rate_fold_results,
      learning_rate_summary =
        learning_rate_summary,
      best_learning_rate =
        best_learning_rate,
      learning_rate_comparison =
        learning_rate_comparison,
      learning_rate_inference =
        learning_rate_inference
    ),
    file.path(
      LEARNING_RATE_DIR,
      "learning_rate_analysis_objects.rds"
    )
  )


  message(
    "Learning-rate sensitivity analysis completed."
  )
}


# ------------------------------------------------------------------------------
# 21. Final reproducibility record
# ------------------------------------------------------------------------------
package_versions <- tibble(
  package = packages,
  version = vapply(
    packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  )
)
write_csv(package_versions, file.path(OUTPUT_DIR, "package_versions.csv"))
writeLines(as.character(SEED), file.path(OUTPUT_DIR, "primary_random_seed.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info_complete.txt"))
message("Complete analysis workflow finished. Output directory: ", OUTPUT_DIR)
