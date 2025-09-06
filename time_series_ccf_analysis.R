# ---- Config ----
TZ_USED <- "UTC"                                   # e.g., "UTC" or "America/New_York"
ts_file_path <- "path/to/timeseries.csv"           # time-series CSV to analyze
window_ref_path <- "path/to/window_reference.csv"  # table with ID + start/end dates
out_dir <- "ccf_outputs"                           # where to save PNGs and CSV
max_lag_hours <- 1.5                                # reporting window (±hours)

# Patterns to identify the two series columns in ts_file_path
# Adjust these to match your dataset's column naming
vector1_pattern <- "vector1"                        # e.g., "temperature.*setpoint"
vector2_pattern <- "vector2"                        # e.g., "valve_feedback|flow"

# Ensure output directory exists
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Build a label from the input filename (e.g., "FBH18C" -> "FBH18C")
label <- tools::file_path_sans_ext(basename(ts_file_path))

# ---- Helpers ----
normalize_id <- function(x) toupper(trimws(x))

parse_window <- function(start_str, end_str, tz = TZ_USED) {
  st_date <- suppressWarnings(lubridate::mdy(start_str, tz = tz))
  if (is.na(st_date)) st_date <- suppressWarnings(lubridate::ymd(start_str, tz = tz))
  et_date <- suppressWarnings(lubridate::mdy(end_str, tz = tz))
  if (is.na(et_date)) et_date <- suppressWarnings(lubridate::ymd(end_str, tz = tz))

  if (any(is.na(c(st_date, et_date))))
    stop("Unable to parse start/end dates for noon-anchored window.")

  st_dt <- lubridate::as_datetime(lubridate::floor_date(st_date, "day") + lubridate::hours(12), tz = tz)  # start 12:00
  et_dt <- lubridate::as_datetime(lubridate::floor_date(et_date, "day") + lubridate::hours(12), tz = tz)  # end   12:00
  list(start = st_dt, end = et_dt)
}

# ---- 1) Load window reference, normalize IDs, locate window ----
ref_tbl <- read.csv(window_ref_path, check.names = FALSE, stringsAsFactors = FALSE)
names(ref_tbl) <- make.unique(names(ref_tbl))

# Expect a column that names the series/window "id"
if (!("id" %in% names(ref_tbl))) stop("Column 'id' not found in window reference.")

# Try to find start/end columns in a general way
start_col <- grep("^(?i).*start.*date", names(ref_tbl), value = TRUE)[1]
end_col   <- grep("^(?i).*end.*date",   names(ref_tbl), value = TRUE)[1]
if (is.na(start_col) || is.na(end_col)) {
  stop("Could not find start/end date columns in window reference (looked for .*start.*date / .*end.*date).")
}

ref_tbl <- dplyr::mutate(ref_tbl, id_norm = normalize_id(id))

# Attempt to infer an ID from the timeseries filename (fallback to whole stem if no token)
id_from_name <- stringr::str_extract(basename(ts_file_path), "[A-Za-z0-9]+")
if (is.na(id_from_name)) id_from_name <- label
id_norm <- normalize_id(id_from_name)

if (!(id_norm %in% ref_tbl$id_norm)) {
  stop("Could not find normalized id '", id_norm,
       "' in window reference. Candidates: ", paste(unique(ref_tbl$id_norm), collapse = ", "))
}

ref_row <- dplyr::filter(ref_tbl, id_norm == id_norm) |> dplyr::slice(1)
win <- parse_window(ref_row[[start_col]], ref_row[[end_col]], tz = TZ_USED)
start_time <- win$start; end_time <- win$end
if (!is.finite(start_time) || !is.finite(end_time) || start_time >= end_time) {
  stop("Invalid or unparsable analysis window for id '", id_norm,
       "'. Start: ", ref_row[[start_col]], " | End: ", ref_row[[end_col]])
}

# ---- 2) Load time series, parse timestamps, and clip to window ----
df <- readr::read_csv(ts_file_path, show_col_types = FALSE)
names(df) <- make.unique(names(df))

# Expect a timestamp column named 't_stamp' (general but explicit)
if (!("t_stamp" %in% names(df))) stop("Column 't_stamp' not found in the time-series file.")

t_try <- lubridate::parse_date_time(
  df$t_stamp,
  orders = c("mdy HM","mdy HMS","mdy IMp","mdy IMS p","ymd HM","ymd HMS","ymd IMp","ymd IMS p"),
  tz = TZ_USED
)
if (all(is.na(t_try))) t_try <- suppressWarnings(lubridate::mdy_hm(df$t_stamp, tz = TZ_USED))

df <- df |>
  dplyr::mutate(t_stamp = t_try) |>
  dplyr::arrange(t_stamp) |>
  dplyr::filter(!is.na(t_stamp), t_stamp >= start_time, t_stamp <= end_time)

if (nrow(df) < 10) stop("After windowing, too few rows to analyze for id '", id_norm, "'.")

# ---- 3) Identify the two series columns (vector1 / vector2) ----
v1_col <- grep(vector1_pattern, names(df), value = TRUE, ignore.case = TRUE)
v2_col <- grep(vector2_pattern, names(df), value = TRUE, ignore.case = TRUE)

if (length(v1_col) != 1) stop("Could not uniquely identify vector1 column. Found: ", paste(v1_col, collapse = ", "))
if (length(v2_col) != 1) stop("Could not uniquely identify vector2 column. Found: ", paste(v2_col, collapse = ", "))

# ---- 4) Interpolate missing values & z-score within the window ----
vector1_z <- scale(zoo::na.approx(df[[v1_col]], na.rm = FALSE))
vector2_z <- scale(zoo::na.approx(df[[v2_col]], na.rm = FALSE))

# ---- 5) Sampling cadence (MINUTES per sample) & ±max_lag_hours window ----
dt_min <- median(
  as.numeric(difftime(df$t_stamp[-1], df$t_stamp[-nrow(df)], units = "mins")),
  na.rm = TRUE
)
if (!is.finite(dt_min) || dt_min <= 0) dt_min <- 1
lag_window_samples <- max(1L, round(max_lag_hours * 60 / dt_min))  # ±max_lag_hours

# ---- 6) CCF (raw) ----
ccf_vals <- ccf(as.numeric(vector2_z), as.numeric(vector1_z),  # order: x leads y interpretation
                lag.max = lag_window_samples, plot = FALSE)

# Constrain report to ±lag_window_samples
in_win <- which(abs(ccf_vals$lag) <= lag_window_samples)
lag_use <- ccf_vals$lag[in_win]
acf_use <- ccf_vals$acf[in_win]

lag_at_max <- lag_use[which.max(abs(acf_use))]
max_corr   <- max(abs(acf_use), na.rm = TRUE)
lag_hours  <- lag_at_max * dt_min / 60

# Save full-range raw CCF plot
png(file.path(out_dir, sprintf("ccf_%s.png", label)), width = 1400, height = 800, res = 150)
plot(ccf_vals,
     main = paste0("CCF (", id_norm, "): ", v2_col, " vs ", v1_col,
                   sprintf("  [report window ±%.1f h]", lag_window_samples*dt_min/60)),
     xlab = "Lag (samples)", ylab = "CCF")
abline(v = lag_at_max, col = "red", lty = 2)
dev.off()

# ---- 7) Prewhitened CCF (AR(1) residuals) ----
x <- as.numeric(vector2_z); y <- as.numeric(vector1_z)
fit_x <- try(stats::arima(x, order = c(1,0,0)), silent = TRUE)
fit_y <- try(stats::arima(y, order = c(1,0,0)), silent = TRUE)
rx <- if (inherits(fit_x, "try-error")) x else resid(fit_x)
ry <- if (inherits(fit_y, "try-error")) y else resid(fit_y)

ccf_pw_full <- ccf(rx, ry, lag.max = lag_window_samples, plot = FALSE)

# Constrain prewhitened report to ±lag_window_samples
in_win_pw <- which(abs(ccf_pw_full$lag) <= lag_window_samples)
lag_use_pw <- ccf_pw_full$lag[in_win_pw]
acf_use_pw <- ccf_pw_full$acf[in_win_pw]

lag_at_max_pw <- lag_use_pw[which.max(abs(acf_use_pw))]
max_corr_pw   <- max(abs(acf_use_pw), na.rm = TRUE)

# 95% band on prewhitened residuals
n_obs_pw <- min(sum(is.finite(rx)), sum(is.finite(ry)))
N_eff_pw <- max(1L, n_obs_pw - abs(lag_at_max_pw))
ci95_pw  <- 1.96 / sqrt(N_eff_pw)
is_sig_pw <- (abs(max_corr_pw) > ci95_pw)

# Save full-range prewhitened CCF plot
png(file.path(out_dir, sprintf("ccf_prewhitened_%s.png", label)), width = 1400, height = 800, res = 150)
plot(ccf_pw_full,
     main = sprintf("Prewhitened CCF: %s vs %s  [report window ±%0.1f h]",
                    v2_col, v1_col, lag_window_samples*dt_min/60),
     xlab = "Lag (samples)", ylab = "CCF")
abline(v = lag_at_max_pw, col = "red", lty = 2)
dev.off()

# ---- 8) Export CSV report ----
report <- data.frame(
  label               = label,
  id_norm             = id_norm,
  window_start_utc    = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
  window_end_utc      = format(end_time,   "%Y-%m-%d %H:%M:%S %Z"),
  vector1_column      = v1_col,
  vector2_column      = v2_col,
  dt_min_per_sample   = dt_min,                                   # MINUTES per sample
  lag_window_hours    = lag_window_samples * dt_min / 60,         # reported window in hours
  # Raw CCF (constrained to window)
  raw_peak_abs_ccf    = as.numeric(max_corr),
  raw_lag_samples     = as.integer(lag_at_max),
  raw_lag_hours       = as.numeric(lag_hours),
  # Prewhitened CCF (constrained to window)
  pw_peak_abs_ccf     = as.numeric(max_corr_pw),
  pw_lag_samples      = as.integer(lag_at_max_pw),
  pw_ci95_band        = as.numeric(ci95_pw),
  pw_significant      = is_sig_pw
)

readr::write_csv(report, file.path(out_dir, sprintf("ccf_report_%s.csv", label)))

# ---- 9) Console summary (optional) ----
cat(
  sprintf("Cadence ≈ %.3f min/sample (%.1f sec)\n", dt_min, dt_min*60),
  "Saved:\n",
  " -", file.path(out_dir, sprintf("ccf_%s.png", label)), "\n",
  " -", file.path(out_dir, sprintf("ccf_prewhitened_%s.png", label)), "\n",
  " -", file.path(out_dir, sprintf("ccf_report_%s.csv", label)), "\n",
  sep = ""
)
cat(sprintf("Raw peak |CCF| = %.3f at lag %d (%.2f h)\n",
            report$raw_peak_abs_ccf, report$raw_lag_samples, report$raw_lag_hours))
cat(sprintf("Prewhitened peak |CCF| = %.3f at lag %d (CI95=±%.3f) => %s\n",
            report$pw_peak_abs_ccf, report$pw_lag_samples, report$pw_ci95_band,
            ifelse(report$pw_significant, "SIGNIFICANT", "not significant")))
