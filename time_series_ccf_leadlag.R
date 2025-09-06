# ---- Parameters (tweak as needed) -----------------------------------------
slope_smooth_minutes <- 5       # smoothing window for slope calc (minutes)
event_z_threshold    <- 1.5     # threshold on |z-score of slope| to call an event
pairing_window_min   <- 60      # max time gap to pair events (minutes)
png_w <- 1400; png_h <- 800; png_res <- 150

# ---- Build smoothed series (reuse standardized series if present) ----------
# Expect plot_df with vector1_z / vector2_z from block 2; otherwise create it.
if (!exists("plot_df") || !"vector1_z" %in% names(plot_df) || !"vector2_z" %in% names(plot_df)) {
  plot_df <- df %>%
    dplyr::transmute(
      t_stamp,
      vector1_raw = .data[[v1_col]],
      vector2_raw = .data[[v2_col]]
    ) %>%
    dplyr::mutate(
      vector1_z = as.numeric(scale(zoo::na.approx(vector1_raw, na.rm = FALSE))),
      vector2_z = as.numeric(scale(zoo::na.approx(vector2_raw, na.rm = FALSE)))
    )
}

# Smoothing window in *samples*
w <- max(3L, round(slope_smooth_minutes / dt_min))
v1_sm <- zoo::rollmean(plot_df$vector1_z, k = w, fill = NA, align = "center")
v2_sm <- zoo::rollmean(plot_df$vector2_z, k = w, fill = NA, align = "center")

# Slopes (per minute): Δvalue / Δtime
v1_slope <- c(NA, diff(v1_sm)) / dt_min
v2_slope <- c(NA, diff(v2_sm)) / dt_min

# Align with timestamps & drop NAs introduced by smoothing/diff
slopedf <- tibble::tibble(
  t_stamp   = plot_df$t_stamp,
  vector1_slope = v1_slope,
  vector2_slope = v2_slope
) %>%
  dplyr::filter(is.finite(vector1_slope), is.finite(vector2_slope))

if (nrow(slopedf) < 20) {
  warning("Too few slope points after smoothing to analyze; consider reducing 'slope_smooth_minutes'.")
}

# ---- (A) CCF on slopes (which series moves first?) ------------------------
# Keep order consistent with earlier blocks: CCF(vector2, vector1)
ccf_slopes <- ccf(as.numeric(slopedf$vector2_slope),
                  as.numeric(slopedf$vector1_slope),
                  lag.max = lag_window_samples, plot = FALSE)

# Constrain report to ±lag_window_samples
in_win <- which(abs(ccf_slopes$lag) <= lag_window_samples)
lag_use <- ccf_slopes$lag[in_win]
acf_use <- ccf_slopes$acf[in_win]

lag_at_max_slope <- lag_use[which.max(abs(acf_use))]
max_corr_slope   <- max(abs(acf_use), na.rm = TRUE)
lag_hours_slope  <- lag_at_max_slope * dt_min / 60

# Quick prewhitened significance on slopes (AR(1) residuals)
x <- as.numeric(slopedf$vector2_slope)
y <- as.numeric(slopedf$vector1_slope)
fit_x <- try(stats::arima(x, order = c(1,0,0)), silent = TRUE)
fit_y <- try(stats::arima(y, order = c(1,0,0)), silent = TRUE)
rx <- if (inherits(fit_x, "try-error")) x else resid(fit_x)
ry <- if (inherits(fit_y, "try-error")) y else resid(fit_y)
ccf_pw_slope <- ccf(rx, ry, lag.max = lag_window_samples, plot = FALSE)

in_win_pw <- which(abs(ccf_pw_slope$lag) <= lag_window_samples)
lag_use_pw <- ccf_pw_slope$lag[in_win_pw]
acf_use_pw <- ccf_pw_slope$acf[in_win_pw]
lag_at_max_slope_pw <- lag_use_pw[which.max(abs(acf_use_pw))]
max_corr_slope_pw   <- max(abs(acf_use_pw), na.rm = TRUE)

n_obs_pw <- min(sum(is.finite(rx)), sum(is.finite(ry)))
N_eff_pw <- max(1L, n_obs_pw - abs(lag_at_max_slope_pw))
ci95_pw  <- 1.96 / sqrt(N_eff_pw)
is_sig_pw_slope <- (abs(max_corr_slope_pw) > ci95_pw)

# Save slope-CCF plot (full range)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(out_dir, sprintf("ccf_slopes_%s.png", label)), width = png_w, height = png_h, res = png_res)
plot(ccf_slopes,
     main = sprintf("Slope CCF (Vector2 vs Vector1) — %s\nWindow: %s → %s  |  Peak lag: %d (%.2f h)",
                    id_norm,
                    format(start_time, "%Y-%m-%d %H:%M %Z"),
                    format(end_time,   "%Y-%m-%d %H:%M %Z"),
                    lag_at_max_slope, lag_hours_slope),
     xlab = "Lag (samples)", ylab = "CCF")
abline(v = lag_at_max_slope, col = "red", lty = 2)
dev.off()

# ---- (B) Event-based lead/lag on slope bursts (fast significance) ---------
# Z-score slopes to standardize event detection
z_v1 <- as.numeric(scale(slopedf$vector1_slope))
z_v2 <- as.numeric(scale(slopedf$vector2_slope))

idx_v1 <- which(abs(z_v1) >= event_z_threshold)
idx_v2 <- which(abs(z_v2) >= event_z_threshold)

t_v1 <- slopedf$t_stamp[idx_v1]
t_v2 <- slopedf$t_stamp[idx_v2]

lead_lags_min <- numeric(0)
if (length(t_v2) >= 1 && length(t_v1) >= 1) {
  # For each Vector2 event, find nearest Vector1 event within pairing_window_min
  for (t2 in t_v2) {
    k <- which.min(abs(as.numeric(difftime(t_v1, t2, units = "mins"))))
    if (length(k) == 1 && is.finite(k)) {
      dmin <- as.numeric(difftime(t_v1[k], t2, units = "mins"))
      if (abs(dmin) <= pairing_window_min) {
        # Positive result => Vector2 leads (event in Vector2 occurs earlier)
        # Negative result => Vector1 leads
        lead_lags_min <- c(lead_lags_min, -dmin)  # negate so positive = Vector2 leads
      }
    }
  }
}

# Fast significance (no bootstrap): t-test & Wilcoxon signed-rank vs 0
t_res <- wilcox_res <- NULL
mean_lead_min <- median_lead_min <- NA_real_
p_t <- p_wilcox <- NA_real_

if (length(lead_lags_min) >= 5) {
  mean_lead_min   <- mean(lead_lags_min, na.rm = TRUE)
  median_lead_min <- median(lead_lags_min, na.rm = TRUE)
  # one-sample t-test against 0
  t_res <- try(stats::t.test(lead_lags_min, mu = 0), silent = TRUE)
  if (!inherits(t_res, "try-error")) p_t <- t_res$p.value
  # Wilcoxon signed-rank (robust to non-normality)
  wilcox_res <- try(stats::wilcox.test(lead_lags_min, mu = 0, exact = FALSE), silent = TRUE)
  if (!inherits(wilcox_res, "try-error")) p_wilcox <- wilcox_res$p.value
}

# ---- Export CSV summary ----------------------------------------------------
summary_csv <- tibble::tibble(
  label                 = label,
  id_norm               = id_norm,
  window_start_utc      = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
  window_end_utc        = format(end_time,   "%Y-%m-%d %H:%M:%S %Z"),
  dt_min_per_sample     = dt_min,
  slope_smooth_minutes  = slope_smooth_minutes,
  event_z_threshold     = event_z_threshold,
  pairing_window_min    = pairing_window_min,
  # Slope CCF
  slope_ccf_peak_abs    = as.numeric(max_corr_slope),
  slope_ccf_lag_samples = as.integer(lag_at_max_slope),
  slope_ccf_lag_hours   = as.numeric(lag_hours_slope),
  slope_ccf_pw_peak_abs = as.numeric(max_corr_slope_pw),
  slope_ccf_pw_ci95     = as.numeric(ci95_pw),
  slope_ccf_pw_sig      = is_sig_pw_slope,
  # Event lead/lag summary
  n_vector2_events      = length(idx_v2),
  n_vector1_events      = length(idx_v1),
  n_pairs               = length(lead_lags_min),
  mean_lead_minutes     = mean_lead_min,     # + => Vector2 leads
  median_lead_minutes   = median_lead_min,   # + => Vector2 leads
  p_t_mean_gt0          = p_t,               # t-test vs 0
  p_wilcox_median_gt0   = p_wilcox           # Wilcoxon signed-rank vs 0
)

readr::write_csv(
  summary_csv,
  file.path(out_dir, sprintf("slope_leadlag_report_%s.csv", label))
)

# ---- Optional: quick overlay plot of slopes (z-scored) ---------------------
p_slopes <- slopedf %>%
  dplyr::mutate(
    z_v1 = z_v1,
    z_v2 = z_v2
  ) %>%
  dplyr::select(t_stamp, z_v1, z_v2) %>%
  tidyr::pivot_longer(-t_stamp, names_to = "series", values_to = "z") %>%
  dplyr::mutate(series = dplyr::recode(
    series,
    z_v1 = sprintf("Vector1 slope (z): %s", v1_col),
    z_v2 = sprintf("Vector2 slope (z): %s", v2_col)
  )) %>%
  ggplot2::ggplot(ggplot2::aes(t_stamp, z, color = series)) +
  ggplot2::geom_line(linewidth = 0.5, alpha = 0.9, na.rm = TRUE) +
  ggplot2::geom_hline(yintercept = c(-event_z_threshold, event_z_threshold), linetype = 2, alpha = 0.5) +
  ggplot2::labs(
    title = sprintf("%s — Slope overlay (z-scored)\nWindow: %s → %s",
                    id_norm,
                    format(start_time, "%Y-%m-%d %H:%M %Z"),
                    format(end_time,   "%Y-%m-%d %H:%M %Z")),
    x = "Time", y = "Slope (z)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "top")

ggplot2::ggsave(file.path(out_dir, sprintf("slopes_overlay_%s.png", label)),
                p_slopes, width = 12, height = 6, dpi = 150)

# ---- Console summary -------------------------------------------------------
cat(
  sprintf("Slope CCF: peak |CCF|=%.3f at lag %d (%.2f h)  |  prewhitened |CCF|=%.3f (CI95=±%.3f) => %s\n",
          max_corr_slope, lag_at_max_slope, lag_hours_slope,
          max_corr_slope_pw, ci95_pw, ifelse(is_sig_pw_slope, "SIGNIFICANT", "not significant"))
)
if (length(lead_lags_min) >= 5) {
  cat(sprintf("Event-based lead/lag (Vector2 leads +): mean=%.2f min, median=%.2f min, t p=%.4f, wilcox p=%.4f (n=%d pairs)\n",
              mean_lead_min, median_lead_min, p_t, p_wilcox, length(lead_lags_min)))
} else {
  cat(sprintf("Event-based lead/lag: insufficient paired events (n_pairs=%d). Try lowering 'event_z_threshold' or increasing 'pairing_window_min'.\n",
              length(lead_lags_min)))
}
