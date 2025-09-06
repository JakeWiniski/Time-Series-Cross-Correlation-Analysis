# ---- Ensure output dir exists (reuse out_dir + label from previous block) ----
if (!exists("out_dir")) out_dir <- "ccf_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!exists("label")) {
  if (exists("ts_file_path")) {
    label <- tools::file_path_sans_ext(basename(ts_file_path))
  } else {
    label <- "timeseries"
  }
}

# ---- Build plotting frame with raw + standardized versions ----
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

# ---- Overlay of standardized series ----
overlay_title <- if (exists("lag_at_max") && exists("lag_hours")) {
  sprintf("%s — Vector1 vs Vector2 (z-scored)\nWindow: %s → %s  |  Peak lag = %d samples (%.2f h)",
          id_norm,
          format(start_time, "%Y-%m-%d %H:%M %Z"),
          format(end_time,   "%Y-%m-%d %H:%M %Z"),
          lag_at_max, lag_hours)
} else {
  sprintf("%s — Vector1 vs Vector2 (z-scored)\nWindow: %s → %s",
          id_norm,
          format(start_time, "%Y-%m-%d %H:%M %Z"),
          format(end_time,   "%Y-%m-%d %H:%M %Z"))
}

p_overlay <- plot_df %>%
  dplyr::select(t_stamp, vector1_z, vector2_z) %>%
  tidyr::pivot_longer(-t_stamp, names_to = "series", values_to = "z") %>%
  ggplot2::ggplot(ggplot2::aes(t_stamp, z, color = series)) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.95, na.rm = TRUE) +
  ggplot2::scale_color_manual(
    values = c(vector1_z = "#1f77b4", vector2_z = "#ff7f0e"),
    labels = c(vector1_z = sprintf("Vector1 (z): %s", v1_col),
               vector2_z = sprintf("Vector2 (z): %s", v2_col))
  ) +
  ggplot2::labs(title = overlay_title, x = "Time", y = "Standardized value (z)", color = NULL) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "top")

# ---- Raw values in facets ----
p_raw <- plot_df %>%
  dplyr::select(t_stamp, vector1_raw, vector2_raw) %>%
  tidyr::pivot_longer(-t_stamp, names_to = "series", values_to = "value") %>%
  dplyr::mutate(series = dplyr::recode(
    series,
    vector1_raw = sprintf("Vector1 (raw): %s", v1_col),
    vector2_raw = sprintf("Vector2 (raw): %s", v2_col)
  )) %>%
  ggplot2::ggplot(ggplot2::aes(t_stamp, value)) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.95, na.rm = TRUE) +
  ggplot2::facet_wrap(~ series, ncol = 1, scales = "free_y") +
  ggplot2::labs(
    title = sprintf("%s — Raw time series\nWindow: %s → %s",
                    id_norm,
                    format(start_time, "%Y-%m-%d %H:%M %Z"),
                    format(end_time,   "%Y-%m-%d %H:%M %Z")),
    x = "Time", y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12)

# ---- Lag-aligned overlay (Vector2 shifted by lag_at_max) ----
p_lag_aligned <- NULL
if (exists("lag_at_max") && lag_at_max != 0) {
  plot_df_lag <- plot_df %>%
    dplyr::mutate(
      vector2_z_shift = if (lag_at_max > 0) {
        dplyr::lag(vector2_z, n = lag_at_max)
      } else {
        dplyr::lead(vector2_z, n = abs(lag_at_max))
      }
    )

  p_lag_aligned <- plot_df_lag %>%
    dplyr::select(t_stamp, vector1_z, vector2_z_shift) %>%
    tidyr::pivot_longer(-t_stamp, names_to = "series", values_to = "z") %>%
    dplyr::mutate(series = dplyr::recode(
      series,
      vector1_z = sprintf("Vector1 (z): %s", v1_col),
      vector2_z_shift = sprintf("Vector2 (z) shifted %d samples (%.2f h): %s",
                                lag_at_max, lag_hours, v2_col)
    )) %>%
    ggplot2::ggplot(ggplot2::aes(t_stamp, z, color = series)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.95, na.rm = TRUE) +
    ggplot2::labs(
      title = sprintf("%s — Lag-aligned overlay", id_norm),
      x = "Time", y = "Standardized value (z)", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}

# ---- Save PNGs ----
ggplot2::ggsave(filename = file.path(out_dir, sprintf("overlay_%s.png", label)),
                plot = p_overlay, width = 12, height = 6, dpi = 150)
ggplot2::ggsave(filename = file.path(out_dir, sprintf("raw_%s.png", label)),
                plot = p_raw, width = 12, height = 7, dpi = 150)
if (!is.null(p_lag_aligned)) {
  ggplot2::ggsave(filename = file.path(out_dir, sprintf("lag_aligned_%s.png", label)),
                  plot = p_lag_aligned, width = 12, height = 6, dpi = 150)
}

# Optionally also print to the Rmd output:
p_overlay
p_raw
if (!is.null(p_lag_aligned)) p_lag_aligned
