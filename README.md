# Cross-Correlation Time Series Workflow

This repository contains an R workflow for **lead–lag analysis between two time series vectors**.  
It uses cross-correlation functions (CCF), visualization, and slope-based event detection.  

The workflow is organized into **three modular scripts** so each step can be reused or adapted for different datasets.

---

## Overview

This workflow helps answer questions such as:

- *Are changes in one series leading or lagging changes in another?*  
- *What is the strength of that relationship, and is it statistically significant?*  
- *How does the lead–lag structure change when looking at slopes (rates of change) rather than raw values?*  

It also generates **publication-quality plots** and **summary CSV reports** for downstream analysis.

---

## Workflow Structure

### 1. Analysis: Cross-Correlation Calculation  
**File:** `time_series_ccf_analysis.R`

- Loads a time-series dataset and a reference file defining analysis windows (start/end times).  
- Extracts two vectors of interest (configured via regex patterns).  
- Standardizes, interpolates, and windows the data.  
- Computes raw and prewhitened cross-correlations within a user-defined lag window.  
- **Outputs:**  
  - PNG plots of the CCF (raw + prewhitened).  
  - CSV report of peak correlation, lag, and significance.  

---

### 2. Visualization: Overlay and Raw Series  
**File:** `time_series_ccf_plots.R`

- Builds a plotting frame with both raw and standardized values.  
- Generates three plots:  
  - **Overlay plot** of z-scored series.  
  - **Raw values plot** faceted by series.  
  - **Lag-aligned overlay** (one series shifted by the detected lag).  
- **Outputs:**  
  - PNG files for each plot.  

---

### 3. Lead–Lag on Slopes: Event-Based Analysis  
**File:** `time_series_ccf_leadlag.R`

- Smooths series and calculates slopes (rate of change).  
- Computes CCF on slopes (raw + prewhitened).  
- Detects “events” where slope magnitude exceeds a z-score threshold.  
- Pairs events across series within a configurable time window to estimate who leads.  
- Performs quick significance tests (t-test and Wilcoxon signed-rank).  
- **Outputs:**  
  - PNG of slope CCF.  
  - PNG of z-scored slope overlay with event thresholds.  
  - CSV summary report of slope-based lead/lag metrics.  

---

## Outputs

Each script writes to the `ccf_outputs/` directory (created automatically).  
Typical outputs include:

ccf_<label>.png
ccf_prewhitened_<label>.png
overlay_<label>.png
raw_<label>.png
lag_aligned_<label>.png
ccf_slopes_<label>.png
slopes_overlay_<label>.png
ccf_report_<label>.csv
slope_leadlag_report_<label>.csv


Where `<label>` is derived from the input file name.

---

## Customization

- **Patterns:** adjust `vector1_pattern` and `vector2_pattern` to match your dataset’s column names.  
- **Windows:** define start/end windows in the reference CSV (flexible date parsing).  
- **Parameters:** tweak lag window size (`max_lag_hours`), smoothing length, event thresholds, and output paths.  

---

## Dependencies

This workflow uses several R packages:

- **tidyverse** (dplyr, ggplot2, tidyr, readr, stringr)  
- **lubridate**  
- **zoo**  
- **stats** (base R)  

Install them with:
install.packages(c("tidyverse", "lubridate", "zoo"))

## Suggested Usage

Rscript time_series_ccf_analysis.R
Rscript time_series_ccf_plots.R
Rscript time_series_ccf_leadlag.R

