# AI coverage and information inequality in the UK press: A natural language processing approach to media attention economics
Data and R code for analysing changes in AI-related coverage in
The Guardian and The Sun around the public release of ChatGPT on
30 November 2022.

## Contents

- `weekly_series_full.csv`: weekly topic shares and article counts.
- `stationarity-test.R`: ADF, Phillips–Perron, and KPSS tests.
- `timeseris.R`: interrupted time-series (ITS) and ARIMAX estimation,
  diagnostics, and figures.

The data cover July 2016–October 2024 and contain seven topics
for each newspaper.

## Requirements

R packages: `tseries`, `urca`, `forecast`, `lmtest`, `sandwich`, and `zoo`.

## Usage

1. Download the R scripts and CSV files into the same directory.
2. Set that directory as the R working directory.
3. Run `stationarity-test.R`, followed by `timeseris.R`.

Ensure that the topic labels in each script match its input CSV.

## Scope

This repository contains only the time-series data and analysis.
Article collection, text preprocessing, and topic-modelling code
are not included.

## Citation

Please cite the associated article when using these materials.
