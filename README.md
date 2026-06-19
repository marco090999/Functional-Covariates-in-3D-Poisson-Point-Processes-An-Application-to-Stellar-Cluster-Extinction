# Functional Spatial Covariates in Three-Dimensional Poisson Point Processes

This repository contains the code and processed data used for the study:

> **Functional Spatial Covariates in Three-Dimensional Poisson Point Processes: An Application to Stellar Cluster Extinction**

The project develops a modelling framework for incorporating function-valued spatial covariates into the first-order intensity of three-dimensional Poisson point processes. The methodology combines Functional Principal Component Analysis (FPCA), spatial completion of FPCA score fields, and log-linear Poisson point process regression.

The empirical application investigates the association between the spatial intensity of young stellar clusters and the within-cluster distribution of stellar extinction.

## Repository structure

```text
.
├── R/
│   ├── 01_simulation_study.R
│   ├── 02_stellar_cluster_application.R
│   └── 03_fixedK_consistency.R
│
├── data/
│   └── processed/
│       ├── stars_with_extinction.rds
│       └── simulation_summary_tables.RData
│
├── .gitignore
├── LICENSE
└── README.md
```

## Code

### `R/01_simulation_study.R`

This script contains the code for the main simulation study. It evaluates the finite-sample performance of the proposed functional Poisson point process model under different data-generating scenarios.

The simulation study considers:

* recovery of the functional coefficient;
* recovery of the Poisson intensity;
* comparison with scalar-summary and reduced alternatives;
* sensitivity to the number of retained FPCA components;
* sensitivity to spatial score-field completion;
* robustness under alternative functional and spatial configurations.

### `R/02_stellar_cluster_application.R`

This script reproduces the application to young stellar clusters in the Solar neighbourhood.

The workflow includes:

* preparation of the star-level extinction data;
* construction of cluster-specific extinction quantile curves;
* Functional Principal Component Analysis of the extinction curves;
* spatial interpolation of FPCA score fields;
* fitting of the three-dimensional functional Poisson point process model;
* comparison with scalar-summary and baseline models;
* reconstruction of the estimated functional coefficient;
* generation of the main application figures, tables, and diagnostic analyses.

### `R/03_fixedK_consistency.R`

This script contains the simulation code supporting the fixed-(K) consistency analysis of the proposed truncated estimator.

It evaluates the empirical behaviour of the estimator as the expected number of events increases and compares the oracle score-field specification with the interpolated score-field specification.

## Data

### `data/processed/stars_with_extinction.rds`

This file contains the minimal processed star-level dataset required to reproduce the stellar-cluster application.

The dataset includes the following variables:

| Variable      | Description                                                                 |
| ------------- | --------------------------------------------------------------------------- |
| `ID_CL_PAPER` | Cluster identifier                                                          |
| `X`           | Heliocentric Cartesian coordinate along the first spatial axis, in parsecs  |
| `Y`           | Heliocentric Cartesian coordinate along the second spatial axis, in parsecs |
| `Z`           | Heliocentric Cartesian coordinate along the third spatial axis, in parsecs  |
| `A0_ML`       | Estimated stellar extinction                                                |

The coordinates stored as `X`, `Y`, and `Z` correspond to the original variables `X1`, `Y1`, and `Z1` in the source star-level dataset.

### `data/processed/simulation_summary_tables.RData`

This file contains the tabular summaries of the simulation results used to produce the tables and figures reported in the manuscript.

It includes only data-frame-like objects and excludes functions, fitted model objects, intermediate workspace objects, and other session-specific quantities.

## Requirements

The analyses were developed in R. The main packages used in the project include:

```r
mgcv
dplyr
tidyr
ggplot2
GET
plot3D
```

Additional packages may be required depending on the selected workflow and the functions used in the scripts.

## Reproducibility

The scripts are designed to be run from the root directory of the repository.

A typical workflow is:

```r
source("R/01_simulation_study.R")
source("R/03_fixedK_consistency.R")
source("R/02_stellar_cluster_application.R")
```

The simulation workflows can be computationally demanding. The file `simulation_summary_tables.RData` is included to provide the tabular results reported in the manuscript without requiring users to rerun all simulations immediately.

## Output

When executed, the scripts generate model estimates, diagnostic quantities, tables, and figures corresponding to the simulation study, the fixed-(K) consistency analysis, and the stellar-cluster application.

Generated output files are not necessarily tracked in the repository and may need to be stored locally depending on their size.

## Licence

This repository is distributed under the licence specified in the `LICENSE` file.

## Authors

Marco Tarantino, Nicoletta D'Angelo, Loredana Prisinzano, Radu S. Stoica, and Giada Adelfio.
