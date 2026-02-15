# Replication Materials: Robust Bayesian Spatio-Temporal Modelling of Conflict Data

This repository contains the data and `R` code required to reproduce the results, tables, and figures presented in the manuscript:

> Title: Robust Bayesian Spatio-Temporal Modelling of Conflict Data

## Repository Structure

The project follows a reproducible, project-oriented workflow.

* `data/`: Contains the processed datasets required for the analysis.
    * `FinalDATASET.rds`: The main event/fatality dataset.
    * `cabo_delgado_boundary.rds`: Shapefiles for the Cabo Delgado districts.
    
* `scripts/`: Contains the R analysis scripts.
    * `Rob_BAMLSS_ST.R`: The primary script that runs the PLE selection, fits the BAMLSS models (ZIPIG, etc.), and generates diagnostics.
    
* `figures/`: Output directory where generated plots will be saved.

* `Robust_conflict_data_analysis.Rproj`: RStudio Project file. This is the entry point for the analysis.

## Software and Dependencies

The analysis was performed using `R` (version 4.3.2). 
To ensure reproducibility, this project uses the `here` package for relative file paths, eliminating the need for absolute paths (e.g., `setwd()`).

**Required Packages:**
The script utilizes `pacman` to handle package management. Key dependencies include `bamlss`, `gamlss.dist`, `sf`, `spdep`, `ggplot2`, and `here`.

## Instructions for Reproduction

To reproduce the analysis, please follow these steps strictly to ensure file paths work correctly:

1.  Download this repository (Click "Code" -> "Download ZIP") and unzip it.
2.  Open the file `Robust_conflict_data_analysis.Rproj` in RStudio. 
    * *Critical:* Opening the project file sets the working directory to the project root automatically. Do not use `setwd()`.
3.  Open the script `scripts/Rob_BAMLSS_ST.R`.
4.  Run the code. 
    * *Note on Computation Time:* The full analysis presented in the paper utilizes 10,000 MCMC iterations (Phase 2), which is computationally intensive and may take several hours on a standard laptop. For a faster replication test, you may reduce `ni` to 1,000 and `nb` to 200 in the script; this will produce qualitatively consistent estimates in a fraction of the time.

## Data Source

The datasets constructed for this analysis integrate information from the following primary sources:

* Conflict Data: Armed Conflict Location & Event Data Project (ACLED).

* Demographics: United States Census Bureau International Data Base.

* Environmental Data: National Oceanic and Atmospheric Administration (NOAA) and Climate Hazards Center.

* Boundaries: Shapefiles were obtained from GADM.

## Contact

For questions regarding the code or data, please contact the corresponding author via the journal submission system to maintain the integrity of the double-blind review process.
