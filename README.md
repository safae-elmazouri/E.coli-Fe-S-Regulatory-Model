# E. coli Fe–S Regulatory Model

## Overview

This repository contains the logical regulatory model, simulation outputs, and R scripts developed to investigate Fe–S cluster biogenesis and stress response regulation in *Escherichia coli*.

The model integrates regulatory mechanisms involved in iron homeostasis, oxidative and nitrosative stress responses, and adaptation to oxygen availability.

It focuses on the regulation of the Isc and Suf Fe–S cluster biogenesis systems under different environmental conditions, including iron availability, oxygen, hydrogen peroxide, and nitric oxide.

The repository provides the model and computational resources used to analyze its dynamical behavior and generate the associated figures.

## Repository structure

- **Model/**: Final logical regulatory model of Fe–S cluster biogenesis and stress response in *E. coli* (.zginml format).
- **Data/**: Simulation outputs and datasets used for the analyses.
- **R_scripts/**: R scripts used to analyze simulation results and generate figures.
  - **01_WT_and_KO/**: Wild-type and knockout mutant analyses.
  - **02_Attractor_Analysis/**: Analysis of attractor sizes and dynamical behavior.
  - **03_Mutant_heatmaps/**: Heatmap generation for mutant simulations.
 
## Requirements

The analyses and figures were generated using R.

The R scripts use the following packages:

- ggplot2
- dplyr
- tidyr
- reshape2
- gridExtra
- patchwork

The logical model is provided in `.zginml` format, compatible with GINsim.

## Reproducing the analyses

The repository provides the logical model, simulation outputs, and R scripts used in this study.

1. Download or clone the repository.
2. Open the `R_scripts/` directory.
3. Install the required R packages.
4. Use the datasets provided in `Data/` to run the corresponding R scripts.

The scripts are organized into three categories:

- **01_WT_and_KO/**: Analysis and visualization of wild-type and knockout simulations.
- **02_Attractor_Analysis/**: Analysis of attractor sizes and dynamical behavior.
- **03_Mutant_heatmaps/**: Generation of heatmaps for mutant simulations.
