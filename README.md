# Acclimated_gs_optimization_models

## Description

Welcome to the GitHub repository for the study: "Incorporating Photosynthetic Acclimation Improves Stomatal Optimisation Models" by Victor Flo, Jaideep Joshi, Manon Sabot, David Sandoval, and Iain Colin Prentice (2024).

### Overview
This repository includes a collection of R scripts and complementary datasets.

### Abstract
Stomatal opening in plant leaves is regulated through a balance of carbon and water exchange under different environmental conditions. Accurate estimation of stomatal regulation is crucial for understanding how plants respond to changing environmental conditions, particularly under climate change. A new generation of optimality-based modelling schemes determines instantaneous stomatal responses from a balance of trade-offs between carbon gains and hydraulic costs, but most such schemes do not account for biochemical acclimation in response to drought. Here, we compare the performance of six instantaneous stomatal optimization models with and without accounting for photosynthetic acclimation. Using experimental data from 37 plant species, we found that accounting for photosynthetic acclimation improves the prediction of carbon assimilation in a majority of the tested models.  Photosynthetic acclimation contributed significantly to the reduction of photosynthesis under drought conditions in all tested models. Drought effects on photosynthesis could not accurately be explained by the hydraulic impairment functions embedded in the stomatal models alone, indicating that photosynthetic acclimation must be considered to improve estimates of carbon assimilation during drought.

### Repository Contents
- **R Scripts:** Each script is crafted to replicate the analyses and model simulations as presented in our study.
- **Data Sets:** Data set used in our analysis: "DATA/drying_experiments_meta-analysis_extended.csv"


## Getting Started

This section guides you through the initial steps to get the repository up and running on your machine.

### Cloning the project in RStudio

#### Step 1: Install Git (if not already installed)
Ensure Git is installed on your system. RStudio uses Git for version control and to clone repositories.

#### Step 2: Configure Git in RStudio
1. Open RStudio.
2. Go to `Tools` > `Global Options` > `Git/SVN`. Ensure that the path to the Git executable is correctly set, and configure your user name and email for Git.

#### Step 3: Start a New Project from Version Control
1. In RStudio, go to `File` > `New Project`.
2. Choose `Version Control` and then `Git`.
3. In the “Repository URL” field, paste the URL of the repository you wish to clone: `https://github.com/vflo/Acclimated_gs_optimization_models.git`.
4. Fill in the other fields such as the Project directory name and the location on your computer where you'd like the project to be saved.
5. Click `Create Project`.

RStudio will then clone the repository, and the project will be opened with the contents of the repository ready for use.

#### Installing Dependencies
The project relies on several R packages. You can install them using the following R command:

```R
install.packages(c("dplyr", "purrr", "tidyr", "ggplot2", "gridExtra", "scales", "zoo", "stringr", "DEoptim", "optimr", "lmerTest", "ggpubr", "ggalt", "grid", "ggpattern", "gridExtra", "MASS", "psych", "effects", "emmeans", "rstatix"))
```

Ensure you have R installed on your system. Run this command in your R environment to set up all required dependencies. This step is crucial for the proper execution of the scripts in the repository.


### Executing the Results

To achieve the desired results, follow these steps:

#### Calibration Phase
Begin with the calibration scripts. These scripts are crucial for setting up the models accurately:
1. **calibrate_kmaxww_alpha.R**: This script corresponds to calibrations 1 and 2 of the study. It's the initial step in the calibration process.
2. **calibrate_kmaxww_alpha_1_3.R**: Similar to the first script but utilizes P50OX. An alternative version of calibrations 1 and 2.
3. **calibrate_kmax_alpha.R**: Corresponds to calibration 3. Be aware that this process is computationally intensive and may take several days on a standard computer.

#### Simulation Phase
After completing the calibrations, proceed with the simulation scripts:
1. **simulate_kmaxww_alpha.R**
2. **simulate_kmaxww_alpha_1_3.R**
3. **simulate_kmax_alpha.R**

These scripts are designed to model g[s] and A[net] for the different calibration scenarios, providing the framework for generating the data needed for analysis.

#### Analysis
Once the simulations are complete, it's time to compile and analyze the data:
- **analysis.R**: This script will process the simulation outputs, synthesizing them into coherent results. Follow the sections of the script step by step.

#### Note on Pre-runned Simulations
For convenience, this repository includes pre-run simulation outputs. This means you can directly proceed to the analysis phase without performing the calibrations and simulations, should you prefer. This option is ideal for those looking to quickly validate results or for environments where computational resources are limited.


## Authors

**Victor Flo**
- Email: [vflosierra@gmail.com](mailto:vflosierra@gmail.com)

Feel free to reach out for any queries or collaborations related to this project. Your contributions and feedback are greatly appreciated.


## License
This project is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
