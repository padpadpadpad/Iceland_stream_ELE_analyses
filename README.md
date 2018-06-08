# Analyses of population- and ecosystem-level stream metabolism as published in:

_Padfield et al. (2017) Metabolic compensation constrains the temperature dependence of gross primary production. Ecology Letters, 20.10: 1250-1260_

DOI of paper:

[doi: 10.1111/ele.12820](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12820)

DOI of repository and dataset:

[![DOI](https://zenodo.org/badge/99806728.svg)](https://zenodo.org/badge/latestdoi/99806728)

### Outline

This repository contains the final datasets, analyses and figures of the above-mentioned paper. It can recreate the population- and ecosystem-level analyses and the plots from the main text (Figure 2 and Figure 3).

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio.
- Firstly run the script `install_packages.R` which will install all the packages necessary to run the analyses (at least one is downloaded from GitHub)
- `pop_level_metab.R` contains the analysis of the metabolic rates of population-level metabolism.
- `eco_level_metab_process.R` contains the code to process rates of oxygen concentration into estimates of stream metabolism. This script follows the methods in the Supplementary Information of the paper. This script does not have to be run for `eco_level_metab_analysis.R` to work.
- `eco_level_metab_analysis.R` contains the analysis of ecosystem-level GPP.
- All of the data needed to run the analyses are stored in `data/`
- All figures produced by the analysis are eventually saved in `figures/` and are labelled as they are in the main text.

__All analses are done in R version 3.4.1, on macOS Sierra 10.12.3. I am unsure whether some older version of R will support all of the packages and whether the analyses will run exactly the same.__

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Iceland_stream_ELE_analyses/issues) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### NB

I am unsure whether every single package is needed for the scripts but they have been streamlined from previous analyses so there might be some unnecessary loading of packages. Apologies for any illegible coding! It may also be the case that I have missed out a package, in which case they can be installed in the same way as in `install_packages.R`.
