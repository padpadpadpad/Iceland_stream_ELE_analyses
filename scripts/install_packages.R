# Install all necessary packages for other scripts
# eco_level_metab_analysis.R
# eco_level_metab_process.R
# pop_level_metab.R

# 1. Install packages from CRAN
# package list from CRAN
pkgs <- c('tidyr',
          'dplyr',
          'ggplot2',
          'nlme',
          'lme4',
          'grid',
          'gridExtra',
          'MuMIn',
          'magrittr',
          'broom',
          'lsmeans',
          'scales',
          'caTools',
          'chron',
          'lubridate',
          'mise')

# check if they are already installed
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

# install packages
if(length(new_pkgs)) install.packages(new_pkgs)

# 2. install packages from GitHub ####

# install devtools
install.packages('devtools') 

# install packages from GitHub
devtools::install_github('padpadpadpad/TeamPhytoplankton',
                         'padpadpadpad/nlsLoop')



