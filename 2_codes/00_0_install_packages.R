# This Script -------------------------------------------------------------
# Author: Laura Zwyssig
# Goal: Install/update all necessary packages
# Last edited: 03.06.2018

# Install/update all packages below
# The specific version numbers used for this project
# can be found in the README.html file
dependencies <- c('broom',
                  'car',
                  'ggplot2',
                  'haven',
                  'labelled',
                  'lmtest',
                  'multcomp',
                  'plm',
                  'lfe',
                  'plyr',
                  'readr',
                  'readxl',
                  'reshape2',
                  'sandwich',
                  'scales',
                  'stargazer',
                  'stringr',
                  'tidyverse',
                  'rdrobust',
                  'dplyr',
                  'RStata',
                  'GGally',
                  'VGAM')

install.packages(dependencies)
