#!/bin/bash
echo $REXE
if [ -z ${REXE}]
then
REXE=Rscript
fi
echo $REXE
## create the project
$REXE -e "install.packages('devtools');devtools::create('path/to/package/annual-forecasts')"
## set directory to the project
cd path/to/package/annual-forecasts
## load in the library required to run forecasts
$REXE -e "install.packages('packrat');packrat::init()"
## put together forecasting data from raw data
$REXE R/make-forecasting-data.R > data/make-forecasting-data-output.Rout
## conduct cross validation, model selection, and test phase forecasts
$REXE R/make-forecasts.R > data/make-forecasts.Rout
## make the figures for the manuscript
$REXE -e "library(knitr);setwd('manuscript');knit2pdf('Lauer-annual-DHF-figures.Rnw')"
## make the supplement
$REXE -e "library(knitr);setwd('manuscript');knit2pdf('Lauer-annual-DHF-supplement.Rnw')"
## make the manuscript
$REXE -e "library(knitr);setwd('manuscript');knit2pdf('Lauer-annual-DHF-manuscript.Rnw')"
