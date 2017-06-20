#!/bin/bash
echo $REXE
if [ -z ${REXE}]
then
REXE=Rscript
fi
echo $REXE
## set directory to the project
cd path/to/package/annual-predictions-paper
## create the project and load in the packages required to run forecasts (this step may take a while)
$REXE -e "install.packages('packrat/src/packrat/packrat_0.4.8-1.tar.gz', repos=NULL, type='source');install.packages('packrat/src/stringi/stringi_1.1.3.tar.gz', repos=NULL, type='source');packrat::init()"
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
