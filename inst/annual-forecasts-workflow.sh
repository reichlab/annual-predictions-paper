#!/bin/bash
echo $REXE
if [ -z ${REXE}]
then
REXE=Rscript
fi
echo $REXE
## set directory to the project
cd path/to/package
## create the project and load in the packages required to run forecasts (this step may take 10+ minutes, especially to install stringi)
$REXE -e "install.packages('packrat/src/packrat/packrat_0.4.8-1.tar.gz', repos=NULL, type='source');packrat::init()"
## put together forecasting data from raw data
$REXE R/make-forecasting-data.R > data/make-forecasting-data-output.Rout
## conduct cross validation, model selection, and test phase forecasts
$REXE R/make-forecasts.R > data/make-forecasts.Rout
## set workspace to manuscript folder
cd manuscript
## make the figures for the manuscript
$REXE -e "library(knitr);knit2pdf('Lauer-annual-DHF-figures.Rnw')"
## make the supplement
$REXE -e "library(knitr);knit2pdf('Lauer-annual-DHF-supplement.Rnw')"
## make the manuscript
$REXE -e "library(knitr);knit2pdf('Lauer-annual-DHF-manuscript.Rnw')"
