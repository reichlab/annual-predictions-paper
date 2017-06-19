## annual-predictions workflow
## Stephen Lauer, April 2017

## put together forecasting data from raw data
source("R/make-forecasting-data.R")
## conduct cross validation, model selection, and test phase forecasts
source("R/make-forecasts.R")

## compile manuscript/Lauer-annual-DHF-figures.Rnw
## compile manuscript/Lauer-annual-DHF-supplement.Rnw
## compile manuscript/Lauer-annual-DHF-manuscript.Rnw