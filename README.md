# annual-predictions-paper
The code and data required to make annual dengue hemorrhagic fever (DHF) incidence forecasts for Thailand.

In order to re-run the experiment, please run the shell script inst/annual-forecasts-workflow.sh, which performs the following:

1) creates an R `packrat` project, which installs all of the packages used to originally run the analysis
2) runs `make-forecasting-data.R`, which puts together the forecasting data from the raw data
3) runs `make-forecasts.R` which conducts cross validation, model selection, and makes out-of-sample forecasts
4) runs `Lauer-annual-DHF-figures.Rnw` which makes the figures for the manuscript
5) runs `Lauer-annual-DHF-supplement.Rnw` which makes the supplement
6) runs `Lauer-annual-DHF-manuscript.Rnw` which makes the manuscript

If any of this process fails to work, please make an issue and I will try to help!
