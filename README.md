# annual-predictions-paper
In order to re-run this experiment, please follow these steps:
1) download this repository by pressing the green button in the upper right-hand corner
2) unzip the downloaded folder
3) replace `path/to/package` in line 9 of the shell script `inst/annual-forecasts-workflow.sh` with the path to the unzipped folder
4) run the shell script (on a Mac in the Terminal write: `sh path/to/package/annual-predictions-paper-master/inst/annual-forecasts-workflow.sh`; the author does not know the PC equivalent of this command)

This shell script:
1) creates an R `packrat` project, which installs all of the packages used to originally run the analysis (note: this step may take 10 minutes or longer, sometimes due to a long installation of the `stringi` package)
2) runs `make-forecasting-data.R`, which puts together the forecasting data from the raw data
3) runs `make-forecasts.R` which conducts cross validation, model selection, and makes out-of-sample forecasts
4) runs `Lauer-annual-DHF-figures.Rnw` which makes the figures for the manuscript
5) runs `Lauer-annual-DHF-supplement.Rnw` which makes the supplement
6) runs `Lauer-annual-DHF-manuscript.Rnw` which makes the manuscript

If any of this process fails to work, please make an issue and I will try to help!
