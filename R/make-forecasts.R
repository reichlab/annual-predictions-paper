library(dplyr)
library(dengueThailand)
library(mgcv)
# devtools::install_github("krlmlr/here")
library(doMC)
registerDoMC(20)

set.seed(1)
VERBOSE <- TRUE
TIMES <- TRUE
FIRST_TEST_YEAR <- 2010
LAST_TEST_YEAR <- 2014
MAX_KNOTS <- 8

setwd(here:here())
source("R/annual-forecasting-utilities.R")

## read in latest data for prediction
most_recent_dat <- find_recent_file(name_start = "annual-pred-data",
                                    path = "data/")
dengue_dat <- read.csv(most_recent_dat)[,-1]

## create a cross-validation set for finding the best model
years <- unique(dengue_dat$year)
cv_years <- years[which(years < FIRST_TEST_YEAR)]
cv_dat <- filter(dengue_dat, year %in% cv_years)

## empty data to be filled in later
cv_train <- c()
best_models <- c()
test_preds <- c()

## label each variable if they are available to group 1, 2, or 3
var_set <- data_frame(var = names(dengue_dat),
                      set = c(rep(1, 17), rep(2, 10), rep(3, 11)))

##### Run cross validation on all single covariate GAMs
#####

## create a data frame of all combinations of variables, years, and knots
cov_df <- expand.grid(var = var_set$var[4:nrow(var_set)],
                      knots = 3:MAX_KNOTS,
                      year = cv_years) %>%
    left_join(var_set, by = "var")
## if a prior cross validation run is read in, remove duplicates
dup_rows <- which(cov_df$var %in% cv_train$model_covs &
                      cov_df$year %in% cv_train$year &
                      cov_df$knots %in% cv_train$knots)
cov_df <- if(length(dup_rows)>0) cov_df[-c(dup_rows),] else cov_df

## predict each row in the data frame in parallel
if(nrow(cov_df)>0){
    tic <- Sys.time()
    each_train <- foreach(i=1:nrow(cov_df), .combine = rbind) %dopar% {
        if(VERBOSE)
            print(paste("simulation", i, "of", nrow(cov_df)))
        predict_single_season(x = cv_dat,
                              vars = cov_df$var[i],
                              test_year = cov_df$year[i],
                              loso_cv = TRUE,
                              knots = cov_df$knots[i],
                              sets=cov_df$set[i],
                              distribution = FALSE)
    }
    cv_train <- bind_rows(cv_train, each_train)
    toc <- Sys.time()
}

if(TIMES) toc-tic
if(nrow(cov_df)>0)
    saveRDS(cv_train, paste0("data/cv-train-", Sys.Date(), ".rds"))

## create a table of all the best single covariate models
## take the model with the least number of knots within one standard deviation of the best model for each covariate
table_dat <- cv_train %>%
    filter(num_cov==1) %>%
    group_by(model_covs, year, knots) %>%
    summarise(year_mae = mean(abs(log(obs_counts)-log(pred_counts)))) %>%
    group_by(model_covs, knots) %>%
    summarise(cv_mae = mean(year_mae),
              cv_se = sd(year_mae)) %>%
    group_by(model_covs) %>%
    filter(knots <= knots[which.min(cv_mae)],
           cv_mae <= cv_mae[which.min(cv_mae)]+cv_se[which.min(cv_mae)]) %>%
    filter(knots == min(knots)) %>%
    left_join(var_set, by = c("model_covs" = "var")) %>%
    arrange(model_covs) %>%
    rename(best_k=knots)

#####  Run cross validation for null model (i.e. only province-specific random effects)
#####
cov_df <- expand.grid(var = "", year = cv_years)
dup_rows <- which(cov_df$var %in% cv_train$model_covs &
                      cov_df$year %in% cv_train$year)
cov_df <- if(length(dup_rows)>0) cov_df[-c(dup_rows),] else cov_df
if(nrow(cov_df)>0){
    tic <- Sys.time()
    null_train <- foreach(i=1:nrow(cov_df), .combine = rbind) %dopar% {
        predict_single_season(x = cv_dat,
                              vars = NULL,
                              test_year = cov_df$year[i],
                              loso_cv = TRUE,
                              knots = 0,
                              sets = 4,
                              distribution = FALSE)
    }
    cv_train <- bind_rows(cv_train, null_train)
    toc <- Sys.time()
}

if(TIMES) toc-tic
null_models <- cv_train %>%
    filter(model_covs=="", sets==4)
null_mae <- mean(abs(log(null_models$obs_counts) - log(null_models$pred_counts)))

#####  Run forward-backward algorithm to find the final models
#####

## take only the single-covariate splines with the best number of knots chosen above (to reduce the size of the data being carried around)
gam_train <- cv_train %>%
    filter(num_cov!=1|
               num_cov==1 &
               knots==table_dat$best_k[match(model_covs, table_dat$model_covs)])

if(is.null(best_models)){
    ## find the model with best MAE thus far
    cv_mae <- gam_train %>%
        group_by(formula, model_covs, num_cov) %>%
        summarise(mae = mean(abs(log(obs_counts)-log(pred_counts)))) %>%
        ungroup() %>%
        arrange(mae, num_cov)

    best_mae <- cv_mae$mae[1]
    cur_mae <- 0
    tic <- Sys.time()
    ## keep running algorithm until the model does not improve across one iteration
    while(cur_mae < best_mae){
        if(cur_mae > 0){
            ## find model with best MAE thus far
            cv_mae <- gam_train %>%
                group_by(formula, model_covs, num_cov) %>%
                summarise(mae = mean(abs(log(obs_counts) -
                                             log(pred_counts)))) %>%
                ungroup() %>%
                arrange(mae, num_cov)
        }
        best_mae <- cv_mae$mae[1]
        best_covs <- strsplit(as.character(cv_mae$model_covs[1]), split = ",")[[1]]
        if(VERBOSE) print(paste(Sys.time(), "Vars:", paste(best_covs, collapse=","), "MAE:", round(best_mae,2)))
        ## find all new models to compete with best existing MAE
        new_var <- unique(table_dat$model_covs)
        new_covs <- c()
        ## not pretty but will put in new variable or remove existing variable into new models
        for(i in 1:length(new_var))
            new_covs[i] = ifelse(new_var[i] %in% best_covs,
                                 paste(sort(best_covs[-which(best_covs==new_var[i])]),
                                       collapse=","),
                                 paste(sort(c(best_covs, new_var[i])), collapse=","))
        ## create a data frame for new model runs, removing duplicates from past runs
        cov_df <- expand.grid(var=new_covs,
                              year=cv_years)
        dup_rows <- which(cov_df$var %in% cv_mae$model_covs)
        cov_df <- if(length(dup_rows)>0) cov_df[-c(dup_rows),] else cov_df
        ## run cross validation on all new models
        if(nrow(cov_df)>0){
            new_train <- foreach(i=1:nrow(cov_df), .combine = rbind) %dopar% {
                if(VERBOSE)
                    print(paste("simulation", i, "of", nrow(cov_df)))
                predict_single_season(x = cv_dat,
                                      vars = cov_df$var[i],
                                      test_year = cov_df$year[i],
                                      loso_cv = TRUE,
                                      table_dat = table_dat,
                                      distribution = FALSE)
            }
            ## find the MAE for each new model
            new_mae <- new_train %>%
                group_by(formula, model_covs, num_cov) %>%
                summarise(mae = mean(abs(log(obs_counts)-
                                             log(pred_counts)))) %>%
                ungroup() %>%
                arrange(mae, num_cov)
            cur_mae <- new_mae$mae[1]
            gam_train <- bind_rows(gam_train, new_train)
        } else cur_mae <- best_mae
    }
    toc <- Sys.time()
}

if(TIMES) toc-tic
if(nrow(cov_df)>0)
    cv_train <- bind_rows(cv_train, gam_train) %>% distinct()
if(nrow(cov_df)>0)
    saveRDS(cv_train, paste0("data/cv-train-", Sys.Date(), ".rds"))
## find the best model for each number of covariates
cv_mae <- cv_train %>%
    filter(num_cov!=1|
               num_cov==1 &
               knots==table_dat$best_k[match(model_covs,table_dat$model_covs)]) %>%
    group_by(year, formula, model_covs, num_cov) %>%
    summarise(year_mae = mean(abs(log(obs_counts)-log(pred_counts)))) %>%
    group_by(formula, model_covs, num_cov) %>%
    summarise(mae = mean(year_mae),
              cv_sd = sd(year_mae)) %>%
    group_by(num_cov) %>%
    filter(mae == min(mae)) %>%
    ungroup() %>%
    arrange(mae, num_cov)

## the formula for the model that had the smallest CV MAE
best_cv_formula <- as.character(cv_mae$formula[1])

best_models <- cv_mae %>%
    select(formula, model_covs, num_cov, cv_mae=mae) %>%
    slice(1)

## find the best model within one standard deviation of the smallest CV MAE
reduced_model <- cv_mae %>%
    ungroup() %>%
    filter(num_cov<=num_cov[which.min(mae)],
           mae<=min(mae)+cv_sd[which.min(mae)]) %>%
    arrange(num_cov) %>%
    slice(1) %>%
    select(formula, model_covs, num_cov, cv_mae=mae)

best_models <- bind_rows(best_models, reduced_model)

best_model_covs <- strsplit(best_models$model_covs[1], ",")[[1]]

##### Find cross-validated MAE for model that predicts province median rate
#####
cv_median_mae <- cv_dat %>%
    group_by(pid) %>%
    mutate(median_rate = median(obs)) %>%
    ungroup() %>%
    summarise(formula = "",
              model_covs = "",
              num_cov = 0,
              cv_mae = mean(abs(log(obs) - log(median_rate))))

best_models <- bind_rows(best_models, cv_median_mae)
best_models$model <- c("Fitted",
                       "Simple",
                       "Baseline")

##### Run final and reduced models through the testing phase
#####

## create a data frame for all test years and models
test_years <- years[which(years >= FIRST_TEST_YEAR & years<=LAST_TEST_YEAR)]
cov_df <- expand.grid(var = best_models$model_covs, year = test_years) %>%
    filter(var != "")
dup_rows <- which(cov_df$var %in% test_preds$model_covs &
                      cov_df$year %in% test_preds$year)
cov_df <- if(length(dup_rows)>0) cov_df[-c(dup_rows),] else cov_df
## make forecasts for testing set
if(nrow(cov_df)>0){
    tic <- Sys.time()
    preds <- foreach(i=1:nrow(cov_df), .combine = rbind) %dopar% {
        if(VERBOSE) print(paste("simulation", i, "of", nrow(cov_df)))
        predict_single_season(x = dengue_dat,
                              vars = cov_df$var[i],
                              test_year = cov_df$year[i],
                              loso_cv = FALSE,
                              table_dat = table_dat,
                              distribution = TRUE,
                              coef_perms = 100)
    }
    test_preds <- bind_rows(test_preds, preds)
    toc <- Sys.time()
}

if(TIMES) toc-tic
if(nrow(cov_df)>0)
    saveRDS(test_preds, paste0("data/test-preds-", Sys.Date(), ".rds"))
