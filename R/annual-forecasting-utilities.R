#' Find recent file
#' Given a start of a name and a path to a folder, will return the name of the most recent file with that name in that folder
#'
#' @param name_start character string, first letters in file name
#' @param path character string, path to folder of interest, end with "/"
#'
#' @return character string, path to most recent file
#' @export
#'
#' @examples
find_recent_file <- function(name_start, path, exclude){
    if(substring(path, nchar(path))!="/")
        warning('Path does not end with a "/", problems may ensue.')
    ## view all files of that name at that path
    file_list <- list.files(path=path,
                            pattern=paste0(name_start, "*"))
    ## remove files with unwanted patterns
    if(!missing(exclude))
        file_list <- file_list[!grepl(pattern = exclude, file_list)]
    if(length(file_list)==0){
        warning('File not found')
        return(NA)
    }
    ## view file info
    file_info <- file.info(paste0(path, file_list))
    ## find most recent file
    most_recent_file <- paste0(path,
                               file_list[which.max(file_info$mtime)])
    return(most_recent_file)
}

#' Multivariate Normal random deviates
#' From mgcv examples in predict.gam()
#'
#' @param n number of draws from MVN
#' @param mu mean estimates of each normal in MVN
#' @param sig covariance matrix of MVN
rmvn <- function(n,mu,sig) { ## MVN random deviates
    L <- mroot(sig);m <- ncol(L);
    t(mu + L%*%matrix(rnorm(m*n),m,n))
}

#' Prediction intervals for negative binomial GAM regression
#' Randomly samples the multivariate normal distribution of the coefficient parameters and draws from a negative binomial distribution to create prediction intervals for a GAM spline of the negative binomial family.
#'
#'
#' @param train_fit gam object from mgcv's gam(), conducted with family = nb()
#' @param test_dat data to make predictions on
#' @param coef_perms number of normally distributed coefficient permutations
#' @param nb_draws number of draws from a negative binomial for each coefficient permutation
pred_dist <- function(train_fit,
                      test_dat,
                      coef_perms,
                      nb_draws){
    if(missing(nb_draws))
        nb_draws <- coef_perms
    ## get the mean and covariance matrix for predictions on the testing data
    test_preds <- predict(train_fit, test_dat, type = "lpmatrix")
    ## replicate parameter vectors
    rep_params <- rmvn(coef_perms,coef(train_fit),train_fit$Vp)
    pred_count <- matrix(0, nrow=nrow(test_dat), ncol=coef_perms)
    ## replicate predictions for each set of parameters
    for (i in 1:coef_perms) {
        rep_preds <- test_preds %*% rep_params[i,]
        pred_count[,i] <- exp(rep_preds)*test_dat$population
    }
    ## find the dispersion parameter for the negative binomial
    r <- train_fit$family$getTheta(trans=TRUE)
    ## sample from the negative binomial distribution
    all_preds <- matrix(rnbinom(nrow(pred_count)*coef_perms*nb_draws,
                                mu=as.vector(pred_count), size=r),
                        nrow=nrow(test_dat)) %>%
        as.data.frame()
    return(all_preds)
}

#' Predict single season
#' Predict DHF incidence for one year by training on  trained model to obtain confidence intervals and probability of outbreak
#'
#' @param x data.frame with training and test data
#' @param vars character vector of variables to be used
#' @param test_year integer, year to be predicted
#' @param loso_cv logical, is this leave-one-season-out cross validation (use all other data to predict one season)? Otherwise treated as a forecast (use only past data to predict one season)
#' @param knots numeric vector of the spline knots for each variable
#' @param sets numeric vector of the set for each variable (1=esrl, 2=ncdc, 3=noaa)
#' @param table_dat data.frame, specify a table with knots and sets for each variable instead of knots and sets separately
#' @param distribution logical, do you want a prediction distribution? Otherwise will just predict the expected value
#' @param coef_perms integer, number of normally distributed coefficient permutations for each GAM distribution
#' @param nb_draws integer, number of draws from a negative binomial for each coefficient permutation
#'
#' @return prediction_dat data.frame with the full formula used, variables, test year, and observed and predicted counts, and prediction interval (if desired)
#'
#' @examples
predict_single_season  <- function(x,
                                   vars,
                                   test_year,
                                   loso_cv=TRUE,
                                   knots,
                                   sets,
                                   table_dat,
                                   distribution=FALSE,
                                   coef_perms,
                                   nb_draws){
    require(dplyr)
    require(mgcv)
    require(tidyr)
    vars <- if(is.null(vars)) NULL else sort(strsplit(as.character(vars), split=",")[[1]])
    if(missing(knots) & missing(table_dat))
        stop("Please specify knots or a table with knots for each variable")
    if(!missing(table_dat))
        knots <- table_dat$best_k[which(table_dat$model_covs %in% vars)]
    if(missing(sets) & missing(table_dat))
        stop("Please specify sets or a table with set for each variable")
    if(!missing(table_dat))
        sets <- table_dat$set[which(table_dat$model_covs %in% vars)]
    if(distribution & missing(coef_perms))
        stop("Please specify the number of coef_perms and nb_draws for the distribution.")

    ## specify columns with variables of interest
    var_col <- which(colnames(x) %in% vars)
    ## make training data
    ## if leave-one-season-out cross validation, take out only test year
    if(loso_cv){
        train_dat <- cv_dat %>%
            filter(year != test_year) %>%
            select(1:4, var_col) %>%
            data.frame()
    } else { ## otherwise remove all data prior to test year
        train_dat <- x %>%
            filter(year < test_year) %>%
            select(1:4, var_col)
    }
    ## make test data
    test_dat <- x %>%
        filter(year == test_year) %>%
        select(1:4, var_col)

    ## set the null equation with an offset for population and province-level random effects
    null_eqn <- "obs ~ offset(log(population)) + s(pid, bs = 're')"
    ## list the variables for each source of the data (noaa, ncdc, or esrl)
    noaa_cov <- vars [which(sets<=3)]
    ncdc_cov <- vars[which(sets<=2)]
    esrl_cov <- vars[which(sets==1)]
    ## build a spline for each covariate
    test_splines <- paste0("s(", vars, ", k = ", knots, ", m=2, bs = 'cr')")

    if(length(noaa_cov)>0){
        noaa_eqn <- paste(null_eqn, paste(test_splines, collapse = "+"), sep="+")
    } else noaa_eqn <- null_eqn

    if(length(ncdc_cov)>0){
        ncdc_eqn <- paste(null_eqn, paste(test_splines[which(sets<=2)],
                                          collapse = "+"), sep="+")
    } else ncdc_eqn <- null_eqn

    if(length(esrl_cov)>0){
        esrl_eqn <- paste(null_eqn, paste(test_splines[which(sets<=1)],
                                          collapse = "+"), sep="+")
    } else esrl_eqn <- null_eqn

    ## remove NA values from noaa training and test sets
    noaa_train <- train_dat
    noaa_test <- test_dat
    if(length(noaa_cov)>0){
        for(i in 1:length(noaa_cov)){
            noaa_train_na_rows <- which(is.na(noaa_train[,noaa_cov[i]]))
            if(length(noaa_train_na_rows)>0)
                noaa_train <- noaa_train[-noaa_train_na_rows,]
            noaa_test_na_rows <- which(is.na(noaa_test[,noaa_cov[i]]))
            if(length(noaa_test_na_rows)>0)
                noaa_test <- noaa_test[-noaa_test_na_rows,]
        }}
    noaa_test <- filter(noaa_test, pid %in% noaa_train$pid)

    if(nrow(noaa_test)>0){
        ## fit model on noaa training data
        train_noaa_fit <- tryCatch(gam(formula(noaa_eqn),
                                       data = noaa_train,
                                       family = nb()),
                                   error=function(e) e)

        ## if the fit is an error, then fill the column with NAs, otherwise make predictions
        if(inherits(train_noaa_fit, "error")){
            noaa_preds <- matrix(NA, nrow = nrow(noaa_test), ncol = 1)
        } else{
            ## if using a distribution, draw from coefficient permutations and negative binomial
            if(distribution){
                noaa_preds <- pred_dist(train_fit=train_noaa_fit,
                                        test_dat=noaa_test,
                                        coef_perms=coef_perms,
                                        nb_draws=nb_draws)
            } else
                noaa_preds <- predict(train_noaa_fit, noaa_test, type = "response")
        }
        ## store data
        if(distribution){
            prediction_dat <- data_frame(year = noaa_test$year,
                                         pid = noaa_test$pid,
                                         formula = noaa_eqn,
                                         group = 3,
                                         model_covs = paste(vars, collapse=","),
                                         num_cov = length(vars),
                                         obs_counts = noaa_test$obs) %>%
                bind_cols(noaa_preds) %>%
                gather("sim", "pred_counts", starts_with("V")) %>%
                filter(!is.na(pred_counts))
        } else{
            prediction_dat <- data_frame(year = noaa_test$year,
                                         pid = noaa_test$pid,
                                         formula = noaa_eqn,
                                         group = 3,
                                         knots = sum(knots),
                                         model_covs = paste(vars, collapse=","),
                                         sets = sum(sets),
                                         num_cov = length(vars),
                                         obs_counts = noaa_test$obs,
                                         pred_counts = noaa_preds)
        }
    } else{
        prediction_dat <- c()
    }

    ## remove the noaa provinces from the ncdc test set (still used for training)
    ncdc_test <- test_dat %>% filter(!pid %in% prediction_dat$pid)
    ncdc_train <- train_dat
    ## remove NA values from ncdc training and test sets
    if(length(ncdc_cov)>0){
        for(i in 1:length(ncdc_cov)){
            ncdc_train_na_rows <- which(is.na(ncdc_train[,ncdc_cov[i]]))
            if(length(ncdc_train_na_rows)>0)
                ncdc_train <- ncdc_train[-ncdc_train_na_rows,]
            ncdc_test_na_rows <- which(is.na(ncdc_test[,ncdc_cov[i]]))
            if(length(ncdc_test_na_rows)>0)
                ncdc_test <- ncdc_test[-ncdc_test_na_rows,]
        }}
    ncdc_test <- filter(ncdc_test, pid %in% ncdc_train$pid)

    ## fit model on ncdc training data
    if(nrow(ncdc_test)>0){
        train_ncdc_fit <- tryCatch(gam(formula(ncdc_eqn), data = ncdc_train,
                                       family = nb()),
                                   error=function(e) e)
        ## if the fit is an error, then fill the column with NAs, otherwise make predictions
        if(inherits(train_ncdc_fit, "error")){
            ncdc_preds <- matrix(NA, nrow = nrow(ncdc_test), ncol = 1)
        } else{
            ## if using a distribution, draw from coefficient permutations and negative binomial
            if(distribution){
                ncdc_preds <- pred_dist(train_fit=train_ncdc_fit,
                                        test_dat=ncdc_test,
                                        coef_perms=coef_perms,
                                        nb_draws=nb_draws)
            } else
                ncdc_preds <- predict(train_ncdc_fit, ncdc_test, type = "response")
        }

        ## store data
        if(distribution){
            ncdc_dat <- data_frame(year = ncdc_test$year,
                                   pid = ncdc_test$pid,
                                   formula = noaa_eqn,
                                   group = 2,
                                   model_covs = paste(vars, collapse=","),
                                   num_cov = length(vars),
                                   obs_counts = ncdc_test$obs) %>%
                bind_cols(ncdc_preds) %>%
                gather("sim", "pred_counts", starts_with("V")) %>%
                filter(!is.na(pred_counts))
        } else{
            ncdc_dat <- data_frame(year = ncdc_test$year,
                                   pid = ncdc_test$pid,
                                   formula = noaa_eqn,
                                   group = 2,
                                   knots = sum(knots),
                                   model_covs = paste(vars, collapse=","),
                                   sets = sum(sets),
                                   num_cov = length(vars),
                                   obs_counts = ncdc_test$obs,
                                   pred_counts = ncdc_preds)
        }
        prediction_dat <- bind_rows(prediction_dat, ncdc_dat)
    }

    ## remove provinces with noaa or ncdc predictions from esrl testing set
    esrl_test <- test_dat %>% filter(!pid %in% prediction_dat$pid)
    esrl_train <- train_dat
    ## remove NA values from ncdc training and test sets
    if(length(esrl_cov)>0){
        for(i in 1:length(esrl_cov)){
            esrl_na_rows <- which(is.na(esrl_train[,esrl_cov[i]]))
            if(length(esrl_na_rows)>0)
                esrl_train <- esrl_train[-esrl_na_rows,]
        }}
    esrl_test <- filter(esrl_test, pid %in% esrl_train$pid)

    ## fit model on esrl training data
    if(nrow(esrl_test)>0){
        train_esrl_fit <- tryCatch(gam(formula(esrl_eqn), data = esrl_train,
                                       family = nb()),
                                   error=function(e) e)
        ## if the fit is an error, then fill the column with NAs, otherwise make predictions
        if(inherits(train_esrl_fit, "error")){
            esrl_preds <- matrix(NA, nrow = nrow(esrl_test), ncol = 1)
        } else{
            ## if using a distribution, draw from coefficient permutations and negative binomial
            if(distribution){
                esrl_preds <- pred_dist(train_fit=train_esrl_fit,
                                        test_dat=esrl_test,
                                        coef_perms=coef_perms,
                                        nb_draws=nb_draws)
            } else
                esrl_preds <- predict(train_esrl_fit, esrl_test, type = "response")
        }

        ## store data
        if(distribution){
            esrl_dat <- data_frame(year = esrl_test$year,
                                   pid = esrl_test$pid,
                                   formula = noaa_eqn,
                                   group = 1,
                                   model_covs = paste(vars, collapse=","),
                                   num_cov = length(vars),
                                   obs_counts = esrl_test$obs) %>%
                bind_cols(esrl_preds) %>%
                gather("sim", "pred_counts", starts_with("V"))
        } else{
            esrl_dat <- data_frame(year = esrl_test$year,
                                   pid = esrl_test$pid,
                                   formula = noaa_eqn,
                                   group = 1,
                                   knots = sum(knots),
                                   model_covs = paste(vars, collapse=","),
                                   sets = sum(sets),
                                   num_cov = length(vars),
                                   obs_counts = esrl_test$obs,
                                   pred_counts = esrl_preds)
        }
        prediction_dat <- bind_rows(prediction_dat, esrl_dat)
    }
    return(prediction_dat)
}
