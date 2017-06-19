## Prepare a data frame to make annual predictions from
## Stephen Lauer, March 2017

library(dengueThailand)
library(dplyr)
# devtools::install_github("krlmlr/here")

data("thai_prov_data")
## adjust thai_prov_data for Bueng Kan as part of Nong Khai
thai_prov_data[thai_prov_data$ISO==43,"Area_km"] <-
    thai_prov_data[thai_prov_data$ISO==43,"Area_km"] +
    thai_prov_data[thai_prov_data$ISO==38,"Area_km"]
data("thai_census_interpolated")

FIRST_YEAR <- 2000
LAST_YEAR <- 2014
FIRST_TEST_YEAR <- 2010
LAST_LOW_MONTH <- 3
FIRST_LOW_MONTH <- 11

setwd(here::here())
source("R/annual-forecasting-utilities.R")

## read in the monthly counts
recent_monthly_counts <- find_recent_file(name_start = "monthly-counts",
                                          path = "data/")
## adjust counts for Bueng Kan as part of Nong Khai
formatted_counts <- read.csv(recent_monthly_counts) %>%
    mutate(ISO=ifelse(ISO==38,43,ISO)) %>%
    left_join(select(thai_prov_data, ISO, FIPS)) %>%
    select(pid=FIPS, date_sick_year, date_sick_month, count)

## read in Thai censes
annual_population <- thai_census_interpolated %>%
    filter(year>=FIRST_YEAR-1,
           year<=LAST_YEAR) %>%
    left_join(select(thai_prov_data, FIPS, Area_km),
              by=c("pid"="FIPS")) %>%
    mutate(population_density = population/Area_km) %>%
    select(-Area_km)

## fill in all zero count months
all_counts <- expand.grid(pid = unique(formatted_counts$pid),
                          year = unique(formatted_counts$date_sick_year),
                          month = unique(formatted_counts$date_sick_month)) %>%
    left_join(formatted_counts,
              by = c("pid", "year"="date_sick_year", "month"="date_sick_month")) %>%
    mutate(count = ifelse(is.na(count), 0, count))


## aggregate counts from the preseason, high season, and postseason for each year, and make all counts into per 100,000 population rates
high_season_counts <- all_counts %>%
    group_by(pid, year) %>%
    summarise(pre_season_counts =
                  sum(count[which(month<=LAST_LOW_MONTH)]),
              high_season_counts =
                  sum(count[which(month>LAST_LOW_MONTH &
                                      month<FIRST_LOW_MONTH)]),
              post_season_counts =
                  sum(count[which(month>=FIRST_LOW_MONTH)])) %>%
    ungroup() %>%
    full_join(annual_population, by = c("pid", "year")) %>%
    group_by(pid) %>%
    mutate(next_pre_rate = lead(pre_season_counts),
           next_pre_rate = 100000*next_pre_rate/population,
           post_season_rate = 100000*post_season_counts/population,
           high_season_rate = 100000*high_season_counts/population,
           pre_season_rate = lag(next_pre_rate),
           last_post_rate = lag(post_season_rate),
           last_high_rate = lag(high_season_rate)) %>%
    select(pid, year, population, population_density, pre_season_rate,
           last_post_rate, last_high_rate)

## find the estimated relative susceptibility for each province
susceptibles <- high_season_counts %>%
    mutate(annual_rate = lead(pre_season_rate + last_post_rate +
                                   last_high_rate)) %>%
    filter(!is.na(annual_rate)) %>%
    select(pid, year, annual_rate)

all_sus <- c()
for(prov in unique(susceptibles$pid)){
    tmp_sus <- filter(susceptibles, pid == prov)
    N_bar <- mean(tmp_sus$annual_rate[which(tmp_sus$year<(FIRST_TEST_YEAR-1))])
    tmp_sus$S <- NA
    tmp_sus$S[1] <- N_bar
    tmp_sus$S[2] <- tmp_sus$S[1] - tmp_sus$annual_rate[1] + N_bar
    tmp_sus$S[3] <- tmp_sus$S[2] - tmp_sus$annual_rate[2] + N_bar
    for(j in 4:nrow(tmp_sus))
        tmp_sus$S[j] <- tmp_sus$S[j-1] - tmp_sus$annual_rate[j-1] + tmp_sus$annual_rate[j-3]
    all_sus <- c(all_sus, tmp_sus$S)
}

susceptibles$sus_pct <- all_sus

susceptibles <- select(susceptibles, -annual_rate)

## put together annual count, high, preseason, and postseason counts, and susceptibles into one dataframe
annual_counts <- all_counts %>%
    group_by(pid, year) %>%
    summarise(annual_count = sum(count[which(month>LAST_LOW_MONTH)])) %>%
    select(pid, year, annual_count) %>%
    left_join(high_season_counts, by = c("pid", "year")) %>%
    left_join(susceptibles, by = c("pid", "year")) %>%
    filter(year >= FIRST_YEAR,
           year <= LAST_YEAR) %>%
    ungroup() %>%
    rename(obs = annual_count)

## read in most recently aggregated weather data
most_recent_weather <- find_recent_file(name_start = "all-weather-dat",
                                        path = "data/")
weather_dat <- read.csv(most_recent_weather)[,-1] %>%
    filter(!is.na(pid))

## join count data with weather data
dengue_dat <- left_join(annual_counts, weather_dat, by = c("pid", "year")) %>%
    ungroup() %>%
    mutate(pid = as.factor(pid)) %>%
    filter(!is.na(obs)) %>%
    mutate(obs = ifelse(obs==0, 1, obs))

write.csv(dengue_dat, file=paste0("data/annual-pred-data-", Sys.Date(),".csv"))
