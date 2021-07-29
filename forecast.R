## Dissolved oxygen forecast for submission to the EFI-NEON Forecast Challenge
## Forecast made using a state space model with atmospheric pressure as a covariate

## Author: Pricilla Ceja
## Minor edits by Alyssa Willson

rm(list = ls())

# Required packages 
library(tidyverse)
library(lubridate)
library(rjags)
library(tidybayes)
library(modelr)
library(aws.s3)
library(EFIstandards)
library(EML)
library(jsonlite)
library('devtools')
library('ncdf4')
library("reticulate")
#remotes::install_github("eco4cast/neon4cast")
library('neon4cast')
neonstore::neon_dir()

## Download meteorlogical driver data

## NEON data used for past time series
## Note: prior to conducting this step, the data must be downlaoded from the
## NEON website and saved to the working directory in a single folder

# Empty matrix to fill in loop
day_av = matrix(, nrow = 1404, ncol = 2)
# Counter for loop
count = 0

# Months and years that we will search for in the file names
files = c('2017-08', '2017-09', '2017-10', '2017-11', '2017-12',
          '2018-01', '2018-02', '2018-03', '2018-04', '2018-05', '2018-06',
          '2018-07', '2018-08', '2018-09', '2018-10', '2018-11', '2018-12',
          '2019-01', '2019-02', '2019-03', '2019-04', '2019-05', '2019-06',
          '2019-07', '2019-08', '2019-09', '2019-10', '2019-11', '2019-12',
          '2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06',
          '2020-07', '2020-08', '2020-09', '2020-10', '2020-11', '2020-12',
          '2021-01', '2021-02', '2021-03', '2021-04', '2021-05', '2021-06')

# Loop through the months and years that are all in different files
for(j in 1:length(files)){
  
  # List all folder names
  lf = list.files('NEON_pressure-air-buoy/.')
  # Find the folder name that matches the month
  num = grep(files[j], lf)
  # Save folder name
  foldnam = lf[num]
  # List files within that specific month's folder
  lf2 = list.files(paste0('NEON_pressure-air-buoy/',foldnam,'/.'))
  # Find the file name within that folder that matches 30 minute resolution
  num2 = grep('30min', lf2)
  # Save file name
  filnam = lf2[num2]
  # String together the folder and file name
  fn = paste0('NEON_pressure-air-buoy/',foldnam, '/', filnam)
  # Read in the CSV
  month = read.csv(fn)
  
  # Get days of the month in proper format
  unique_days = unique(substr(month$startDateTime, 1, 10))
  # Save number of days we are looping through
  ndays = length(unique_days)
  
  # Loop through each day in that specific CSV
  for(i in 1:length(unique_days)){
    # Increment counter for saving all averages in one matrix
    count = count + 1
    # Subset the CSV for one day of the month
    sub = subset(month, substr(month$startDateTime, 1, 10) == unique_days[i])
    # Find daily mean pressure
    day_av[count,] = cbind(unique_days[i], mean(as.numeric(sub$staPresMean), na.rm = T))
  }
}

## Formatting
# make data frame
day_av = as.data.frame(day_av)
# Change column names
colnames(day_av) = c('Date', 'Pressure')
# Make numeric
day_av$Pressure = as.numeric(day_av$Pressure)
# Change missing values from NaN to NA (easier to work with)
is.na(day_av$Pressure) = is.na(day_av$Pressure)


## NOAAGEFS Forecast of meteorological driver
siteID <- "BARC"
lat <- 29.675982
long <- -82.008414
dte <- as.Date("2021-06-30")
time_interval <- "1hr"
cycle <- "00"
download_noaa(siteID=siteID,interval=time_interval,date=dte,cycle=cycle,dir='../documents')

allpressure <- matrix(nrow=30,ncol=841,NA) #Each row is one 35-day ensemble member for different time values on the columns

for(i in 1:30){
  if(i<10){
    i <- paste("0",as.character(i),sep="")
  }
  if(as.numeric(i)>0){
    fileName <- paste("../documents/noaa/noaa/NOAAGEFS_",time_interval,"/",siteID,"/",dte,"/",cycle,
                      "/NOAAGEFS_",time_interval,"_",siteID,"_",dte,"T",cycle,"_",dte+35
                      ,"T",cycle,"_ens",i,".nc",sep="")
  }else{
    fileName <- paste("../documents/noaa/noaa/NOAAGEFS_",time_interval,"/",siteID,"/",dte,"/",cycle,
                      "/NOAAGEFS_",time_interval,"_",siteID,"_",dte,"T",cycle,"_",dte+16
                      ,"T",cycle,"_ens",i,".nc",sep="")
  }
  nc <- nc_open(fileName) #opening the file and putting it in nc
  time <- as.integer(ncdf4::ncvar_get(nc, "time")) #getting the column time from nc and putting it in a variable time
  tustr <- lubridate::as_datetime(strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value
                                           , " ")[[1]][3])
  time <- as.POSIXct.numeric((time*60*60), origin = tustr,tz = "UTC")
  pressure <- ncvar_get(nc,"air_pressure") #getting the air pressure section and putting it in pressure
  allpressure[as.numeric(i),] <- pressure #putting it in allpressure as numbers?
}


## setting up meteorological driver for use in model
dataATM <- (day_av$Pressure) * 0.00986923 # converting units in historical data
forecastATM <- (allpressure) * 9.86923e-6 # converting units in forecast
forecastMEAN <- apply(forecastATM, 1, mean)  # mean of forecast
# note that by taking the mean of the ensemble, 
# uncertainty in the meteorological driver is ignored

cov <- c(dataATM,forecastMEAN) # concatenating driver data and forecast into one time series


# set the random number for reproducible MCMC runs
set.seed(329)

# Generate plot to visualized forecast
generate_plots <- TRUE
# Is the forecast run on the Ecological Forecasting Initiative Server?
# Setting to TRUE published the forecast on the server
# efi_server <- TRUE

# List of team members. Used in the generation of the metadata
team_list <- list(list(individualName = list(givenName = "Pricilla", surName = "Ceja"),
                  list(individualName = list(givenName = "Alyssa",  surName ="Willson"))))


# Team name code
team_name <- "LSAMP REU AWPC (Team Name ID: LSAMP_AWPC)"

# Download target file from the server
download.file("https://data.ecoforecast.org/targets/aquatics/aquatics-targets.csv.gz",
              "aquatics-targets.csv.gz")

# Read in target file.  The guess_max is specified because there could be a lot of
# NA values at the beginning of the file
targets <- read_csv("aquatics-targets.csv.gz", guess_max = 10000)

# Focal sites
site_names <- c("BARC", "POSE")

# Forecast horizon 
f_days = 7 

## State space model modified from random walk model provided by the design team

state_space = "
model{
  # Priors
  x[1] ~ dnorm(x_ic,tau_init)
  tau_add ~ dgamma(0.1,0.1)
  tau_init ~ dgamma(0.1,0.1)
  beta_1 ~ dnorm(0, 0.001)
  beta_0 ~ dnorm(0, 0.001)
  beta_x ~ dnorm(0, 0.001)

  # Process Model
  for(t in 2:n){
 mu[t] <- x[t-1] + beta_0 + beta_1 * cov[t] + beta_x * x[t-1]
x[t] ~ dnorm(mu[t],tau_add) 
x_obs[t] ~ dnorm(x[t],tau_obs[t])

  }
  # Data Model
  for(i in 1:nobs){
    y[i] ~ dnorm(x[y_wgaps_index[i]], tau_obs[y_wgaps_index[i]])
  }
}
"

# Create variable for combined forecasts across sites
# Note: this is leftover from the provided code
# Only one site is considered here due to limited time duirng the REU program
forecast_saved_oxygen <- NULL
oxygen_figures <- list()

# Select BARC
site_data_var <- targets %>%
  filter(siteID == site_names[1])
  
# Last day in the observed data and add one day for the start of the forecast
start_forecast <- max(site_data_var$time) + days(1)
  
# The forecast horizon added to the end of the data
full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(f_days), by = "1 day"))
  
# Join the full time with the site_data_var so there aren't gaps in the time column
site_data_var <- left_join(full_time, site_data_var)
  
# observed oxygen: Full time series with gaps
y_wgaps <- site_data_var$oxygen
sd_wgaps <- imputeTS::na_interpolation(site_data_var$oxygen_sd,option = "linear")
time <- c(site_data_var$time)
# observed oxygen: time series without gaps
y_nogaps <- y_wgaps[!is.na(y_wgaps)]
# Index:number of the time series with gaps
y_wgaps_index <- 1:length(y_wgaps)
# Index: the index of the non-NA values in time series with gaps
y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
# Covarite Interpolation
cov_nogaps = imputeTS::na_interpolation(cov, option = 'linear')
# Note: ideally, this should be included in the model
# This step was skipped in the modeling process due to time constraints during the REU
  
# Generate starting initial conditions for latent states
init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
# Create a list of the data for use in JAGS.
data <- list(y = y_nogaps,
             y_wgaps_index = y_wgaps_index,
             nobs = length(y_wgaps_index),
             tau_obs = 1/(sd_wgaps ^ 2),
             n = length(y_wgaps),
             x_ic = 0.0,
             cov=cov_nogaps)
  
# Initialize parameters 
nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                    tau_init = mean( 1/var(diff(y_nogaps)), na.rm = TRUE),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i],
                    x = init_x)
}
  
# Initialize JAGS model
j.model   <- jags.model (file = textConnection(state_space),
                         data = data,
                         inits = init,
                         n.chains = 3)
  
# Run JAGS model as the burn-in
jags.out   <- coda.samples(model = j.model, 
                           variable.names = c("tau_add","tau_init"), 
                           n.iter = 10000)
  
# Run JAGS model again and sample from the posteriors
m   <- coda.samples(model = j.model,
                    variable.names = c("x","tau_add","tau_init", "x_obs"),
                    n.iter = 10000,
                    thin = 5)
  
# Use TidyBayes package to clean up the JAGS output. Better clean up
model_output <- m %>%
  spread_draws(x_obs[day]) %>%
  filter(.chain == 1) %>%
  rename(ensemble = .iteration) %>%
  mutate(time = full_time$time[day]) %>%
  ungroup() %>%
  select(time, x_obs, ensemble)
  
## Plotting
#Pull in the observed data for plotting
obs <- tibble(time = full_time$time, 
              obs = y_wgaps)
    
    
# Plots were divided into historical and future time series
# for a cleaner presentation of the model output

# Past plot
model_output %>% 
  group_by(time) %>% #group days together and able to plot
  summarise(mean = mean(x_obs), #good function, 
            upper = quantile(x_obs, 0.975),
            lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
  ggplot(aes(x = time, y = mean)) + #the plot
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
  geom_point(data = obs, aes(x = time, y = obs), color = "red") +
  labs(x = "Date", y = "oxygen (ATM)") +
  xlim(as.Date('2017-10-21'), as.Date('2020-12-31'))
   
# Present plot (includes only dates available in 2021)
model_output %>% 
  group_by(time) %>% #group days together and able to plot
  summarise(mean = mean(x_obs), #good function, 
            upper = quantile(x_obs, 0.975),
            lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
  ggplot(aes(x = time, y = mean)) + #the plot
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
  geom_point(data = obs, aes(x = time, y = obs), color = "red") +
  labs(x = "Date", y = "oxygen (ATM)") +
  xlim(as.Date('2021-01-01'), as.Date('2021-07-07'))
    
# Forecast plot
model_output %>% 
  group_by(time) %>% #group days together and able to plot
  summarise(mean = mean(x_obs), #good function, 
            upper = quantile(x_obs, 0.975),
            lower = quantile(x_obs, 0.025),.groups = "drop") %>% 
  ggplot(aes(x = time, y = mean)) + #the plot
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
  geom_point(data = obs, aes(x = time, y = obs), color = "red") +
  labs(x = "Date", y = "oxygen (ATM)") +
  xlim(as.Date('2021-07-01'), as.Date('2021-07-07'))

  
# Filter only the forecasted dates and add columns for required variable
# the forecast based on the missing data for each ensemble member
forecast_saved_tmp <- model_output %>%  #saving the model out put into new frame
  filter(time > start_forecast) %>% #going to reduce it and filter
  rename(oxygen = x_obs) %>% #literally renaming
  mutate(data_assimilation = 0,  #
         forecast = 1,
         obs_flag = 2,
         siteID = site_names[1]) %>%
  mutate(forecast_iteration_id = start_forecast) %>%
  mutate(forecast_project_id = team_name)
  
# Combined with the previous sites
forecast_saved_oxygen <- rbind(forecast_saved_oxygen, forecast_saved_tmp)

# Combined the oxygen and temperature forecasts together and re-order column
forecast_saved <- cbind(forecast_saved_oxygen)%>% 
  select(time, ensemble, siteID, oxygen, obs_flag, forecast, data_assimilation)

# Save file as CSV in the
# [theme_name]-[year]-[month]-[date]-[team_name].csv
# [aquatic_ecosystems]-[2021]-[July/August?]-[LSAMP REU AWPC (Team Name ID: LSAMP_AWPC)].csv
forecast_file_name_base <- paste0("aquatics-",as_date(start_forecast),"-",team_name)
forecast_file <- paste0(forecast_file_name_base, ".csv.gz")
write_csv(forecast_saved, forecast_file)

## Generate metadata

# Get system time for setting the issue time of the forecast
curr_time <- with_tz(Sys.time(), "UTC")
# Forecast_issue_time <- format(curr_time,format = "%Y-%m-%d %H:%M:%SZ", usetz = F)
forecast_issue_time <- as_date(curr_time)
forecast_iteration_id <- start_forecast

# The team name is the `forecast_model_id`
forecast_model_id <- team_name

source("metadata/generate_metadata.R")

meta_data_filename <- generate_metadata(forecast_file =  forecast_file,
                                        metadata_yaml = "metadata/metadata.yml",
                                        forecast_issue_time = as_date(with_tz(Sys.time(), "UTC")),
                                        forecast_iteration_id = start_forecast,
                                        forecast_file_name_base = forecast_file_name_base)


# Publish the forecast automatically.  Run only on EFI Challenge server
if(efi_server){
  source("../neon4cast-shared-utilities/publish.R")
  publish(code = "03_generate_null_forecast_aquatics.R",
          data_in = "aquatics-targets.csv.gz",
          data_out = forecast_file,
          meta = meta_data_filename,
          prefix = "aquatics/",
          bucket = "forecasts")
}