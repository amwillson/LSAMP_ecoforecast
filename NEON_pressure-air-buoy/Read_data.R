# Remove stuff from global environment
rm(list = ls())

# Set working directory
# CHANGE THIS LINE TO THE FOLDER YOU STORED ALL THOSE SUBDIRECTORIES FOR EACH MONTH
setwd('Downloads/NEON_pressure-air-buoy/')

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
  lf = list.files('.')
  # Find the folder name that matches the month
  num = grep(files[j], lf)
  # Save folder name
  foldnam = lf[num]
  # List files within that specific month's folder
  lf2 = list.files(paste0(foldnam,'/.'))
  # Find the file name within that folder that matches 30 minute resolution
  num2 = grep('30min', lf2)
  # Save file name
  filnam = lf2[num2]
  # String together the folder and file name
  fn = paste0(foldnam, '/', filnam)
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

# Formatting
# make data frame
day_av = as.data.frame(day_av)
# Change column names
colnames(day_av) = c('Date', 'Pressure')
# Make numeric
day_av$Pressure = as.numeric(day_av$Pressure)
# Change missing values from NaN to NA (easier to work with)
is.na(day_av$Pressure) = is.na(day_av$Pressure)
