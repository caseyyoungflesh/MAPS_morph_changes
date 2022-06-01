#######################
# process daymet data
#
# https://daymet.ornl.gov/
# 
#######################


# input -------------------------------------------------------------------

dir <- 'XXXX'
MAPS_process_date <- '2021-05-31'
daymet_process_date <- '2021-04-02'


# load packages -----------------------------------------------------------

library(dplyr)


# read in MAPS data -------------------------------------------------------

setwd(paste0(dir, '/Data'))
maps_data <- readRDS(paste0('MAPS-master-', MAPS_process_date, '.rds'))


# download point values for daymet from ORNL ------------------------------

# #run bash script to DL daymet data from individual points: https://github.com/ornldaac/daymet-single-pixel-batch/tree/master/bash
ifelse(!dir.exists(paste0(dir, 'Data/daymet/RAW/', daymet_process_date)),
       dir.create(paste0(dir, 'Data/daymet/RAW/', daymet_process_date)),
       FALSE)

setwd(paste0(dir, 'Data/daymet/RAW/', daymet_process_date))
system(paste0(dir, 'Scripts/3-process-daymet/daymet_spt_batch.sh ',
              dir, 'Scripts/3-process-daymet/daymet_query.txt'))


# combine csv files -------------------------------------------------------

setwd(paste0(dir, 'Data/daymet/RAW/', daymet_process_date))
files <- list.files()
N <- NROW(read.csv(files[1], skip = 7)) * length(files)
fout <- data.frame(year = rep(NA, N), 
                   yday = NA,
                   tmax = NA,
                   lat = NA, 
                   lng = NA)
counter <- 1
for (i in 1:length(files))
{
  #i <- 1
  print(paste0('file ', i, ' of ', length(files)))
  
  tt <- read.csv(files[i], skip = 7)
  t_idx <- counter:(counter + NROW(tt) - 1)
  
  fout$year[t_idx] <- tt$year
  fout$yday[t_idx] <- tt$yday
  fout$tmax[t_idx] <- tt$tmax..deg.c.
  
  splt_nm <- strsplit(files[i], '_')
  fout$lat[t_idx] <- as.numeric(splt_nm[[1]][3])
  fout$lng[t_idx] <- as.numeric(splt_nm[[1]][5])
  
  counter <- counter + NROW(tt)
}


# save daymet RDS ----------------------------------------------------------------

setwd(paste0(dir, '/Data/daymet/processed'))
saveRDS(fout, paste0('daymet-processed-', daymet_process_date, '.rds'))


# Mean abiotic data --------------------------------------------------------

#average prcp, tmin, and tmax for each station
ust <- unique(maps_data[,c('station', 'lng', 'lat')])
env_out <- data.frame(station = rep(NA, (NROW(ust) * length(unique(fout$year)))),
                      year = NA,
                      May_tmax = NA,
                      June_tmax = NA,
                      July_tmax = NA)
counter <- 1
for (i in 1:NROW(ust))
{
  #i <- 1
  print(paste0(i, ' of ', NROW(ust)))
  
  #filter env by lat/lng of station
  tenv <- dplyr::filter(fout, lat == ust$lat[i], lng == ust$lng[i])
  
  #May 1 - May 30
  drng_May <- 121:151
  #June 1 - June 30
  drng_June <- 152:181
  #July 1 - July 31
  drng_July <- 182:212
 
  yrs <- unique(tenv$year)
  for (j in 1:length(yrs))
  {
    #j <- 1
    tval_May <- dplyr::filter(tenv, yday %in% drng_May, year == yrs[j])
    tval_June <- dplyr::filter(tenv, yday %in% drng_June, year == yrs[j])
    tval_July <- dplyr::filter(tenv, yday %in% drng_July, year == yrs[j])
    
    #fill df
    env_out$station[counter] <- ust$station[i]
    env_out$year[counter] <- yrs[j]
    env_out$May_tmax[counter] <- mean(tval_May$tmax)
    env_out$June_tmax[counter] <- mean(tval_June$tmax)
    env_out$July_tmax[counter] <- mean(tval_July$tmax)

    counter <- counter + 1
  }
}


# save RDS object -----------------------------------------------------

setwd(paste0(dir, '/Data/daymet/processed'))
saveRDS(env_out, paste0('daymet-master-', daymet_process_date, '.rds'))

