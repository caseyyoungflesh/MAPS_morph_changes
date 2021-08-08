####################################
# Backtransform effect of time, lat, elev from IDX to raw mass and wing length
#
# Produce data for Table S2, S3
####################################

# top-level dir
dir <- 'XXXX'

si_tle_date <- '2021-07-28'
wi_tle_date <- '2021-07-28'
maps_master_date <- '2021-07-28'
daymet_process_date <- '2021-04-02'
run_date <- '2021-08-05'

#load packages
library(dplyr)
library(MCMCvis)

#read in data
setwd(paste0(dir, '/Data'))
mdata <- readRDS(paste0('MAPS-master-', maps_master_date, '.rds'))
usp <- sort(unique(mdata$sci_name))

#load tle model results and extract posteriors
setwd(paste0(dir, 'Results/si-tle-', si_tle_date))
si_fit <- readRDS(paste0('si-tle-fit-', si_tle_date, '.rds'))
si_data <- readRDS(paste0('si-tle-data-', si_tle_date, '.rds'))

si_eta_ch <- MCMCvis::MCMCchains(si_fit, params = 'eta')
si_gamma_ch <- MCMCvis::MCMCchains(si_fit, params = 'gamma')
si_theta_ch <- MCMCvis::MCMCchains(si_fit, params = 'theta')

rm(si_fit)
gc()

setwd(paste0(dir, 'Results/wi-tle-', wi_tle_date))
wi_fit <- readRDS(paste0('wi-tle-fit-', wi_tle_date, '.rds'))
wi_data <- readRDS(paste0('wi-tle-data-', wi_tle_date, '.rds'))

wi_eta_ch <- MCMCvis::MCMCchains(wi_fit, params = 'eta')
wi_gamma_ch <- MCMCvis::MCMCchains(wi_fit, params = 'gamma')
wi_theta_ch <- MCMCvis::MCMCchains(wi_fit, params = 'theta')

rm(wi_fit)
gc()


# backtransform -----------------------------------------------------------

#FROM 1-process-data.R
sl <- 0.3333609

#rotate points by slope of log(wl) log(mass) relationship
#get angle to rotate points (in radians)
angle <- -atan(sl)

#create rotation matrix
#rotm = [[cos(theta) -sin(theta]
#        [sin(theta)  cos(theta)]]
#90 degrees
# rotm <- matrix(c(0, -1, 
#                  1, 0), ncol = 2)

#forward rotation matrix
rotm <- matrix(c(cos(angle), sin(angle), 
                 -sin(angle), cos(angle)), ncol = 2)

#backward rotation matrix
brotm <- t(rotm)


# backtransform indices to traits -----------------------------------------

#function to log, rotate, then scale data (to check backtransform)
rotate_fun <- function(input)
{
  dm <- cbind(log(input$weight), log(input$wing_chord))
  M <- t(rotm %*% t(dm))
  M_sc <- apply(M, 2, function(x) scale(x, scale = TRUE))
  
  return(M_sc)
}

#function to derotate data
derotate_fun <- function(si_input, wi_input)
{
  M2 <- cbind(si_input, wi_input)
  s_dm2 <- t(brotm %*% t(M2))
  return(s_dm2)
}

#scale factors
#si_data$scf_YR #scale factor for year used in model
#si_data$scf_LAT #scale factor for lat used in model
#si_data$scf_ELEV #scale factor for elev used in model
if (si_data$scf_YR == wi_data$scf_YR &
    si_data$scf_LAT == wi_data$scf_LAT &
    si_data$scf_ELEV == wi_data$scf_ELEV)
{
  scf_YR <- si_data$scf_YR
  scf_LAT <- si_data$scf_LAT
  scf_ELEV <- si_data$scf_ELEV
} else {
  print("WARNING: scf don't match for SI and WI")
}

# l_res_param = change in log response for one unit change in covariate
# res_rng_param = ((e^param * L) - 1) * 100 = percent change in trait for L unit change in covariate, where L is range of cov covered for that species in dataset

# idx_param = change in index for each 1 unit change in cov

# eta = year
# gamma = lat
# theta = elev

rdf <- data.frame(sci_name = rep(NA, length(usp)),
                  sp_id = NA,
                  si_eta = NA,
                  si_gamma = NA,
                  si_theta = NA,
                  wi_eta = NA,
                  wi_gamma = NA,
                  wi_theta = NA,
                  l_mass_year = NA,
                  l_mass_lat = NA,
                  l_mass_elev = NA,
                  l_wing_year = NA,
                  l_wing_lat = NA,
                  l_wing_elev = NA,
                  mass_rng_year = NA,
                  mass_rng_lat = NA,
                  mass_rng_elev = NA,
                  wing_rng_year = NA,
                  wing_rng_lat = NA,
                  wing_rng_elev = NA)

#array for posterior realizations for index change in percent change in metric over time
rs_year_ch <- array(NA, dim = c(NROW(si_eta_ch), 2, length(usp)))
rs_lat_ch <- array(NA, dim = c(NROW(si_gamma_ch), 2, length(usp)))
rs_elev_ch <- array(NA, dim = c(NROW(si_theta_ch), 2, length(usp)))

for (i in 1:length(usp))
{
  #i <- 1
  print(paste0('processing species: ', i, ' of ', length(usp)))
  #filter by species
  temp <- dplyr::filter(mdata, sp_id == i)
  
  ###
  # #TO CHECK BACKTRANSFORM
  # #log, rotate, scale
  # rpd <- rotate_fun(temp)
  # #unscale
  # x_prime <- rpd[,1] * temp$sd_x[1] + temp$mn_x[1]
  # y_prime <- rpd[,2] * temp$sd_y[1] + temp$mn_y[1]
  # urpd <- cbind(x_prime, y_prime)
  # #derotate
  # lmet <- t(brotm %*% t(urpd))
  # #unlog (exponentiate)
  # met <- exp(lmet)
  # #values near zero confirm method is working as expected
  # hist(met[,1] - temp$weight)
  # hist(met[,2] - temp$wing_chord)
  ###
  
  #for each iteration in chain
  #divide by covariate scale (to get to 1 year, 1 deg lat, 1 m) and multiple by sd of idx
  #derotate data
  #effect of cov on log(response) - X unit change in log response for one unit change in covariate
  us_si_eta <- (si_eta_ch[,i] / scf_YR) * temp$sd_x[1]
  us_wi_eta <- (wi_eta_ch[,i] / scf_YR) * temp$sd_y[1]
  us_si_gamma <- (si_gamma_ch[,i] / scf_LAT) * temp$sd_x[1]
  us_wi_gamma <- (wi_gamma_ch[,i] / scf_LAT) * temp$sd_y[1]
  us_si_theta <- (si_theta_ch[,i] / scf_ELEV) * temp$sd_x[1]
  us_wi_theta <- (wi_theta_ch[,i] / scf_ELEV) * temp$sd_y[1]
  
  tdr_year_ch <- derotate_fun(us_si_eta, us_wi_eta)
  tdr_lat_ch <- derotate_fun(us_si_gamma, us_wi_gamma)
  tdr_elev_ch <- derotate_fun(us_si_theta, us_wi_theta)
  
  #INTERPRETATION
  #((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
  #((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
  
  #range of covariate for each species (except year which covers entire 30 year period)
  L_YR <- diff(range(mdata$year)) + 1
  L_LAT <- diff(range(temp$lat))
  L_ELEV <- diff(range(temp$GMTED_elev))
  
  mwc_year <- (exp(tdr_year_ch * L_YR) - 1) * 100
  mwc_lat <- (exp(tdr_lat_ch * L_LAT) - 1) * 100
  mwc_elev <- (exp(tdr_elev_ch * L_ELEV) - 1) * 100
  
  #fill arrays
  rs_year_ch[,,i] <- mwc_year
  rs_lat_ch[,,i] <- mwc_lat
  rs_elev_ch[,,i] <- mwc_elev
  
  #fill DF
  rdf$sp_id[i] <- i
  rdf$sci_name[i] <- usp[i]
  
  #effect of covariate on log response
  rdf$l_mass_year[i] <- mean(tdr_year_ch[,1])
  rdf$l_wing_year[i] <- mean(tdr_year_ch[,2])
  rdf$l_mass_lat[i] <- mean(tdr_lat_ch[,1])
  rdf$l_wing_lat[i] <- mean(tdr_lat_ch[,2])
  rdf$l_mass_elev[i] <- mean(tdr_elev_ch[,1])
  rdf$l_wing_elev[i] <- mean(tdr_elev_ch[,2])
  
  #percent change in response over range of cov for that species
  rdf$mass_rng_year[i] <- mean(mwc_year[,1])
  rdf$wing_rng_year[i] <- mean(mwc_year[,2])
  rdf$mass_rng_lat[i] <- mean(mwc_lat[,1])
  rdf$wing_rng_lat[i] <- mean(mwc_lat[,2])
  rdf$mass_rng_elev[i] <- mean(mwc_elev[,1])
  rdf$wing_rng_elev[i] <- mean(mwc_elev[,2])
  
  #raw param estimate
  rdf$si_eta[i] <- (mean(si_eta_ch[,i]) / scf_YR)
  rdf$si_gamma[i] <- (mean(si_gamma_ch[,i]) / scf_LAT)
  rdf$si_theta[i] <- (mean(si_theta_ch[,i]) / scf_ELEV)
  rdf$wi_eta[i] <- (mean(wi_eta_ch[,i]) / scf_YR)
  rdf$wi_gamma[i] <- (mean(wi_gamma_ch[,i]) / scf_LAT)
  rdf$wi_theta[i] <- (mean(wi_theta_ch[,i]) / scf_ELEV)
}


# save results to RDS -----------------------------------------------------

saveRDS(rdf, paste0(dir, 'Results/per-ch-table-data-', run_date, '.rds'))

saveRDS(rs_year_ch, paste0(dir, 'Results/per-ch-year-array-', run_date, '.rds'))
saveRDS(rs_lat_ch, paste0(dir, 'Results/per-ch-lat-array-', run_date, '.rds'))
saveRDS(rs_elev_ch, paste0(dir, 'Results/per-ch-elev-array-', run_date, '.rds'))
