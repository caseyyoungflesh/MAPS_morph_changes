####################
# Analyze data and create figs
#
# 
####################


# top-level dir --------------------------------------------------------------

dir <- 'XXXX'

#SI over time, lat, elev
si_tle_date <- '2021-07-28'
#WI over time, lat, elev
wi_tle_date <- '2021-07-28'
#SI in response to spatial variation in temp
si_temp_space_date <- '2022-05-25'
#SI in response to temporal variation in temp
si_temp_date <- '2022-05-25'
#Spatial responses as a function of cov
si_temp_sp_cov_date <- '2022-06-01'
#Date Daymet data were processed
daymet_process_date <- '2021-04-02'
#Date MAPS data were processed
maps_master_date <- '2021-07-28'
#date per ch backtransform was run
per_ch_date <- '2021-08-05'

fig_dir <- paste0(dir, '/Figures/paper_figs/')


# load packages -----------------------------------------------------------

library(dplyr)
library(MCMCvis)
library(ggplot2)
library(viridis)
library(tidyr)
library(raster)
library(tmap)
library(sf)


# load data ---------------------------------------------------------------

#read in data
mdata <- readRDS(paste0(dir, 'Data/MAPS-master-', maps_master_date, '.rds'))
usp <- sort(unique(mdata$sci_name))


# correlation between log metrics and indices -----------------------------

si_cor <- rep(NA, length(usp))
wi_cor <- rep(NA, length(usp))
for (i in 1:length(usp)) 
{
  #i <- 1
  tt <- dplyr::filter(mdata, sci_name == usp[i])
  si_cor[i] <- cor(log(tt$weight), tt$size_idx)
  wi_cor[i] <- cor(log(tt$wing_chord), tt$wingy_idx)
}

mean(si_cor)
range(si_cor)
mean(wi_cor)
range(wi_cor)


# load model results ------------------------------------------------------

#SI ~ TLE
setwd(paste0(dir, 'Results/si-tle-', si_tle_date))
si_fit <- readRDS(paste0('si-tle-fit-', si_tle_date, '.rds'))
si_data <- readRDS(paste0('si-tle-data-', si_tle_date, '.rds'))

#WI ~ TLE
setwd(paste0(dir, 'Results/wi-tle-', wi_tle_date))
wi_fit <- readRDS(paste0('wi-tle-fit-', wi_tle_date, '.rds'))
wi_data <- readRDS(paste0('wi-tle-data-', wi_tle_date, '.rds'))

#SI ~ temp space
setwd(paste0(dir, 'Results/si-temp-space-', si_temp_space_date))
si_temp_space_fit <- readRDS(paste0('si-temp-space-fit-', si_temp_space_date, '.rds'))
si_temp_space_data <- readRDS(paste0('si-temp-space-data-', si_temp_space_date, '.rds'))

#Temp response space ~ temp experienced
setwd(paste0(dir, 'Results/temp-sp-cov-', si_temp_sp_cov_date))
si_temp_sp_cov_fit <- readRDS(paste0('temp-sp-cov-fit-', si_temp_sp_cov_date, '.rds'))
si_temp_sp_cov_data <- readRDS(paste0('temp-sp-cov-data-', si_temp_sp_cov_date, '.rds'))

#SI ~ temp lag 0
setwd(paste0(dir, 'Results/si-temp-l0-', si_temp_date))
si_temp_l0_fit <- readRDS(paste0('si-temp-l0-fit-', si_temp_date, '.rds'))
si_temp_l0_data <- readRDS(paste0('si-temp-l0-data-', si_temp_date, '.rds'))

#SI ~ temp lag 1
setwd(paste0(dir, 'Results/si-temp-l1-', si_temp_date))
si_temp_l1_fit <- readRDS(paste0('si-temp-l1-fit-', si_temp_date, '.rds'))
si_temp_l1_data <- readRDS(paste0('si-temp-l1-data-', si_temp_date, '.rds'))

#SI ~ temp lag 2
setwd(paste0(dir, 'Results/si-temp-l2-', si_temp_date))
si_temp_l2_fit <- readRDS(paste0('si-temp-l2-fit-', si_temp_date, '.rds'))
si_temp_l2_data <- readRDS(paste0('si-temp-l2-data-', si_temp_date, '.rds'))

#change over X number of units
SYR <- 10
SLAT <- 10
SELEV <- 1000

#read in backtransformation
per_ch_year <- readRDS(paste0(dir, 'Results/per-ch-year-array-', per_ch_date, '.rds'))
per_ch_lat <- readRDS(paste0(dir, 'Results/per-ch-lat-array-', per_ch_date, '.rds'))
per_ch_elev <- readRDS(paste0(dir, 'Results/per-ch-elev-array-', per_ch_date, '.rds'))


# SI TIME param est -------------------------------------------------------

#SI TIME
me_ch_SI <- (MCMCvis::MCMCchains(si_fit, params = 'mu_eta') / si_data$scf_YR) * SYR
round(mean(me_ch_SI), 2)
round(apply(me_ch_SI, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(me_ch_SI < 0) / NROW(me_ch_SI), 2)


# SI TEMP lag param est -------------------------------------------------------

#SI TEMP LAG 0 change
mg_ch_SI_L0 <- (MCMCvis::MCMCchains(si_temp_l0_fit, params = 'mu_gamma')) / si_temp_l0_data$scf_temp
round(mean(mg_ch_SI_L0), 2)
round(apply(mg_ch_SI_L0, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mg_ch_SI_L0 < 0) / NROW(mg_ch_SI_L0), 2)

#SI TEMP LAG 0 effect of temp as function of mn temp
mt_ch_SI_L0 <- (MCMCvis::MCMCchains(si_temp_l0_fit, params = 'mu_theta')) / si_temp_l0_data$scf_temp
round(mean(mt_ch_SI_L0), 2)
round(apply(mt_ch_SI_L0, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mt_ch_SI_L0 < 0) / NROW(mt_ch_SI_L0), 2)

#SI TEMP LAG 1
mg_ch_SI_L1 <- (MCMCvis::MCMCchains(si_temp_l1_fit, params = 'mu_gamma')) / si_temp_l0_data$scf_temp
round(mean(mg_ch_SI_L1), 2)
round(apply(mg_ch_SI_L1, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mg_ch_SI_L1 < 0) / NROW(mg_ch_SI_L1), 2)

#SI TEMP LAG 1 change as function of mn temp
mt_ch_SI_L1 <- (MCMCvis::MCMCchains(si_temp_l1_fit, params = 'mu_theta')) / si_temp_l0_data$scf_temp
round(mean(mt_ch_SI_L1), 2)
round(apply(mt_ch_SI_L1, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mt_ch_SI_L1 < 0) / NROW(mt_ch_SI_L1), 2)


#SI TEMP LAG 2
mg_ch_SI_L2 <- (MCMCvis::MCMCchains(si_temp_l2_fit, params = 'mu_gamma')) / si_temp_l0_data$scf_temp
round(mean(mg_ch_SI_L2), 2)
round(apply(mg_ch_SI_L2, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mg_ch_SI_L2 < 0) / NROW(mg_ch_SI_L2), 2)

#SI TEMP LAG 2 change as function of mn temp
mt_ch_SI_L2 <- (MCMCvis::MCMCchains(si_temp_l2_fit, params = 'mu_theta')) / si_temp_l0_data$scf_temp
round(mean(mt_ch_SI_L2), 2)
round(apply(mt_ch_SI_L2, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mt_ch_SI_L2 < 0) / NROW(mt_ch_SI_L2), 2)

#species-specific - per 1 degree C
g_ch_SI_L0 <- (MCMCvis::MCMCchains(si_temp_l0_fit, params = 'gamma')) / si_temp_l0_data$scf_temp
g_ch_SI_L1 <- (MCMCvis::MCMCchains(si_temp_l1_fit, params = 'gamma')) / si_temp_l1_data$scf_temp
sum(apply(g_ch_SI_L0, 2, mean) < 0) / NCOL(g_ch_SI_L0)


# WI TIME param est -------------------------------------------------------

#WI TIME
me_ch_WI <- (MCMCvis::MCMCchains(wi_fit, params = 'mu_eta') / wi_data$scf_YR) * SYR
round(mean(me_ch_WI), 2)
round(apply(me_ch_WI, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(me_ch_WI > 0) / NROW(me_ch_WI), 2)


# SI LAT param est --------------------------------------------------------

#SI LAT
mg_ch_SI <- (MCMCvis::MCMCchains(si_fit, params = 'mu_gamma') / si_data$scf_LAT) * SLAT
round(mean(mg_ch_SI), 2)
round(apply(mg_ch_SI, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mg_ch_SI > 0) / NROW(mg_ch_SI), 2)


# SI TEMP space param est -------------------------------------------------

#SI TEMP space change
mb_ch_SI_space <- (MCMCvis::MCMCchains(si_temp_space_fit, params = 'mu_beta'))
round(mean(mb_ch_SI_space), 2)
round(apply(mb_ch_SI_space, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mb_ch_SI_space < 0) / NROW(mb_ch_SI_space), 2)

#SI TEMP space cov
g_ch_SI_sp_cov <- (MCMCvis::MCMCchains(si_temp_sp_cov_fit, params = 'beta'))
round(mean(g_ch_SI_sp_cov), 2)
round(apply(g_ch_SI_sp_cov, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(g_ch_SI_sp_cov < 0) / NROW(g_ch_SI_sp_cov), 2)

#species-specific
b_ch_SI_space <- (MCMCvis::MCMCchains(si_temp_space_fit, params = 'beta'))

            
# compare response to space temp to temporal temp -------------------------

#space is 10 degrees, time is 1 degree
#subtract means and add variances
delta_l0_mn <- mean(mb_ch_SI_space) - mean(mg_ch_SI_L0 * 10)
delta_l0_sd <- as.numeric(sqrt(var(mb_ch_SI_space) + var(mg_ch_SI_L0 * 10)))

#species-specific
space_samp <- sample(1:NROW(b_ch_SI_space), 5000)
l0_samp <- sample(1:NROW(g_ch_SI_L0), 5000)

delta_ss_mn <- apply(b_ch_SI_space, 2, mean) - apply(g_ch_SI_L0 * 10, 2, mean)
delta_ss_sd <- sqrt(apply(b_ch_SI_space, 2, var) + apply(g_ch_SI_L0 * 10, 2, var))

#neg = larger negative impact of temp over space
sum(delta_ss_mn < 0) / length(delta_ss_mn)

#plot spatial response vs. temporal response
b_mn <- apply(b_ch_SI_space, 2, mean)
b_sd <- apply(b_ch_SI_space, 2, sd)
g_mn <- apply(g_ch_SI_L0 * 10, 2, mean)
g_sd <- apply(g_ch_SI_L0 * 10, 2, sd)

st_data <- data.frame(b_mn, b_sd, g_mn, g_sd)
st2_data <- data.frame(mb_mn = mean(mb_ch_SI_space[,1]),
                       mb_sd = sd(mb_ch_SI_space[,1]),
                       mg_mn = mean(mg_ch_SI_L0[,1]*10),
                       mg_sd = sd(mg_ch_SI_L0[,1]*10))

#spatial vs. temporal
st_temp_plt <- ggplot(data = st_data, 
                      aes(b_mn, g_mn)) +
  geom_vline(xintercept = 0, linetype = 'dashed', 
             col = 'black', size = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed', 
             col = 'black', size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = g_mn - g_sd,
                    ymax = g_mn + g_sd),
                width = 0, 
                size = 0.7,
                color = 'black', alpha = 0.2) +
  geom_errorbarh(aes(xmin = b_mn - b_sd,
                     xmax = b_mn + b_sd),
                 height = 0, 
                 size = 0.7,
                 color = 'black', alpha = 0.2) +
  geom_point(color = 'black', size = 3, alpha = 0.3) +
  geom_point(data = st2_data, 
             aes(mb_mn, mg_mn), 
             color = 'red', alpha = 0.8, size = 30,
             pch = '*') +
theme_bw() +
  xlab('Temp sens over space (SI / 10 degree C)') +
  ylab('Temp sens over time (SI / 10 degree C)') +
  ylim(c(-0.4, 0.2)) +
  xlim(c(-2, 1.5)) +
  # geom_abline(slope = 1, #linetype = 'dotted',
  #             col = 'red', size = 1, alpha = 0.8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
        axis.ticks.length = unit(0.2, 'cm')) #length of axis tick

#spatial - temporal
tplt <- data.frame(mn = delta_ss_mn)
tplt2 <- data.frame(dmn = delta_l0_mn,
                    dsd = delta_l0_sd)
pnt_plt <- ggplot(tplt, aes(0, mn)) +
  geom_jitter(shape = 1, alpha = 0.5, size = 4, height = 0,
              width = 0.8) +
  geom_hline(yintercept = 0,
             linetype = 'dashed',
             size = 2,
             alpha = 0.4) +
  geom_point(aes(0, tplt2$dmn), inherit.aes = FALSE, color = 'red', pch = '*',
             size = 40, alpha = 0.05) +
  xlim(-2, 2) +
  ylab('Spatial temp response - temporal temp response') +
  xlab('') +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
        axis.ticks.length = unit(0.2, 'cm')) #length of axis tick
            

# WI LAT param est --------------------------------------------------------------------

#WI LAT
mg_ch_WI <- (MCMCvis::MCMCchains(wi_fit, params = 'mu_gamma') / wi_data$scf_LAT) * SLAT
round(mean(mg_ch_WI), 2)
round(apply(mg_ch_WI, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mg_ch_WI > 0) / NROW(mg_ch_WI), 2)


# SI ELEV param est -------------------------------------------------------

#SI ELEV
mt_ch_SI <- (MCMCvis::MCMCchains(si_fit, params = 'mu_theta') / si_data$scf_ELEV) * SELEV
round(mean(mt_ch_SI), 2)
round(apply(mt_ch_SI, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mt_ch_SI < 0) / NROW(mt_ch_SI), 2)


# WI ELEV param est -------------------------------------------------------

#WI ELEV
mt_ch_WI <- (MCMCvis::MCMCchains(wi_fit, params = 'mu_theta') / wi_data$scf_ELEV) * SELEV
round(mean(mt_ch_WI), 2)
round(apply(mt_ch_WI, 2, function(x) quantile(x, probs = c(0.055, 0.945))), 2)
round(sum(mt_ch_WI > 0) / NROW(mt_ch_WI), 2)


# Data metrics -----------------------------------------------------------------

#DATA STATS
# num samples - 253,488
NROW(mdata)
# num stations - 1124
length(unique(mdata$station))
# num species - 105
length(unique(mdata$sci_name))
# range years - 1989-2018
range(mdata$year)
# range lat - 26.1 - 69.4
range(mdata$lat)
# range elev - 0 - 2996
range(mdata$GMTED_elev)


# mass ~ time back transformation -----------------------------------------

#posterior for mean percent change over all species
mn_mass_per_time_ch <- apply(per_ch_year[,1,], 1, mean)
round(mean(mn_mass_per_time_ch), 2)
round(quantile(mn_mass_per_time_ch, probs = c(0.055, 0.945)), 2)
round(sum((mn_mass_per_time_ch < 0)) / length(mn_mass_per_time_ch), 2)

#max change species
mass_per_time_sp <- apply(per_ch_year[,1,], 2, mean)
mass_per_time_sp_qt <- apply(per_ch_year[,1,], 2, function(x) quantile(x, probs = c(0.055, 0.945)))
midx <- which.min(mass_per_time_sp)
usp[midx]
round(mass_per_time_sp[midx], 2)
round(mass_per_time_sp_qt[,midx], 2)
round(sum(per_ch_year[,1,midx] < 0) / NROW(per_ch_year[,1,]), 2)


# wing ~ time back transformation -----------------------------------------

mn_wing_per_time_ch <- apply(per_ch_year[,2,], 1, mean)
round(mean(mn_wing_per_time_ch), 2)
round(quantile(mn_wing_per_time_ch, probs = c(0.055, 0.945)), 2)
round(sum((mn_wing_per_time_ch < 0)) / length(mn_wing_per_time_ch), 2)


# mass ~ lat back transformation -----------------------------------------

mn_mass_per_lat_ch <- apply(per_ch_lat[,1,], 1, mean)
round(mean(mn_mass_per_lat_ch), 2)
round(quantile(mn_mass_per_lat_ch, probs = c(0.055, 0.945)), 2)
round(sum((mn_mass_per_lat_ch > 0)) / length(mn_mass_per_lat_ch), 2)


# wing ~ lat back transformation -----------------------------------------

mn_wing_per_lat_ch <- apply(per_ch_lat[,2,], 1, mean)
round(mean(mn_wing_per_lat_ch), 2)
round(quantile(mn_wing_per_lat_ch, probs = c(0.055, 0.945)), 2)
round(sum((mn_wing_per_lat_ch > 0)) / length(mn_wing_per_lat_ch), 2)


# mass ~ elev back transformation -----------------------------------------

mn_mass_per_elev_ch <- apply(per_ch_elev[,1,], 1, mean)
round(mean(mn_mass_per_elev_ch), 2)
round(quantile(mn_mass_per_elev_ch, probs = c(0.055, 0.945)), 2)
round(sum((mn_mass_per_elev_ch < 0)) / length(mn_mass_per_elev_ch), 2)


# wing ~ elev back transformation -----------------------------------------

mn_wing_per_elev_ch <- apply(per_ch_elev[,2,], 1, mean)
round(mean(mn_wing_per_elev_ch), 2)
round(quantile(mn_wing_per_elev_ch, probs = c(0.055, 0.945)), 2)
round(sum((mn_wing_per_elev_ch > 0)) / length(mn_wing_per_elev_ch), 2)


# simulate SI and WI over time,  lat,  elev -------------------------------

#si - posterior means
si_aegt_mn <- MCMCvis::MCMCpstr(si_fit, 
                                params = c('alpha', 'eta', 'gamma', 'theta'))
wi_aegt_mn <- MCMCvis::MCMCpstr(wi_fit, 
                                params = c('alpha', 'eta', 'gamma', 'theta'))

#mean and sd for cov
yr_center <- attr(scale(si_data$pro_data$year, scale = FALSE), 'scaled:center')
lat_cn_id <- unique(si_data$pro_data[,c('cn_id', 'lat')])
lat_center <- attr(scale(lat_cn_id$lat, scale = FALSE), 'scaled:center')
elev_cn_id <- unique(si_data$pro_data[,c('cn_id', 'GMTED_elev')])
elev_center <- attr(scale(elev_cn_id$GMTED_elev, scale = FALSE), 'scaled:center')

#simulate x data
L <- 25
sim_df <- data.frame(species = rep(NA, L * length(usp)),
                     syear = NA,
                     year = NA,
                     si_sim_year = NA,
                     wi_sim_year = NA,
                     slat = NA,
                     lat = NA,
                     si_sim_lat = NA,
                     wi_sim_lat = NA,
                     selev = NA,
                     elev = NA,
                     si_sim_elev = NA,
                     wi_sim_elev = NA)

counter <- 1
for (i in 1:length(usp))
{
  #i <- 1
  #at obs
  sp_idx1 <- which(si_data$sp == i)
  #at cn_id
  sp_idx2 <- which(si_data$cn_sp == i)
  
  sc2_year <- si_data$year[sp_idx1]
  sc2_lat <- si_data$lat[sp_idx2]
  sc2_elev <- si_data$elev[sp_idx2]
  
  mn_yr <- mean(sc2_year)
  mn_lat <- mean(sc2_lat)
  mn_elev <- mean(sc2_elev)
  
  SYR <- seq(range(sc2_year)[1], range(sc2_year)[2], length.out = L)
  SLAT <- seq(range(sc2_lat)[1], range(sc2_lat)[2], length.out = L)
  SELEV <- seq(range(sc2_elev)[1], range(sc2_elev)[2], length.out = L)
  
  sim_df$species[counter:(counter + L - 1)] <- usp[i]
  
  #alpha + eta * syr + gamma * mn_lat + theta * mn_elev
  sim_df$syear[counter:(counter + L - 1)] <- SYR
  sim_df$year[counter:(counter + L - 1)] <- (SYR * si_data$scf_YR) + yr_center
  sim_df$si_sim_year[counter:(counter + L - 1)] <- si_aegt_mn$alpha[i] + si_aegt_mn$eta[i] * SYR +
    si_aegt_mn$gamma[i] * mn_lat + si_aegt_mn$theta[i] * mn_elev
  sim_df$wi_sim_year[counter:(counter + L - 1)] <- wi_aegt_mn$alpha[i] + wi_aegt_mn$eta[i] * SYR +
    wi_aegt_mn$gamma[i] * mn_lat + wi_aegt_mn$theta[i] * mn_elev
  
  #alpha + gamma * slat + eta * mn_yr + theta * mn_elev
  sim_df$slat[counter:(counter + L - 1)] <- SLAT
  sim_df$lat[counter:(counter + L - 1)] <- (SLAT * si_data$scf_LAT) + lat_center
  sim_df$si_sim_lat[counter:(counter + L - 1)] <- si_aegt_mn$alpha[i] + si_aegt_mn$gamma[i] * SLAT +
    si_aegt_mn$eta[i] * mn_yr + si_aegt_mn$theta[i] * mn_elev
  sim_df$wi_sim_lat[counter:(counter + L - 1)] <- wi_aegt_mn$alpha[i] + wi_aegt_mn$gamma[i] * SLAT +
    wi_aegt_mn$eta[i] * mn_yr + wi_aegt_mn$theta[i] * mn_elev
  
  #alpha + theta * selev + eta * mn_yr + gamma * mn_lat
  sim_df$selev[counter:(counter + L - 1)] <- SELEV
  sim_df$elev[counter:(counter + L - 1)] <- (SELEV * si_data$scf_ELEV) + elev_center
  sim_df$si_sim_elev[counter:(counter + L - 1)] <- si_aegt_mn$alpha[i] + si_aegt_mn$theta[i] * SELEV +
    si_aegt_mn$eta[i] * mn_yr + si_aegt_mn$gamma[i] * mn_lat
  sim_df$wi_sim_elev[counter:(counter + L - 1)] <- wi_aegt_mn$alpha[i] + wi_aegt_mn$theta[i] * SELEV +
    wi_aegt_mn$eta[i] * mn_yr + wi_aegt_mn$gamma[i] * mn_lat
  
  counter <- counter + L
}


# simulate SI and WI over time,  lat,  elev (mean response) -------------------------------

#si - posterior means
si_mu_aegt_mn <- MCMCvis::MCMCpstr(si_fit, 
                                   params = c('mu_alpha', 'mu_eta', 'mu_gamma', 'mu_theta'))
wi_mu_aegt_mn <- MCMCvis::MCMCpstr(wi_fit, 
                                   params = c('mu_alpha', 'mu_eta', 'mu_gamma', 'mu_theta'))

#mean and sd for cov
yr_center <- attr(scale(si_data$pro_data$year, scale = FALSE), 'scaled:center')
lat_cn_id <- unique(si_data$pro_data[,c('cn_id', 'lat')])
lat_center <- attr(scale(lat_cn_id$lat, scale = FALSE), 'scaled:center')
elev_cn_id <- unique(si_data$pro_data[,c('cn_id', 'GMTED_elev')])
elev_center <- attr(scale(elev_cn_id$GMTED_elev, scale = FALSE), 'scaled:center')

#simulate x data
L <- 25
sim_mu_df <- data.frame(syear = rep(NA, L),
                        year = NA,
                        si_sim_year = NA,
                        wi_sim_year = NA,
                        slat = NA,
                        lat = NA,
                        si_sim_lat = NA,
                        wi_sim_lat = NA,
                        selev = NA,
                        elev = NA,
                        si_sim_elev = NA,
                        wi_sim_elev = NA)

sc2_year <- si_data$year
sc2_lat <- si_data$lat
sc2_elev <- si_data$elev

mn_yr <- 0
mn_lat <- 0
mn_elev <- 0

SYR <- seq(range(sc2_year)[1], range(sc2_year)[2], length.out = L)
SLAT <- seq(range(sc2_lat)[1], range(sc2_lat)[2], length.out = L)
SELEV <- seq(range(sc2_elev)[1], range(sc2_elev)[2], length.out = L)

#mu_alpha + mu_eta * syr + mu_gamma * mn_lat + mu_theta * mn_elev
sim_mu_df$syear <- SYR
sim_mu_df$year <- (SYR * si_data$scf_YR) + yr_center
sim_mu_df$si_sim_year <- si_mu_aegt_mn$mu_alpha + 
  si_mu_aegt_mn$mu_eta * SYR +
  si_mu_aegt_mn$mu_gamma * mn_lat + 
  si_mu_aegt_mn$mu_theta * mn_elev
sim_mu_df$wi_sim_year <- wi_mu_aegt_mn$mu_alpha + 
  wi_mu_aegt_mn$mu_eta * SYR +
  wi_mu_aegt_mn$mu_gamma * mn_lat + 
  wi_mu_aegt_mn$mu_theta * mn_elev

#mu_alpha + mu_gamma * slat + mu_eta * mn_yr + mu_theta * mn_elev
sim_mu_df$slat <- SLAT
sim_mu_df$lat <- (SLAT * si_data$scf_LAT) + lat_center
sim_mu_df$si_sim_lat <- si_mu_aegt_mn$mu_alpha + 
  si_mu_aegt_mn$mu_gamma * SLAT +
  si_mu_aegt_mn$mu_eta * mn_yr + 
  si_mu_aegt_mn$mu_theta * mn_elev
sim_mu_df$wi_sim_lat <- wi_mu_aegt_mn$mu_alpha + 
  wi_mu_aegt_mn$mu_gamma * SLAT +
  wi_mu_aegt_mn$mu_eta * mn_yr + 
  wi_mu_aegt_mn$mu_theta * mn_elev

#mu_alpha + mu_theta * selev + mu_eta * mn_yr + mu_gamma * mn_lat
sim_mu_df$selev <- SELEV
sim_mu_df$elev <- (SELEV * si_data$scf_ELEV) + elev_center
sim_mu_df$si_sim_elev <- si_mu_aegt_mn$mu_alpha + 
  si_mu_aegt_mn$mu_theta * SELEV +
  si_mu_aegt_mn$mu_eta * mn_yr + 
  si_mu_aegt_mn$mu_gamma * mn_lat
sim_mu_df$wi_sim_elev <- wi_mu_aegt_mn$mu_alpha + 
  wi_mu_aegt_mn$mu_theta * SELEV +
  wi_mu_aegt_mn$mu_eta * mn_yr + 
  wi_mu_aegt_mn$mu_gamma * mn_lat


# simulate si ~ temp for plotting -----------------------------------------

#lines and CI as ribbon for each temp lag
mu_kappa_l0_ch <- MCMCvis::MCMCchains(si_temp_l0_fit, params = 'mu_kappa')
mu_gamma_l0_ch <- MCMCvis::MCMCchains(si_temp_l0_fit, params = 'mu_gamma')
mu_kappa_l1_ch <- MCMCvis::MCMCchains(si_temp_l1_fit, params = 'mu_kappa')
mu_gamma_l1_ch <- MCMCvis::MCMCchains(si_temp_l1_fit, params = 'mu_gamma')
mu_kappa_l2_ch <- MCMCvis::MCMCchains(si_temp_l2_fit, params = 'mu_kappa')
mu_gamma_l2_ch <- MCMCvis::MCMCchains(si_temp_l2_fit, params = 'mu_gamma')

#sim temp
L <- 25

#get range of temp anomalies
si_temp_l0_data$pro_data %>%
  group_by(cn_id) %>%
  summarize(min_t = min(scale(MJJ_tmax, scale = FALSE)[,1] / si_temp_l0_data$scf_temp),
            max_t = max(scale(MJJ_tmax, scale = FALSE)[,1] / si_temp_l0_data$scf_temp)) -> l0_rng

si_temp_l1_data$pro_data %>%
  group_by(cn_id) %>%
  summarize(min_t = min(scale(MJJ_tmax, scale = FALSE)[,1] / si_temp_l0_data$scf_temp),
            max_t = max(scale(MJJ_tmax, scale = FALSE)[,1]) / si_temp_l0_data$scf_temp) -> l1_rng

si_temp_l2_data$pro_data %>%
  group_by(cn_id) %>%
  summarize(min_t = min(scale(MJJ_tmax, scale = FALSE)[,1] / si_temp_l0_data$scf_temp),
            max_t = max(scale(MJJ_tmax, scale = FALSE)[,1]) / si_temp_l0_data$scf_temp) -> l2_rng

tc <- c(min(l0_rng$min_t), max(l0_rng$max_t), 
        min(l1_rng$min_t), max(l1_rng$max_t),
        min(l2_rng$min_t), max(l2_rng$max_t))

sim_t <- seq(range(tc)[1], range(tc)[2], length = L)

#simulate x data
tanom_sim_df <- data.frame(stanom = sim_t,
                           tanom = sim_t * si_temp_l0_data$scf_temp,
                           l0_sim_mn = NA,
                           l0_sim_LCI = NA,
                           l0_sim_UCI = NA,
                           l1_sim_mn = NA,
                           l1_sim_LCI = NA,
                           l1_sim_UCI = NA,
                           l2_sim_mn = NA,
                           l2_sim_LCI = NA,
                           l2_sim_UCI = NA)

for (i in 1:NROW(tanom_sim_df))
{
  #i <- 1
  
  #mu_kappa + mu_gamma + mu_theta * sc_temp_anom
  pp_l0 <- mu_kappa_l0_ch[,1] + mu_gamma_l0_ch[,1] * tanom_sim_df$stanom[i]
  pp_l1 <- mu_kappa_l1_ch[,1] + mu_gamma_l1_ch[,1] * tanom_sim_df$stanom[i]
  pp_l2 <- mu_kappa_l2_ch[,1] + mu_gamma_l2_ch[,1] * tanom_sim_df$stanom[i]
  
  #L0
  tanom_sim_df$l0_sim_mn[i] <- mean(pp_l0)
  l0_qnt <- quantile(pp_l0, probs = c(0.055, 0.945))
  tanom_sim_df$l0_sim_LCI[i] <- l0_qnt[1]
  tanom_sim_df$l0_sim_UCI[i] <- l0_qnt[2]
  
  #L1
  tanom_sim_df$l1_sim_mn[i] <- mean(pp_l1)
  l1_qnt <- quantile(pp_l1, probs = c(0.055, 0.945))
  tanom_sim_df$l1_sim_LCI[i] <- l1_qnt[1]
  tanom_sim_df$l1_sim_UCI[i] <- l1_qnt[2]
  
  #L2
  tanom_sim_df$l2_sim_mn[i] <- mean(pp_l2)
  l2_qnt <- quantile(pp_l2, probs = c(0.055, 0.945))
  tanom_sim_df$l2_sim_LCI[i] <- l2_qnt[1]
  tanom_sim_df$l2_sim_UCI[i] <- l2_qnt[2]
}


# scaled posteriors for plotting ------------------------------------------

#posterior estimates for SI and WI trends over time, lat, elev

#function to scale back to years, degrees lat, and m
#args to scale to X years, X degrees lat, and X m
sc_fun <- function(input, syr = 10, slat = 10, selev = 1000)
{
  si_ch <- MCMCvis::MCMCchains(input, params = c('mu_eta', 'mu_gamma', 'mu_theta'))
  si_sc_fit <- cbind(mu_eta = (si_ch[,'mu_eta'] / si_data$scf_YR) * syr, 
                     mu_gamma = (si_ch[,'mu_gamma'] / si_data$scf_LAT) * slat,
                     mu_theta = (si_ch[,'mu_theta'] / si_data$scf_ELEV) * selev)#,
  return(si_sc_fit)
}

#apply to si - 10 years, 10 degrees lat, 1km elev
si_sc <- sc_fun(si_fit, syr = 10, slat = 10, selev = 1000)

#apply to wi
wi_sc <- sc_fun(wi_fit, syr = 10, slat = 10, selev = 1000)

#scale for each species
sc_fun2 <- function(input, param, syr = 10, slat = 10, selev = 1000)
{
  if (param == 'eta')
  {
    tt_eta_ch <- (MCMCvis::MCMCchains(input, params = 'eta') / si_data$scf_YR) * syr
    return(tt_eta_ch)
  }
  if (param == 'gamma')
  {
    tt_gamma_ch <- (MCMCvis::MCMCchains(input, params = 'gamma') / si_data$scf_LAT) * slat
    return(tt_gamma_ch)
  }
  if (param == 'theta')
  {
    tt_theta_ch <- (MCMCvis::MCMCchains(input, params = 'theta') / si_data$scf_ELEV) * selev
    return(tt_theta_ch)
  }
}

si_eta_sc <- sc_fun2(si_fit, param = 'eta')
si_gamma_sc <- sc_fun2(si_fit, param = 'gamma')
si_theta_sc <- sc_fun2(si_fit, param = 'theta')
wi_eta_sc <- sc_fun2(wi_fit, param = 'eta')
wi_gamma_sc <- sc_fun2(wi_fit, param = 'gamma')
wi_theta_sc <- sc_fun2(wi_fit, param = 'theta')


# SI ~ time --------------------------------------------

#plotting params - size ~ time
#line size
SZ <- 2
#text size
TSZ <- 14
#text title size
TTSZ <- 18
#bounding box size
LSZ <- 2
#tick size
TISZ <- 1.5
#transparency
#ALPHA <- 0.15
ALPHA <- 0.25
#colors
COLOR_SI <- 'black'
COLOR_WI <- 'black'

#si ~ year
si_year <- ggplot(sim_df, aes(year, si_sim_year, group = factor(species))) +
  geom_line(alpha = ALPHA, size = SZ, color = 'grey') +
  geom_line(data = sim_mu_df, inherit.aes = FALSE, 
            aes(year, si_sim_year), alpha = 0.7, size = SZ + 1.5, color = COLOR_SI) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        legend.position = 'none') +
  xlab('Year') +
  ylab('Size Index')

#plot figs
# SI ~ year
ggsave(paste0(fig_dir, 'si_year_lines.pdf'), si_year, height = 4, width = 4)

#species-specific estimates as points, overall estimate as caterpillar
pdf(paste0(fig_dir, 'si_year_cat.pdf'), height = 3, width = 3)
MCMCvis::MCMCplot(si_sc,
                  params = 'mu_eta',
                  col = COLOR_SI,
                  ci = c(50, 89),
                  labels = NULL,
                  sz_med = 2,
                  sz_thick = 7,
                  sz_thin = 5,
                  xlim = c(-0.15, 0.1))

#species-specific estimates as circles, overall effect as caterpillar
pmed_si_eta <- apply(si_eta_sc, 2, median)
yv <- rep(1, length(pmed_si_eta)) + runif(length(pmed_si_eta), -0.2, 0.2)
points(pmed_si_eta, yv, 
       col = alpha('grey', alpha = 0.1), cex = 1.4, pch = 19)
dev.off()


# SI ~ temp (temporal) ---------------------------------------------

#plotting params - size ~ temp
#line transparency
L_ALPHA <- 0.8
#ribbon transparency
R_ALPHA <- 0.2
#line colors
l0_col <- '#FFC30F'
l1_col <- '#C70039'
l2_col <- '#581845'
#line size
SZ <- 2

#create plot
tanom_plt <- ggplot(tanom_sim_df, aes(tanom, l0_sim_mn)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ)) +
  geom_line(size = SZ, alpha = L_ALPHA, col = l0_col) +
  geom_ribbon(aes(ymin = l0_sim_LCI, ymax = l0_sim_UCI),
              alpha = R_ALPHA, fill = l0_col) +
  geom_line(aes(tanom, l1_sim_mn), col = l1_col,
            size = SZ, alpha = L_ALPHA) +
  geom_ribbon(aes(ymin = l1_sim_LCI, ymax = l1_sim_UCI),
              alpha = R_ALPHA, fill = l1_col) +
  geom_line(aes(tanom, l2_sim_mn), col = l2_col,
            size = SZ, alpha = L_ALPHA) +
  geom_ribbon(aes(ymin = l2_sim_LCI, ymax = l2_sim_UCI),
              alpha = R_ALPHA, fill = l2_col) +
  xlab('Temp anomaly') +
  ylab('Mean SI')


#save RDS
# SI ~ temp
ggsave(paste0(fig_dir, 'si_temp_times_lines.pdf'), tanom_plt, height = 4, width = 4)

#legend
pdf(paste0(fig_dir, 'si_temp_times_lines.pdf_legend.pdf'), height = 3, width = 3)
plot(1:10, col = 'white')
abline(h = 2, col = l0_col, lwd = 3)
abline(h = 4, col = l1_col, lwd = 3)
abline(h = 6, col = l2_col, lwd = 3)
dev.off()


# SI ~ lat ---------------------------------------------------------------

#si ~ lat
si_lat <- ggplot(sim_df, aes(lat, si_sim_lat, group = factor(species))) +
  geom_line(alpha = ALPHA, size = SZ, color = 'grey') +
  geom_line(data = sim_mu_df, inherit.aes = FALSE, 
            aes(lat, si_sim_lat), alpha = 0.7, size = SZ + 1.5, color = COLOR_SI) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        legend.position = 'none') +
  xlab('Latitude') +
  ylab('Size Index')

#plot figs
# SI ~ year
ggsave(paste0(fig_dir, 'si_lat_lines.pdf'), si_lat, height = 4, width = 4)

#species-specific estimates as points, overall estimate as caterpillar
pdf(paste0(fig_dir, 'si_lat_cat.pdf'), height = 3, width = 3)
MCMCvis::MCMCplot(si_sc,
                  params = 'mu_gamma',
                  col = COLOR_SI,
                  ci = c(50, 89),
                  labels = NULL,
                  sz_med = 2,
                  sz_thick = 7,
                  sz_thin = 5,
                  xlim = c(-0.6, 2))
#species-specific estimates as circles, overall effect as caterpillar
pmed_si_gamma <- apply(si_gamma_sc, 2, median)
yv <- rep(1, length(pmed_si_gamma)) + runif(length(pmed_si_gamma), -0.2, 0.2)
points(pmed_si_gamma, yv, 
       col = alpha(COLOR_SI, alpha = 0.1), cex = 1.4, pch = 19)
dev.off()


# SI ~ temp (spatial) ---------------------------------------------

#to get absolute temp not anom
setwd(paste0(dir, '/Data/daymet/processed/'))
daymet <- readRDS(paste0('daymet-master-', daymet_process_date, '.rds')) %>%
  mutate(MJJ_tmax = (May_tmax + June_tmax + July_tmax) / 3) #using mean here does not do what might be expected
si_tsd <- dplyr::select(si_temp_space_data$pro_data, sci_name, sp_id, station, cn_id, lat, 
                        lng, year, GMTED_elev, size_idx) %>%
  dplyr::left_join(daymet, by = c('station', 'year')) %>%
  dplyr::arrange(sp_id, cn_id) %>%
  dplyr::group_by(sp_id, cn_id) %>%
  dplyr::summarize(mn_MJJ_tmax = mean(MJJ_tmax)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sp_id) %>%
  dplyr::summarize(smt = mean(mn_MJJ_tmax))
mean(si_tsd$smt)
range(si_tsd$smt)

#change to absolute temp experienced by species
#effect of temp ~ sp_mean_temp
tplt <- data.frame(y = si_temp_sp_cov_data$y,
                   sd_y = si_temp_sp_cov_data$sd_y,
                   mn_temp = si_temp_sp_cov_data$mn_temp,
                   mn_temp_sc = si_temp_sp_cov_data$mn_temp * 10,
                   mn_temp_act = (si_temp_sp_cov_data$mn_temp * 10) + mean(si_tsd$smt))

#int
alpha_ch <- MCMCvis::MCMCchains(si_temp_sp_cov_fit, params = 'alpha')
#slope
beta_ch <- MCMCvis::MCMCchains(si_temp_sp_cov_fit, params = 'beta')

xsim <- seq(range(si_temp_sp_cov_data$mn_temp)[1], range(si_temp_sp_cov_data$mn_temp)[2], length.out = 50)
mf <- matrix(NA, nrow = NROW(alpha_ch), ncol = length(xsim))
for (i in 1:NROW(alpha_ch))
{
  #i <- 1
  mf[i,] <- alpha_ch[i,] + beta_ch[i,] * xsim
}

fplt <- data.frame(xsim,
                   xsim_sc = xsim * 10,
                   xsim_act = (xsim * 10) + mean(si_tsd$smt),
                   mn_fit = apply(mf, 2, mean),
                   lci_fit = apply(mf, 2, function(x) quantile(x, probs = 0.025)),
                   uci_fit = apply(mf, 2, function(x) quantile(x, probs = 0.975)))

p1 <- ggplot(tplt, aes(mn_temp_act, y)) +
  geom_point(size = 4, alpha = 0.35) +
  geom_errorbar(aes(x = mn_temp_act, ymin = y - sd_y, 
                    ymax = y + sd_y),
                alpha = 0.4) +
  theme_bw() +
  theme(axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
  ylab('Effect of temp on SI over space') +
  xlab('Mean temp experienced by species')

#add grey box instead of line
p2 <- p1 + geom_hline(yintercept = 0, lty = 2, size = 2, col = 'red')

p3 <- p1 + geom_ribbon(data = fplt, inherit.aes = FALSE,
                       aes(x = xsim_act, 
                           ymax = uci_fit, ymin = lci_fit),
                       alpha = 0.25) +
  geom_line(data = fplt, inherit.aes = FALSE,
            aes(x = xsim_act, y = mn_fit), color = 'red',
            alpha = 0.6, size = 2.5)

# #plot figs
ggsave(filename = paste0(fig_dir, 'si_temp_space_mn_sp_temp_p1.pdf'), p1, height = 5, width = 5) 
ggsave(filename = paste0(fig_dir, 'si_temp_space_mn_sp_temp_p2.pdf'), p2, height = 5, width = 5) 
ggsave(filename = paste0(fig_dir, 'si_temp_space_mn_sp_temp_p3.pdf'), p3, height = 5, width = 5) 


# SI ~ elev ---------------------------------------------------------------

#si ~ elev
si_elev <- ggplot(sim_df, aes(elev, si_sim_elev, group = factor(species))) +
  geom_line(alpha = ALPHA, size = SZ, color = 'grey') +
  geom_line(data = sim_mu_df, inherit.aes = FALSE, 
            aes(elev, si_sim_elev), alpha = 0.7, size = SZ + 1.5, color = COLOR_SI) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        legend.position = 'none') +
  xlab('Elevation') +
  ylab('Size Index')

#plot figs
# SI ~ year
ggsave(paste0(fig_dir, 'si_elev_lines.pdf'), si_elev, height = 4, width = 4)

#species-specific estimates as points, overall estimate as caterpillar
pdf(paste0(fig_dir, 'si_elev_cat.pdf'), height = 3, width = 3)
MCMCvis::MCMCplot(si_sc,
                  params = 'mu_theta',
                  col = COLOR_SI,
                  ci = c(50, 89),
                  labels = NULL,
                  sz_med = 2,
                  sz_thick = 7,
                  sz_thin = 5,
                  xlim = c(-1, 1.5))
#species-specific estimates as circles, overall effect as caterpillar
pmed_si_theta <- apply(si_theta_sc, 2, median)
yv <- rep(1, length(pmed_si_theta)) + runif(length(pmed_si_theta), -0.2, 0.2)
points(pmed_si_theta, yv, 
       col = alpha(COLOR_SI, alpha = 0.1), cex = 1.4, pch = 19)
dev.off()


# WI ~ elev ---------------------------------------------------------------

wi_elev <- ggplot(sim_df, aes(elev, wi_sim_elev, group = factor(species))) +
  geom_line(alpha = ALPHA, size = SZ, color = 'grey') +
  geom_line(data = sim_mu_df, inherit.aes = FALSE, 
            aes(elev, wi_sim_elev), alpha = 0.7, size = SZ + 1.5, color = COLOR_WI) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = LSZ),
        axis.ticks = element_line(size = TISZ),
        axis.text.x = element_text(size = TSZ),
        axis.title.x = element_text(size = TTSZ),
        axis.text.y = element_text(size = TSZ),
        axis.title.y = element_text(size = TTSZ),
        legend.position = 'none') +
  xlab('Elevation') +
  ylab('Wing Index')

#plot figs
# SI ~ year
ggsave(paste0(fig_dir, 'wi_elev_lines.pdf'), wi_elev, height = 4, width = 4)

#species-specific estimates as points, overall estimate as caterpillar
pdf(paste0(fig_dir, 'wi_elev_cat.pdf'), height = 3, width = 3)
MCMCvis::MCMCplot(wi_sc,
                  params = 'mu_theta',
                  col = COLOR_WI,
                  ci = c(50, 89),
                  labels = NULL,
                  sz_med = 2,
                  sz_thick = 7,
                  sz_thin = 5,
                  xlim = c(-0.2, 1))
#species-specific estimates as circles, overall effect as caterpillar
pmed_wi_theta <- apply(wi_theta_sc, 2, median)
yv <- rep(1, length(pmed_wi_theta)) + runif(length(pmed_wi_theta), -0.2, 0.2)
points(pmed_wi_theta, yv, 
       col = alpha(COLOR_WI, alpha = 0.1), cex = 1.4, pch = 19)
dev.off()


# temporal coverage data --------------------------------------------------

#lat/lon for distinct stations
mdata %>%
  dplyr::distinct(station, lat, lng, GMTED_elev, year) -> ustyr

#length of line
LL <- 1
plt2 <- ustyr %>%
  mutate(x = year - (LL/2), xend = year + (LL/2), yend = factor(lat))

#sort by LL
temp_coverage <- ggplot(data = plt2) + 
  geom_segment(aes(x = x, xend = xend, y = yend, yend = yend)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks.x = element_line(size = 1.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 20),
        legend.position = 'none') +
  xlab('Year') +
  ylab('Station')

ggsave(paste0(fig_dir, 'temp_coverage.pdf'), temp_coverage, height = 8, width = 8)


# hillshade map of study area ---------------------------------------------------------------------

ustyr %>%
  dplyr::group_by(station) %>%
  dplyr::summarize(year = n(),
                   lat = mean(lat),
                   lng = mean(lng),
                   GMTED_elev = mean(GMTED_elev)) -> cll

#get elev data and merge
alt_US <- getData('alt', country = 'USA', path = paste0(dir, 'Data'))
alt_Can <- getData('alt', country = 'Canada', path = paste0(dir, 'Data'))
alt_Mex <- getData('alt', country = 'Mexico', path = paste0(dir, 'Data'))
alt_Gr <- getData('alt', country = 'Greenland', path = paste0(dir, 'Data'))
alt_mrg <- raster::merge(alt_US[[1]], alt_US[[2]], alt_Can, alt_Mex, 
                         alt_Gr)

#coarsen
alt_mrg2 <- raster::aggregate(alt_mrg, fact = 12)

#exaggerate elevation for relief
alt_mrg3 <- alt_mrg2 * 30

#get hillshade
slope <- raster::terrain(alt_mrg3, opt = 'slope')
aspect <- raster::terrain(alt_mrg3, opt = 'aspect')
hill <- raster::hillShade(slope, aspect, 50, 270)

#color palette taken from map here: https://gis.stackexchange.com/questions/25099/choosing-colour-ramp-to-use-for-elevation
#using this tool: https://imagecolorpicker.com/en
cp <- c('#649847', '#e2dca4', '#fdc66c', '#bb9f8b', '#f2eeeb')

cll_sf <- sf::st_as_sf(cll, coords = c('lng', 'lat'))

#Lambert Conformal Conic proj4string
#https://proj.org/operations/projections/lcc.html
lcc <- '+proj=lcc +lon_0=-90 +lat_0=60 +lat_1=33 +lat_2=45'

#project to lcc
hill_lcc <- raster::projectRaster(hill, crs = sp::CRS(lcc))
alt_mrg2_lcc <- raster::projectRaster(alt_mrg2, crs = sp::CRS(lcc))
cll_sf_lcc <- st_set_crs(cll_sf, value = 4326) %>% 
  st_transform(crs = st_crs(lcc))

#xmin, ymin, xmax, ymax
#works: -150, -50, 15, 75
bb_lcc <- sf::st_bbox(c(xmin = -150, xmax = -50, ymin = 15, ymax = 75), crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  st_transform(crs = st_crs(lcc)) %>%
  st_bbox()

#plot hillshade
sa_map_legend <- tm_shape(hill_lcc, bbox = bb_lcc) +
  tm_raster(palette = gray(0:100 / 100), n = 100, legend.show = FALSE)  +
  tm_shape(alt_mrg2_lcc) +
  tm_raster(alpha = 0.6, 
            breaks = c(-200, 500, 1000, 2000, 3000, 10000),
            palette = cp,
            legend.show = TRUE) + 
  tm_shape(cll_sf_lcc) +
  tm_dots(size = 0.15, alpha = 0.4, col = 'black') +
  tm_graticules(alpha = 0.7) +
  tm_layout(legend.title.size = 1,
            legend.text.size = 0.6,
            legend.position = c("left","bottom"),
            legend.bg.color = "white",
            legend.bg.alpha = 1)

sa_map <- sa_map_legend + 
  tm_layout(legend.show = FALSE)


#import into illustrator: 
#-move legend around
tmap_save(sa_map, filename = paste0(fig_dir, 'sa_map.pdf'), width = 8, height = 8)
tmap_save(sa_map_legend, filename = paste0(fig_dir, 'sa_map_legend.pdf'), width = 8, height = 8)


# SI, WI, and temp maps ---------------------------------------------------------

#SI and WI over space - response of SI to temp over space
si_alpha_ch <- MCMCvis::MCMCchains(si_fit, params = 'alpha')
si_gamma_ch <- MCMCvis::MCMCchains(si_fit, params = 'gamma')
si_theta_ch <- MCMCvis::MCMCchains(si_fit, params = 'theta')

wi_alpha_ch <- MCMCvis::MCMCchains(wi_fit, params = 'alpha')
wi_gamma_ch <- MCMCvis::MCMCchains(wi_fit, params = 'gamma')
wi_theta_ch <- MCMCvis::MCMCchains(wi_fit, params = 'theta')

#add mn station temp and effect of temp on SI
si_temp_l0_data$pro_data %>%
  dplyr::group_by(cn_id) %>%
  dplyr::summarize(mn_st_temp = mean(MJJ_tmax)) %>%
  dplyr::ungroup() -> mn_t

dplyr::distinct(si_temp_l0_data$pro_data, sci_name, station, cn_id, lat, lng) %>%
  dplyr::left_join(mn_t) -> mn_t2

mn_t2$si_temp_beta_mn <- MCMCvis::MCMCpstr(si_temp_l0_fit, params = 'beta')[[1]] / si_temp_l0_data$scf_temp

#read in temperature data to get MJJ tmax for 2018 - for contours
dm_2018_M_tm <- raster::raster(paste0(dir, 'Data/daymet/RAW/daymet_v4_tmax_monavg_na_2018.tif'), band = 5)
dm_2018_June_tm <- raster::raster(paste0(dir, 'Data/daymet/RAW/daymet_v4_tmax_monavg_na_2018.tif'), band = 6)
dm_2018_July_tm <- raster::raster(paste0(dir, 'Data/daymet/RAW/daymet_v4_tmax_monavg_na_2018.tif'), band = 7)

dm_2018_MJJ_tm <- calc(raster::stack(c(dm_2018_M_tm, dm_2018_June_tm, 
                                       dm_2018_July_tm)), fun = mean)
#coarsen to reduce plot time
dm_2018_MJJ_tm_ag <- raster::aggregate(dm_2018_MJJ_tm, fact = 100)

#need to include underscore between genus and species name
#creates pdf of early year idx (SI or WI or temp), late year idx, and plot for legend
map_fun <- function(species, out_dir, idx = 'WI', station = FALSE,
                    orientation, xlim, ylim, N = 1e7, PT_SZ = 0.5, 
                    ST_ALPHA = 0.5, ST_PT_SZ = 1, WW = 8, HH = 8, MN = TRUE, 
                    elev_limit = NA, COARSEN = 25, low_col = '#C7522B', 
                    mid_col = '#FBF2C4', high_col = '#3C5941')
{
  #species <- 'Setophaga_americana'
  #reference key for species synonyms
  setwd(paste0(dir, '/Data/BirdLife_range_maps/metadata'))
  sp_key <- read.csv('species_filenames_key.csv')
  
  if (idx == 'SI')
  {
    alpha_ch <- si_alpha_ch
    gamma_ch <- si_gamma_ch
    theta_ch <- si_theta_ch
    data <- si_data
  }
  if (idx == 'WI')
  {
    alpha_ch <- wi_alpha_ch
    gamma_ch <- wi_gamma_ch
    theta_ch <- wi_theta_ch
    data <- wi_data
  }
  if (idx == 'temp')
  {
    data <- si_data
  }
  
  td <- dplyr::filter(data$pro_data, sci_name == gsub('_', ' ', species))
  sp_idx <- td$sp_id[1]
  
  print('Parsing range')
  #change dir to shp files
  setwd(paste0(dir, '/Data/BirdLife_range_maps/shapefiles/'))
  
  #filter by breeding/migration cells
  #match species name to shp file name
  g_ind <- grep(species, sp_key$file_names_2016)
  
  #check for synonyms if there are no matches
  if (length(g_ind) == 0)
  {
    g_ind2 <- grep(species, sp_key$BL_Checklist_name)
  } else {
    g_ind2 <- g_ind
  }
  
  if (length(g_ind2) > 0)
  {
    
    #get filename and read in
    fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', 
                                                         sp_key[g_ind2, 'filenames'])]) 
    spdf_rng <- rgdal::readOGR(fname[1], verbose = FALSE)
    
    #filter by resident (1) and breeding range (2) - need to convert spdf to sp
    spdf_brng <- spdf_rng[which(spdf_rng$SEASONAL == 1 | spdf_rng$SEASONAL == 2),]
    
    #plotting species breeding range
    spdf_brng@data$id <- rownames(spdf_brng@data)
    #add buffer to deal with intersection warnings
    spdf_brng2 <- rgeos::gBuffer(spdf_brng, width = 0)
    spdf_brng.points <- ggplot2::fortify(spdf_brng2, region = "id")
    spdf_brng.df <- plyr::join(spdf_brng.points, spdf_brng@data, by = "id")
    #convert spdf to sp
    sp_brng <- sp::SpatialPolygons(spdf_brng@polygons)
    sp::proj4string(sp_brng) <- sp::CRS('+init=epsg:4326')
    
    if (idx == 'SI' | idx == 'WI')
    {
      # extract elevation at points 
      setwd(paste0(dir, 'Data/GMTED2010_mean/'))
      
      print('Extracting elevation')
      # read in mosaic tif
      elev_mosaic <- raster::raster('elev_mosaic.tif')
      
      #coarsen, mask, and crop elev raster
      #make sure only NA
      ext <- raster::extent(sp_brng)
      if (ext[3] < 25)
      {
        ext[3] <- 25
      }
      elev_coarse <- raster::crop(raster::mask(raster::aggregate(elev_mosaic, fact = COARSEN), sp_brng, inverse = FALSE), ext)
      pp_rast <- elev_coarse
      
      #plot(pp_rast)
      pp_rast[,] <- NA
      
      #[lon, lat]
      r_lat <- (raster::coordinates(pp_rast)[,2] -  mean(data$pro_data$lat)) / data$scf_LAT
      r_elev <- (raster::getValues(elev_coarse) - mean(data$pro_data$GMTED_elev)) / data$scf_ELEV
      
      print('Predicting idx')
      pp <- data.frame(lat = r_lat, 
                       r_elev = r_elev,
                       mn_pred = NA)
      for (i in 1:length(r_lat))
      {
        #i <- 119431
        print(paste0('processing: ', i, ' of ', length(r_lat)))
        
        t_pred <- alpha_ch[,sp_idx] +
          (gamma_ch[,sp_idx] * r_lat[i]) +
          (theta_ch[,sp_idx] * r_elev[i])
        
        pp$mn_pred[i] <- mean(t_pred)
      }
      
      if (!is.na(elev_limit))
      {
        to.na <- which(raster::getValues(elev_coarse) > elev_limit)
        pp$mn_pred[to.na] <- NA
      }
      
      mn_rng <- mean(range(pp$mn_pred, na.rm = TRUE))

      # plot
      print('Creating plots')
      
      #load map
      worldmap <- ggplot2::map_data("world")
      
      if (MN == TRUE)
      {
        #range of values
        
        pp_rast[,] <- pp$mn_pred - mn_rng
        
        rng_pred <- range(raster::getValues(pp_rast), na.rm = TRUE)
        plt_rast <- data.frame(raster::rasterToPoints(pp_rast))
        
        morph_out <- ggplot() +
          #map
          geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), 
                       inherit.aes = FALSE, fill = alpha('black', 0.1), 
                       color = NA) +
          coord_map("lambert", lat0 = 60, lat1 = 33,
                    xlim = xlim, ylim = ylim,
                    orientation = c(90, 0, -90)) +
          theme_void() +
          #grid lines
          theme(panel.grid.major = element_line(color = alpha('black', 
                                                              0.2),
                                                size = 0.5),
                panel.ontop = FALSE,
                panel.background = element_rect(fill = '#FFFFFF', color = '#FFFFFF')) +
          geom_tile(data = plt_rast, aes(x, y, fill = elev_mosaic)) +
          scale_fill_gradient2(low = low_col, mid = mid_col, high = high_col,
                               limits = c(rng_pred[1], rng_pred[2]), midpoint = 0) +
          labs(fill = 'Predicted Idx') +
          ggtitle(paste0(species, ' predicted ', idx)) +
          xlab('Longitude') +
          ylab('Latitude')
      } else {
        #range of values
        pp_rast[,] <- pp$mn_pred
        
        rng_pred <- range(raster::getValues(pp_rast), na.rm = TRUE)
        plt_rast <- data.frame(raster::rasterToPoints(pp_rast))
        
        morph_out <- ggplot() +
          #map
          geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), 
                       inherit.aes = FALSE, fill = alpha('black', 0.1), 
                       color = NA) +
          coord_map("lambert", lat0 = 60, lat1 = 33,
                    xlim = xlim, ylim = ylim,
                    orientation = c(90, 0, -90)) +
          theme_void() +
          #grid lines
          theme(panel.grid.major = element_line(color = alpha('black', 
                                                              0.2),
                                                size = 0.5),
                panel.ontop = FALSE,
                panel.background = element_rect(fill = '#FFFFFF', color = '#FFFFFF')) +
          geom_tile(data = plt_rast, aes(x, y, fill = elev_mosaic)) +
          scale_fill_gradient2(low = low_col, mid = mid_col, high = high_col,
                               limits = c(rng_pred[1], rng_pred[2]), midpoint = 0) +
          labs(fill = 'Predicted Idx') +
          ggtitle(paste0(species, ' predicted ', idx)) +
          xlab('Longitude') +
          ylab('Latitude')
      }
      
      morph_out_nl <- morph_out +
        theme(legend.position = 'none')
      
      ggsave(paste0(out_dir, species, '_predicted_',idx, '.pdf'), morph_out_nl, width = WW, height = HH)
      ggsave(paste0(out_dir, species, '_predicted_', idx,' _w_legend.pdf'), morph_out, width = WW, height = HH)
    } else {
      #plot effect of temp on SI
      mn_t3 <- dplyr::filter(mn_t2, sci_name == gsub('_', ' ', species))
      
      #simplify breeding range
      spdf_brng_s <- rgeos::gSimplify(spdf_brng, tol = 0.1)
      
      #buffer to make sure contour lines aren't broken
      spdf_brng_sb <- rgeos::gBuffer(spdf_brng_s, width = 1)
      
      print('Creating plots')
      
      #reproject temp
      dm_2018_MJJ_tm_ag_tr <- raster::projectRaster(from = dm_2018_MJJ_tm_ag, 
                                                     crs = sp::CRS(sp::proj4string(spdf_brng_sb)))
      #mask over just breeding/resident range
      temp_msk <- raster::crop(raster::mask(dm_2018_MJJ_tm_ag_tr, spdf_brng_sb), extent(spdf_brng_sb))
      
      temp_mask_spdf <- as.data.frame(as(temp_msk, 'SpatialPixelsDataFrame'))
      colnames(temp_mask_spdf) <- c('value', 'x', 'y')
      
      worldmap <- ggplot2::map_data("world")
      
      morph_out <- ggplot() +
        #map
        geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), 
                     inherit.aes = FALSE, fill = alpha('black', 0.1), 
                     #color = alpha('black', 0.5)) +
                     color = NA) +
        #species range
        geom_polygon(data = spdf_brng_s, aes(long, lat, group = group),
                     inherit.aes = FALSE, fill = alpha('black', 1)) +
        coord_map("lambert", lat0 = 60, lat1 = 33,
                  xlim = xlim, ylim = ylim,
                  orientation = c(90, 0, -90)) +
        theme_void() +
        #grid lines
        theme(panel.grid.major = element_line(color = alpha('black', 
                                                            0.2),
                                              size = 0.5),
              panel.ontop = FALSE,
              panel.background = element_rect(fill = '#FFFFFF', color = '#FFFFFF')) +
        #points
        geom_point(data = mn_t3, aes(lng, lat, col= si_temp_beta_mn),
                   alpha = 0.9, inherit.aes = FALSE,
                   size = PT_SZ) +
        scale_color_gradient2(low = '#C7522B', mid = '#FBF2C4', high = '#3C5941',
                              limits = range(mn_t3$si_temp_beta_mn), midpoint = 0) +
        #temp contours
        stat_contour(data = temp_mask_spdf, aes(x, y, z = value), 
                     color = 'white', 
                     alpha = 1, size = 0.8,
                     breaks = c(0, 18, 22, 26, 30, 34, 50)) +
        ggtitle(paste0(species, ' - temp effect SI')) +
        xlab('Longitude') +
        ylab('Latitude')
      
      ggsave(paste0(out_dir, species, '_temp_effect.pdf'), morph_out, width = WW, height = HH)
    }
  } else {
    print(paste0('Range files not available for: ', species))
  }
}

#SI
#import into illustrator: 
#-remove grey background
map_fun(species = 'Vireo_olivaceus', 
        out_dir = fig_dir, #fig_dir
        xlim = c(-150, -50),
        ylim = c(26, 73),
        #IDX
        idx = 'SI',
        WW = 8,
        HH = 8,
        MN = FALSE, 
        elev_limit = 2000,
        COARSEN = 10, 
        low_col = '#354b99',
        mid_col = '#eaebcc',
        high_col = '#a50026')

#WI - increase over lat and elev
map_fun(species = 'Setophaga_americana',
        out_dir = fig_dir, #fig_dir
        xlim = c(-150, -50),
        ylim = c(26, 73),
        #IDX
        idx = 'WI',
        WW = 8,
        HH = 8,
        MN = FALSE,
        COARSEN = 10, 
        low_col = '#354b99',
        mid_col = '#eaebcc',
        high_col = '#a50026')


#TEMP
#import into illustrator: 
#-remove grey background
#-add text to isoclines
map_fun(species = 'Vireo_olivaceus', 
        out_dir = fig_dir, 
        xlim = c(-150, -50),
        ylim = c(26, 73),
        #temp
        idx = 'temp',
        PT_SZ = 3,
        WW = 8,
        HH = 8)


# morphological indices ---------------------------------------------------

set.seed(1)
N <- 25
x <- scale(runif(N, -10, 10), scale = FALSE)
eps <- rnorm(N, 0, 1)
y <- scale(0.33 * x + eps, scale = FALSE)

tplt <- data.frame(x, y)
p1 <- ggplot(tplt, aes(x, y)) +
  geom_point(size = 5) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
  ylim(range(y)[1] - 1, range(y)[2] + 1) +
  xlim(range(x)[1] - 1, range(x)[2] + 1)

#fit linear regression
ff <- summary(lm(y ~ x))
xsim <- seq(range(x)[1], range(x)[2], length.out = 20)
fplt <- data.frame(xsim = xsim,
                   fit = ff$coefficients[1,1] + xsim * ff$coefficients[2,1])

p2 <- p1 + geom_line(data = fplt, aes(xsim, fit),
                     color = 'red', alpha = 0.8, size = 3)

#apply rotation matrix
angle <- -atan(ff$coefficients[2,1])
rotm <- matrix(c(cos(angle), sin(angle), 
                 -sin(angle), cos(angle)), ncol = 2)

M <- as.data.frame(t(rotm %*% t(cbind(x, y))))
Mf <- as.data.frame(t(rotm %*% t(cbind(fplt$xsim, fplt$fit))))

p3 <- ggplot(M, aes(V1, V2)) +
  geom_point(size = 5) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
  ylim(range(y)[1] - 1, range(y)[2] + 1) +
  xlim(range(x)[1] - 1, range(x)[2] + 1) +
  geom_line(data = Mf, aes(V1, V2), 
            color = 'red', alpha = 0.8, size = 3)


ggsave(p1, filename = paste0(fig_dir, 'rotation1.pdf'), width = 5, height = 5)
ggsave(p2, filename = paste0(fig_dir, 'rotation2.pdf'), width = 5, height = 5)
ggsave(p3, filename = paste0(fig_dir, 'rotation3.pdf'), width = 5, height = 5)


# cross species wing ~ mass -----------------------------------------------

sp_mn_morph <- mdata %>%
  dplyr::distinct(sci_name, mean_l_weight, mean_l_wing_chord)

#log mass ~ log wing
lmlw_plt <- ggplot(sp_mn_morph, aes(mean_l_weight, mean_l_wing_chord)) +
  geom_point(alpha = 0.5, size = 4) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.ticks = element_line(size = 1.5),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) +
  xlab('log(mass)') +
  ylab('log(wing length)')


spp1 <- dplyr::filter(sp_mn_morph, sci_name == 'Vireo olivaceus')
spp2 <- dplyr::filter(sp_mn_morph, sci_name == 'Setophaga americana')

lmlw_plt + 
  geom_point(data = spp1, inherit.aes = FALSE, 
             aes(mean_l_weight, mean_l_wing_chord),
             col = 'red', size = 4) +
  geom_point(data = spp2, inherit.aes = FALSE, 
             aes(mean_l_weight, mean_l_wing_chord),
             col = 'green', size = 4) +  
  
  ggsave(lmlw_plt, filename = paste0(fig_dir, 'lmasslwing.pdf'), width = 5, height = 5)


# caterpillar plots - SI ~ TLE --------------------------------------------------

#10 years, 10 degrees lat, 1km elev
pdf(paste(fig_dir, 'si_year_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_eta_sc, 
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'SI ~ year',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-0.4, 0.6))
dev.off()

pdf(paste(fig_dir, 'si_lat_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_gamma_sc, 
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'SI ~ lat',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-1, 3))
dev.off()

pdf(paste(fig_dir, 'si_elev_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_theta_sc, 
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'SI ~ elev',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-2, 2))
dev.off()


# caterpillar plots - WI ~ TLE --------------------------------------------------

pdf(paste(fig_dir, 'wi_year_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(wi_eta_sc, 
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'WI ~ year',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-0.6, 0.4))
dev.off()

pdf(paste(fig_dir, 'wi_lat_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(wi_gamma_sc, 
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'WI ~ lat',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-1, 1))
dev.off()

pdf(paste(fig_dir, 'wi_elev_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(wi_theta_sc, 
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'WI ~ elev',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-0.5, 1.5))
dev.off()


# # caterpillar plots - SI ~ temp time --------------------------------------------------

si_temp_l0_sc_gamma <- MCMCvis::MCMCchains(si_temp_l0_fit, params = 'gamma') / si_temp_l0_data$scf_temp
si_temp_l1_sc_gamma <- MCMCvis::MCMCchains(si_temp_l1_fit, params = 'gamma') / si_temp_l1_data$scf_temp
si_temp_l2_sc_gamma <- MCMCvis::MCMCchains(si_temp_l2_fit, params = 'gamma') / si_temp_l2_data$scf_temp

pdf(paste(fig_dir, 'gamma_si_temp_l0_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_temp_l0_sc_gamma,
                  params = 'gamma',
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'gamma - SI ~ temp (lag 0)',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-0.15, 0.05))
dev.off()

pdf(paste(fig_dir, 'gamma_si_temp_l1_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_temp_l1_sc_gamma,
                  params = 'gamma',
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'gamma - SI ~ temp (lag 1)',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-0.15, .1))
dev.off()

pdf(paste(fig_dir, 'gamma_si_temp_l2_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_temp_l2_sc_gamma,
                  params = 'gamma',
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'gamma - SI ~ temp (lag 2)',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-0.1, 0.05))
dev.off()

# caterpillar plots - SI ~ temp space --------------------------------------------------

pdf(paste(fig_dir, 'T_effect_mn_temp_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(si_temp_space_fit,
                  params = 'beta',
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'T effect ~ mean temp',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-3, 2))
dev.off()


# caterpillar plots - absolute TLE --------------------------------------------------

pdf(paste(fig_dir, 'abs_mass_year_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(per_ch_year[,1,],
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'absolute mass ~ year',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-10, 20))
dev.off()

pdf(paste(fig_dir, 'abs_mass_lat_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(per_ch_lat[,1,],
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'absolute mass ~ lat',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-40, 80))
dev.off()

pdf(paste(fig_dir, 'abs_mass_elev_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(per_ch_elev[,1,],
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'absolute mass ~ elev',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-30, 30))
dev.off()

pdf(paste(fig_dir, 'abs_wing_year_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(per_ch_year[,2,],
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'absolute wing ~ year',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-4, 4))
dev.off()

pdf(paste(fig_dir, 'abs_wing_lat_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(per_ch_lat[,2,],
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'absolute wing ~ lat',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-10, 20))
dev.off()

pdf(paste(fig_dir, 'abs_wing_elev_full_cat.pdf'), height = 16, width = 8)
MCMCvis::MCMCplot(per_ch_elev[,2,],
                  labels = usp,
                  guide_lines = TRUE,
                  ci = c(50, 89),
                  sz_labels = 0.75,
                  main = 'absolute wing ~ elev',
                  sz_med = 1.75,
                  sz_thick = 5,
                  sz_thin = 3,
                  xlim = c(-10, 15))
dev.off()

