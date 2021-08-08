####################
# Morph change rate
#
####################


# set dir ----------------------------------------------------------------

dir <- 'XXXX'
fig_dir <- paste0(dir, '/Figures/paper_figs/')
si_tle_date <- '2021-07-28'
per_ch_date <- '2021-08-02'


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)
library(MCMCvis)


# input data ------------------------------------------------------

#read in morph fit data
setwd(paste0(dir, '/Results/si-tle-', si_tle_date))
morph_data <- readRDS(paste0('si-tle-data-', si_tle_date, '.rds'))

#Simulate trait value at start and end of period, extract gen time from Bird et al. data, calculate haldanes and compare to Lynch and Burger and Hendry and Kinnison 1999.

#From Bird et al. 2020 Con Bio
gen_time <- read.csv(paste0(dir, 'Data/Gen_time_Bird_et_al_2020_table_4.csv'))

#rate of change
rdf <- readRDS(paste0(dir, 'Results/per-ch-table-data-', per_ch_date, '.rds'))

# Hendy et al. 2008 Mol Eco
hea_2008 <- read.csv(paste0(dir, 'Data/Hendry_et_al_2008_Mol_Eco.csv'))


# process data ------------------------------------------------------------

#station ll
morph_data$pro_data %>%
  dplyr::distinct(station, lat) -> st_ll

#calculate s_p for each trait for each species (mean of station-specific standard deviations)
morph_data$pro_data %>%
  dplyr::group_by(sci_name, station) %>%
  dplyr::summarize(sd_l_wing = sd(log(wing_chord)), sd_l_mass = sd(log(weight))) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(st_ll, by = 'station') -> met_per_sp_st

met_per_sp_st %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(mn_sd_l_wing = mean(sd_l_wing, na.rm = TRUE), 
                   mn_sd_l_mass = mean(sd_l_mass, na.rm = TRUE)) -> mn_sd_sp

#rate of change
dplyr::mutate(rdf, ch_mass = (l_mass_year * 15), ch_wing = (l_wing_year * 15)) %>%
  dplyr::select(sci_name, ch_mass, ch_wing) -> ch_met

morph_data$pro_data %>%
  dplyr::group_by(sci_name) %>%
  dplyr::summarize(mean_l_mass = mean(log(weight)), mean_l_wing = mean(log(wing_chord))) %>%
  dplyr::left_join(ch_met, by = 'sci_name') %>%
  dplyr::mutate(x1_l_mass = mean_l_mass - ch_mass, x2_l_mass = mean_l_mass + ch_mass,
                x1_l_wing = mean_l_wing - ch_wing, x2_l_wing = mean_l_wing + ch_wing) %>%
  dplyr::left_join(mn_sd_sp) -> pro_df

#gen length (GenLength)
gen_time$sci_name <- gen_time$Scientific.name

# #which names don't match
# pro_df %>%
#   dplyr::left_join(tibble(gen_time), by = c('sci_name' = 'Scientific.name')) -> tt
# 
# tt$sci_name[which(is.na(tt$Sequence))]
#Dryobates villosus -> Leuconotopicus villosus
#Icterus bullockii -> Icterus bullockiorum
#Oreothlypis celata -> Leiothlypis celata
#Oreothlypis luciae -> Leiothlypis luciae
#Oreothlypis peregrina -> Leiothlypis peregrina
#Oreothlypis ruficapilla -> Leiothlypis ruficapilla
#Oreothlypis virginiae -> Leiothlypis virginiae

study_names <- c('Dryobates villosus',
                 'Icterus bullockii',
                 'Oreothlypis celata',
                 'Oreothlypis luciae',
                 'Oreothlypis peregrina',
                 'Oreothlypis ruficapilla',
                 'Oreothlypis virginiae')
bird_ea_names <- c('Leuconotopicus villosus',
                   'Icterus bullockiorum',
                   'Leiothlypis celata',
                   'Leiothlypis luciae',
                   'Leiothlypis peregrina',
                   'Leiothlypis ruficapilla',
                   'Leiothlypis virginiae')

gen_time$sci_name[which(gen_time$sci_name %in% bird_ea_names)] <- study_names

pro_df %>%
  dplyr::left_join(gen_time, by = 'sci_name') -> pro_df2


# calculate haldanes -------------------------------------------------------

#from Hendry and Kinnison 1999 (developed by Gingerich 1983)
#h = ((x_2 / s_p) - (x_1 / s_p)) / g
#h = haldane
#x_2 = mean value for trait at end time point
#x_1 = mean value for trait at start time point
#s_p = standard deviation of trait (pooled across time)
#g = number of generations between two time points

pro_df2 %>%
  dplyr::mutate(haldane_l_mass = ((x2_l_mass / mn_sd_l_mass) - (x1_l_mass / mn_sd_l_mass)) / (30 / GenLength),
                haldane_l_wing = ((x2_l_wing / mn_sd_l_wing) - (x1_l_wing / mn_sd_l_wing)) / (30 / GenLength)) -> pro_df3

#in situ antho disturbance only for Hendry et al. 2008 data
hea_2008 %>%
  dplyr::filter(Disturbance == 2, !is.na(Haldanes..absolute.value.)) -> hea_2008_f1

#save figure 
pdf(paste0(fig_dir, 'haldanes.pdf'), height = 4, width = 4)
hist(hea_2008_f1$Haldanes..absolute.value., lwd = 2, xlim = c(0, 0.2), 
     main = '', breaks = 60, ylim = c(0, 50), col = rgb(0,0,1,0.5), xlab = 'Rate of phenotypic change (haldanes)')
hist((abs(pro_df3$haldane_l_mass)), col = rgb(1,0,0,0.5), lwd = 2, add = TRUE, breaks = 10)
dev.off()


#number species haldane > 0.1 (considered max rate of sustained change)
sum(pro_df3$haldane_l_mass >= 0.1)

