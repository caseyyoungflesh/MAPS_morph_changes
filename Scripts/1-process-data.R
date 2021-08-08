####################
# Process MAPS data
#
####################


# set dirs ----------------------------------------------------------------

RUN_DATE <- '2021-07-28'
PHYLO_DATE <- '2021-05-19'

dir <- 'XXXX'


# load packages -----------------------------------------------------------

library(tidyverse)
library(ape)
library(caper)


# input data -------------------------------------------------------------

#read in data
#only AHY (age codes: 5, 1, 6, 7, or 8) males at stations where species was breeding in that year, between day 121 and 220, and - or + standard effort codes. 
#all entries with weight or wing chord NA or 0 excluded as were stations with lat 0
setwd(paste0(dir, '/Data'))
mdata3 <- read.csv('MAPS_data.csv')
usp3 <- sort(unique(mdata3$sci_name))


# remove anomalous values --------------------------------------------------------------

#create empty df
mdata4 <- mdata3
mdata4[,] <- NA 
counter <- 1
#fill df
for (i in 1:length(usp3))
{
  #i <- 1
  print(paste0('species ', i, ' of ', length(usp3)))
  temp <- dplyr::filter(mdata3, sci_name == usp3[i])
  
  #remove all values more than 5 MAD away (v conservative) from median - Leys et al. 2013
  med_wc <- median(temp$wing_chord)
  med_weight <- median(temp$weight)
  
  MAD_wc <- mad(temp$wing_chord)
  MAD_weight <- mad(temp$weight)
  
  low_wc <- med_wc - 5 * MAD_wc
  high_wc <- med_wc + 5 * MAD_wc
  low_weight <- med_weight - 5 * MAD_weight
  high_weight <- med_weight + 5 * MAD_weight
  
  to.rm <- which(temp$wing_chord < low_wc | temp$wing_chord > high_wc |
                   temp$weight < low_weight | temp$weight > high_weight)
  if (length(to.rm) > 0)
  {
    td <- temp[-to.rm,]
    mdata4[counter:(counter + NROW(td) - 1),] <- td
    counter <- counter + NROW(td)
  } else {
    mdata4[counter:(counter + NROW(temp) - 1),] <- temp
    counter <- counter + NROW(temp)
  }
}

#trim excess 0s from df
#only first capture of an individual in a season
mdata5 <- mdata4 %>%
  dplyr::filter(!is.na(sci_name)) %>%
  dplyr::group_by(cp_id) %>%
  dplyr::slice_min(cap_order) %>%
  dplyr::ungroup() %>%
  dplyr::select(-cap_order, -cp_id)
  

# filter by species -------------------------------------------------------

#only species with >= 375 data points
d_cnt <- plyr::count(mdata5, 'sci_name')
sp <- dplyr::filter(d_cnt, freq >= 375)[,1]

#arrange and add sp_id
mdata6 <- dplyr::filter(mdata5, sci_name %in% sp) %>%
  dplyr::arrange(sci_name, station, year) %>%
  dplyr::mutate(sp_id = as.numeric(factor(sci_name)))


# estimate scaling exponent ----------------------------------------------------

# calculate species-level log(mass) and log(wing length)
mn_morph <- dplyr::group_by(mdata6, sci_name) %>%
  dplyr::summarize(mean_l_weight = mean(log(weight)),
                   mean_l_wing_chord = mean(log(wing_chord))) %>%
  dplyr::arrange(desc(mean_l_weight))

#relationship between wing length and body size follows a power law; wing length = mass^s. 
#Log both wing length and mass to linearize the relationship and estimate slope
#data from birdtree.org
setwd(paste0(dir, "/Data/bird_phylo/tree-pruner-d8509cd1-17ae-46f6-ae43-dc270e4deff2"))
ptree <- ape::read.nexus('output.nex')

setwd(paste0(dir, "/Data/bird_phylo/"))
pnk <- read.csv(paste0('phylo_names_key-', PHYLO_DATE, '.csv')) %>%
  dplyr::distinct(phylo_sci_name, .keep_all = TRUE)

mn_morph2 <- dplyr::left_join(mn_morph, pnk, by = 'sci_name') %>%
  dplyr::arrange(sci_name) %>%
  as.data.frame()
mn_morph2$phylo_sci_name <- gsub(' ', '_', mn_morph2$phylo_sci_name)

# run phylogenetic regression for each tree realization
phy_reg_int <- rep(NA, length(ptree))
phy_reg_sl <- rep(NA, length(ptree))
for (i in 1:length(ptree))
{
  #i <- 22
  tree_n <- ptree[[i]]
  
  td <- caper::comparative.data(tree_n, mn_morph2, phylo_sci_name, vcv = TRUE, vcv.dim = 3)
  tf <- caper::pgls(mean_l_wing_chord ~ mean_l_weight, td)
  phy_reg_int[i] <- as.numeric(tf$model$coef[1])
  phy_reg_sl[i] <- as.numeric(tf$model$coef[2])
}

#mean intercept of regression from all trees
int <- mean(phy_reg_int)
#mean slope of regressions from all trees
sl <- mean(phy_reg_sl)
sd_sl <- sd(phy_reg_sl)


# reproject data points ---------------------------------------------------

# Rotate so that x-axis is the size of the bird (size index), and y axis is the wingy-ness of the bird (wingy-ness index).

#rotate points by slope of log(wl) log(mass) relationship
#get angle to rotate points (in radians)
angle <- -atan(sl)

#create rotation matrix
#R = \begin{bmatrix} cos(\theta) & -sin(\theta) \\ sin(\theta) & cos(\theta) \end{bmatrix}
#\begin{bmatrix} x' \\ y' \end{bmatrix} = \begin{bmatrix} cos(\theta) & -sin(\theta) \\ sin(\theta) & cos(\theta) \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix}
#rotm = [[cos(theta) -sin(theta]
#        [sin(theta)  cos(theta)]]
#90 degrees
# rotm <- matrix(c(0, -1, 
#                  1, 0), ncol = 2)
rotm <- matrix(c(cos(angle), sin(angle), 
                 -sin(angle), cos(angle)), ncol = 2)

#function to log, rotate, then scale data
# std = TRUE -> standardized
# std = FALSE -> percent deviation from mean
# percent deviation suffers from not being so interpretatble bc morph metrics were logged - could use though
rotate_fun <- function(input)
{
  dm <- cbind(log(input$weight), log(input$wing_chord))
  M <- t(rotm %*% t(dm))
  
  mn_x <- mean(M[,1])
  sd_x <- sd(M[,1])
  mn_y <- mean(M[,2])
  sd_y <- sd(M[,2])
  
  M_sc <- apply(M, 2, function(x) scale(x, scale = TRUE))
  
  #mean and sd for each rotated axis
  ol <- list(M_sc = M_sc,
             mn_x = mn_x,
             sd_x = sd_x,
             mn_y = mn_y,
             sd_y = sd_y)
  
  return(ol)
}

#apply rotation matrix to data for each species individually and scale (center and std) to get 'size index' (x-axis of rotated data) and 'wing index' (y-axis of rotated data)
usp4 <- sort(unique(mdata6$sci_name))
mdata7 <- mdata6
mdata7$size_idx <- NA
mdata7$wingy_idx <- NA
mn_morph2$mn_x <- NA
mn_morph2$sd_x <- NA
mn_morph2$mn_y <- NA
mn_morph2$sd_y <- NA
for (i in 1:length(usp4))
{
  #candidates: 37
  #i <- 98
  idx <- which(mdata7$sci_name == usp4[i])
  temp <- mdata7[idx, ]
  
  #rotate
  rpd <- rotate_fun(temp)
  
  mdata7$size_idx[idx] <- rpd$M_sc[,1]
  mdata7$wingy_idx[idx] <- rpd$M_sc[,2]
  
  mn_idx <- which(mn_morph2$sci_name == usp4[i])
  mn_morph2$mn_x[mn_idx] <- rpd$mn_x
  mn_morph2$sd_x[mn_idx] <- rpd$sd_x
  mn_morph2$mn_y[mn_idx] <- rpd$mn_y
  mn_morph2$sd_y[mn_idx] <- rpd$sd_y
}

# apply rotation matrix to species-level means
M_sp <- t(rotm %*% t(cbind(mn_morph2$mean_l_weight, mn_morph2$mean_l_wing_chord)))
M_sp_sc <- apply(M_sp, 2, function(x) scale(x, scale = TRUE))
mn_morph2$sp_size_idx <- M_sp_sc[,1]
mn_morph2$sp_wingy_idx <- M_sp_sc[,2]


# locations for daymet extraction -----------------------------------------

#write locations to file
#arrange according to: https://github.com/ornldaac/daymet-single-pixel-batch
setwd(paste0(dir, "Scripts/3-process-daymet"))
#range(unique(mdata7$year))
write.table(unique(data.frame(lat = mdata7$lat, lon = mdata7$lng)),
            'daymet_query.txt', sep = ',', row.names = FALSE, col.names = FALSE)
#add variables to top of txt file
system(paste0("echo 'Variables: tmin, tmax, prcp' | cat - daymet_query.txt > temp && mv temp daymet_query.txt"))


#DATA STATS
# num samples - 253488
NROW(mdata7)
# num stations - 1124
length(unique(mdata7$station))
# num species - 105
length(unique(mdata7$sci_name))
# range years - 1989-2018
range(mdata7$year)
# range lat - 26.1 - 69.4
range(mdata7$lat)

#merge with species-level data
mdata7 %>%
  dplyr::left_join(mn_morph2, by = c('sci_name')) -> mdata8


# save rds ----------------------------------------------------------------

setwd(paste0(dir, '/Data'))
saveRDS(mdata8, paste0('MAPS-processed-', RUN_DATE, '.rds'))

