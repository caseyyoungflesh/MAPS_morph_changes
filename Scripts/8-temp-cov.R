####################
# temp effect space ~ mn_temp + phylo
####################


# set dir ----------------------------------------------------------------

run_date <- '2022-06-01'
maps_master_date <- '2021-07-28'
daymet_process_date <- '2021-04-02'
temp_space_date <- '2022-05-25'
phylo_date <- '2021-05-19'

dir <- 'XXXX'


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)
library(MCMCvis)
library(ape)
library(phytools)


# process input data space -----------------------------------------------------

#read in temp space data
setwd(paste0(dir, '/Results/si-temp-space-', temp_space_date))
temp_space_data <- readRDS(paste0('si-temp-space-data-', temp_space_date, '.rds'))
temp_space_fit <- readRDS(paste0('si-temp-space-fit-', temp_space_date, '.rds'))

#extract posteriors for temp effect
beta_mn_space <- MCMCvis::MCMCpstr(temp_space_fit, params = 'beta')[[1]]
beta_sd_space <- MCMCvis::MCMCpstr(temp_space_fit, params = 'beta', fun = sd)[[1]]

#read in daymet data
setwd(paste0(dir, '/Data/daymet/processed/'))
daymet <- readRDS(paste0('daymet-master-', daymet_process_date, '.rds')) %>%
  mutate(MJJ_tmax = (May_tmax + June_tmax + July_tmax) / 3) #using mean here does not do what might be expected

#merge morph and temp data
mdata <- temp_space_data$pro_data
mdata2 <- dplyr::select(mdata, sci_name, sp_id, station, cn_id, lat, 
                        lng, year, GMTED_elev, size_idx) %>%
  dplyr::left_join(daymet, by = c('station', 'year')) %>%
  dplyr::arrange(sp_id, cn_id)

#scale factor for temp
scf_temp <- 10

#mean temp for each species, centered
mn_st_temp <- mdata2 %>%
  dplyr::group_by(sp_id, cn_id) %>%
  dplyr::summarize(mn_MJJ_tmax = mean(MJJ_tmax)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sp_id) %>%
  dplyr::summarize(smt = mean(mn_MJJ_tmax)) %>%
  dplyr::mutate(sp_mean_temp = scale(smt, scale = FALSE)[,1] / scf_temp) %>%
  dplyr::ungroup()

#read in phylo species names
setwd(paste0(dir, '/Data/bird_phylo'))
bn <- read.table(paste0('phylo_names-', phylo_date, '.txt'), sep = ',')[,1]
bn2 <- gsub(' ', '_', bn)
bn3 <- data.frame(name = bn2, num = 1:length(bn2))

setwd(paste0(dir, "/Data/bird_phylo/"))
pnk <- read.csv(paste0('phylo_names_key-', phylo_date, '.csv'))

#read in phylo tree
setwd(paste0(dir, "/Data/bird_phylo/tree-pruner-d8509cd1-17ae-46f6-ae43-dc270e4deff2"))
ptree <- ape::read.nexus('output.nex')

#consensus tree
tree <- phytools::consensus.edges(ptree)

#calculate covariance matrix
V <- ape::vcv.phylo(tree)

#scale by max variance -> correlation matrix
R <- V[bn2, bn2] / max(V)


# Data space --------------------------------------------------------------

#number of data points
N <- length(beta_mn_space)

DATA_sp <- list(N = N,
                y = beta_mn_space,
                sd_y = beta_sd_space,
                R = R,
                I = diag(1, nrow = N, ncol = N),
                mn_temp = mn_st_temp$sp_mean_temp,
                mn_st_temp = mn_st_temp,
                scf_temp = scf_temp,
                pro_data = mdata2)


# Call Stan model space --------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.93
TREE_DEPTH <- 12
STEP_SIZE <- 0.1
CHAINS <- 4
ITER <- 5000

fit_sp <- rstan::stan(paste0(dir, 'Scripts/Model_files/temp-sp-cov.stan'),
                      data = DATA_sp,
                      chains = CHAINS,
                      iter = ITER,
                      cores = CHAINS,
                      pars = c('alpha',
                               'beta',
                               'lambda',
                               'sigma',
                               'mu_y',
                               'y_rep'), 
                      control = list(adapt_delta = DELTA,
                                     max_treedepth = TREE_DEPTH,
                                     stepsize = STEP_SIZE))


# save summary space ------------------------------------------------------------

setwd(paste0(dir, '/Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit_sp, 
                  round = 4,
                  file_name = paste0('temp-sp-cov-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('temp-sp-cov-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('temp-sp-cov-fit-', run_date),
                  add_obj = list(DATA_sp),
                  add_obj_names = paste0('temp-sp-cov-data-', run_date),
                  cp_file = c('Model_files/temp-sp-cov.stan', 
                              '8-temp-cov.R'),
                  cp_file_names = c(paste0('temp-sp-cov-', run_date, '.stan'),
                                    paste0('8-temp-cov-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit_l0)


# PPC ---------------------------------------------------------------------

y_val <- DATA_sp$y
y_rep <- MCMCvis::MCMCchains(fit_sp, params = 'y_rep')

pdf(paste0(dir, 'Results/temp-sp-cov-', run_date, '/PPC.pdf'), height = 5, width = 5)
plot(density(y_val), lwd = 2, main = 'PPC', xlab = 'Value', ylim = c(0, 1.1))
for (i in 1:100)
{
  lines(density(y_rep[i,]), col = rgb(1,0,0,0.1))
}
dev.off()


# PPO ---------------------------------------------------------------------

#alpha ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit_sp,
                   params = 'alpha',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('temp-sp-cov-trace-alpha-', run_date, '.pdf'))

#beta ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit_sp,
                   params = 'beta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('temp-sp-cov-trace-beta-', run_date, '.pdf'))

#sigma ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_sp,
                   params = 'sigma',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('temp-sp-cov-trace-sigma-', run_date, '.pdf'))

#lambda ~ U(0, 1)
PR <- runif(10000, 0, 1)
MCMCvis::MCMCtrace(fit_sp,
                   params = 'lambda',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('temp-sp-cov-trace-lambda-', run_date, '.pdf'))
