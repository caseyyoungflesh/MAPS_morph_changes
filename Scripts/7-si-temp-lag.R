####################
# SI ~ temp on lag 0, lag 1, and lag 2
####################


# set dir ----------------------------------------------------------------

run_date <- '2022-05-25'
maps_master_date <- '2021-07-28'
idx_tle_date <- '2021-07-28'
daymet_process_date <- '2021-04-02'

dir <- 'XXXX'


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)
library(MCMCvis)


# process input data -------------------------------------------------------------

#read in daymet data
setwd(paste0(dir, '/Data/daymet/processed/'))
daymet <- readRDS(paste0('daymet-master-', daymet_process_date, '.rds')) %>%
  mutate(MJJ_tmax = (May_tmax + June_tmax + July_tmax) / 3) #using mean here does not do what might be expected

#read in morph fit data
setwd(paste0(dir, '/Results/si-tle-', idx_tle_date))
morph_data <- readRDS(paste0('si-tle-data-', idx_tle_date, '.rds'))


# lag temp data ----------------------------------------------------------------

#data from 4-si-tle
mdata <- morph_data$pro_data

#scale factor
#center temp and scale so param estimates are near unit scale for sampling efficiency
scf_temp <- 10

#lag temp data by 0 to 2 years
YR_LAG <- 0
daymet_l0 <- daymet
daymet_l0$actual_year <- daymet_l0$year
daymet_l0$year <- daymet_l0$year + YR_LAG
mdata_l0 <- dplyr::left_join(mdata, daymet_l0, by = c('station', 'year')) %>%
  dplyr::mutate(sc_temp = scale(MJJ_tmax, scale = FALSE)[,1] / scf_temp)

YR_LAG <- 1
daymet_l1 <- daymet
daymet_l1$actual_year <- daymet_l1$year
daymet_l1$year <- daymet_l1$year + YR_LAG
mdata_l1 <- dplyr::left_join(mdata, daymet_l1, by = c('station', 'year')) %>%
  dplyr::mutate(sc_temp = scale(MJJ_tmax, scale = FALSE)[,1] / scf_temp)

YR_LAG <- 2
daymet_l2 <- daymet
daymet_l2$actual_year <- daymet_l2$year
daymet_l2$year <- daymet_l2$year + YR_LAG
mdata_l2 <- dplyr::left_join(mdata, daymet_l2, by = c('station', 'year')) %>%
  dplyr::mutate(sc_temp = scale(MJJ_tmax, scale = FALSE)[,1] / scf_temp)

#sp_id for each cn_id
cn_sp <- unique(mdata_l0[,c('sp_id', 'cn_id')])


# Data L0 --------------------------------------------------------------

#calculate mean temp each cn_id, then scale within each species
mn_st_temp_df_l0 <- mdata_l0 %>%
  dplyr::group_by(sp_id, cn_id) %>%
  dplyr::summarize(mn_st_temp = mean(MJJ_tmax)) %>%
  dplyr::mutate(sc_mn_st_temp = scale(mn_st_temp, scale = FALSE)[,1] / scf_temp) %>%
  dplyr::ungroup()

DATA_l0 <- list(N = NROW(mdata_l0), #number of temp data points
             Nsp = length(unique(mdata_l0$sp_id)), #number of species
             Nsc = length(unique(mdata_l0$cn_id)), #number of y data points
             y = mdata_l0$size_idx,
             temp = mdata_l0$sc_temp,
             mn_st_temp = mn_st_temp_df_l0$sc_mn_st_temp, #mean temp at each station, centered for each species
             sp = mdata_l0$sp_id,
             cn_id = mdata_l0$cn_id, 
             cn_sp = cn_sp$sp_id, 
             daymet = daymet,
             scf_temp = scf_temp,
             mn_st_temp_df_l0 = mn_st_temp_df_l0,
             pro_data = mdata_l0) 


# Call Stan model L0 --------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.93
TREE_DEPTH <- 12
STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 6000

fit_l0 <- rstan::stan(paste0(dir, 'Scripts/Model_files/morph-temp-ss2.stan'),
                   data = DATA_l0,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma',
                            'theta',
                            'gamma2',
                            'theta2',
                            'mu_gamma',
                            'mu_theta',
                            'mu_gamma2',
                            'mu_theta2',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma_gt',
                            'sigma_gt2',
                            'Rho_gt',
                            'Rho_gt2',
                            'sigma_proc',
                            'lambda_proc',
                            'kappa_proc',
                            'nu'),
                   #'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# PPC sim --------------------------------------------------------

#returns list with y, cn_id, sp, y_rep, and mu
PPC_fun <- function(obj,
                    DD,
                    Nsim = 100, 
                    seed = 1)
{
  print('Extracting posteriors')
  alpha_ch <- MCMCvis::MCMCchains(obj, params = 'alpha')
  beta_ch <- MCMCvis::MCMCchains(obj, params = 'beta')
  nu_ch <- MCMCvis::MCMCchains(obj, params = 'nu')
  sigma_proc_ch <- MCMCvis::MCMCchains(obj, params = 'sigma_proc')
  
  #random iterations to draw from posteriors
  idx <- sample(1:NROW(alpha_ch), Nsim)
  
  print('Simulating data')
  y_rep <- matrix(NA, nrow = Nsim, ncol = DD$N)
  mu <- matrix(NA, nrow = Nsim, ncol = DD$N)
  
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = (DD$N * Nsim), style = 3)
  counter <- 0
  for (i in 1:length(idx))
  {
    for (j in 1:DD$N)
    {
      #i <- 1
      #j <- 1
      mu[i,j] <- alpha_ch[idx[i], DD$cn_id[j]] + 
        beta_ch[idx[i], DD$cn_id[j]] * DD$temp[j]
      
      eps <- rt(n = 1, df = nu_ch[idx[i]]) * sigma_proc_ch[idx[i], DD$sp[j]]
      y_rep[i,j] <- mu[i,j] + eps
      
      #progress bar
      setTxtProgressBar(pb, counter)
      counter <- counter + 1
    }
  }
  
  return(list(y = DD$y, 
              sci_name = DD$pro_data$sci_name, 
              cn_id = DD$cn_id, 
              y_rep = y_rep, 
              mu = mu))
}

#run function
ppc_sim_l0 <- PPC_fun(obj = fit_l0,
                      DD = DATA_l0,
                      Nsim = 100)

# save summary L0 ------------------------------------------------------------

setwd(paste0(dir, '/Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit_l0, 
                  round = 4,
                  file_name = paste0('si-temp-l0-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('si-temp-l0-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('si-temp-l0-fit-', run_date),
                  add_obj = list(DATA_l0, ppc_sim_l0),
                  add_obj_names = c(paste0('si-temp-l0-data-', run_date),
                                    paste0('ppc-sim-l0-', run_date)),
                  cp_file = c('Model_files/morph-temp-ss2.stan', 
                              '7-si-temp-lag.R'),
                  cp_file_names = c(paste0('morph-temp-ss2-', run_date, '.stan'),
                                    paste0('7-si-temp-lag-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)


# Data L1 --------------------------------------------------------------

#calculate mean temp each cn_id, then scale within each species
mn_st_temp_df_l1 <- mdata_l1 %>%
  dplyr::group_by(sp_id, cn_id) %>%
  dplyr::summarize(mn_st_temp = mean(MJJ_tmax)) %>%
  dplyr::mutate(sc_mn_st_temp = scale(mn_st_temp, scale = FALSE)[,1] / scf_temp) %>%
  dplyr::ungroup()

DATA_l1 <- list(N = NROW(mdata_l1), #number of temp data points
                Nsp = length(unique(mdata_l1$sp_id)), #number of species
                Nsc = length(unique(mdata_l1$cn_id)), #number of y data points
                y = mdata_l1$size_idx,
                temp = mdata_l1$sc_temp,
                mn_st_temp = mn_st_temp_df_l1$sc_mn_st_temp, #mean temp at each station, centered for each species
                sp = mdata_l1$sp_id,
                cn_id = mdata_l1$cn_id,
                cn_sp = cn_sp$sp_id,
                daymet = daymet,
                scf_temp = scf_temp,
                mn_st_temp_df_l1 = mn_st_temp_df_l1,
                pro_data = mdata_l1)


# Call Stan model l1 --------------------------------------------------------------

fit_l1 <- rstan::stan(paste0(dir, 'Scripts/Model_files/morph-temp-ss2.stan'),
                   data = DATA_l1,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma',
                            'theta',
                            'gamma2',
                            'theta2',
                            'mu_gamma',
                            'mu_theta',
                            'mu_gamma2',
                            'mu_theta2',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma_gt',
                            'sigma_gt2',
                            'Rho_gt',
                            'Rho_gt2',
                            'sigma_proc',
                            'lambda_proc',
                            'kappa_proc',
                            'nu'),
                   #'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# PPC sim -----------------------------------------------------------------

ppc_sim_l1 <- PPC_fun(obj = fit_l1,
                      DD = DATA_l1,
                      Nsim = 100)


# save summary l1 ------------------------------------------------------------

setwd(paste0(dir, '/Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit_l1, 
                  round = 4,
                  file_name = paste0('si-temp-l1-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('si-temp-l1-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('si-temp-l1-fit-', run_date),
                  add_obj = list(DATA_l1, ppc_sim_l1),
                  add_obj_names = c(paste0('si-temp-l1-data-', run_date),
                                    paste0('ppc-sim-l1-', run_date)),
                  cp_file = c('Model_files/morph-temp-ss2.stan', 
                              '7-si-temp-lag.R'),
                  cp_file_names = c(paste0('morph-temp-ss2-', run_date, '.stan'),
                                    paste0('7-si-temp-lag-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)


# Data l2 --------------------------------------------------------------

#calculate mean temp each cn_id, then scale within each species
mn_st_temp_df_l2 <- mdata_l2 %>%
  dplyr::group_by(sp_id, cn_id) %>%
  dplyr::summarize(mn_st_temp = mean(MJJ_tmax)) %>%
  dplyr::mutate(sc_mn_st_temp = scale(mn_st_temp, scale = FALSE)[,1] / scf_temp) %>%
  dplyr::ungroup()

DATA_l2 <- list(N = NROW(mdata_l2), #number of temp data points
                Nsp = length(unique(mdata_l2$sp_id)), #number of species
                Nsc = length(unique(mdata_l2$cn_id)), #number of y data points
                y = mdata_l2$size_idx,
                temp = mdata_l2$sc_temp,
                mn_st_temp = mn_st_temp_df_l2$sc_mn_st_temp, #mean temp at each station, centered for each species
                sp = mdata_l2$sp_id,
                cn_id = mdata_l2$cn_id, 
                cn_sp = cn_sp$sp_id, 
                daymet = daymet,
                scf_temp = scf_temp,
                mn_st_temp_df_l2 = mn_st_temp_df_l2,
                pro_data = mdata_l2)


# Call Stan model l2 --------------------------------------------------------------

fit_l2 <- rstan::stan(paste0(dir, 'Scripts/Model_files/morph-temp-ss2.stan'),
                   data = DATA_l2,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma',
                            'theta',
                            'gamma2',
                            'theta2',
                            'mu_gamma',
                            'mu_theta',
                            'mu_gamma2',
                            'mu_theta2',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma_gt',
                            'sigma_gt2',
                            'Rho_gt',
                            'Rho_gt2',
                            'sigma_proc',
                            'lambda_proc',
                            'kappa_proc',
                            'nu'),
                   #'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# PPC sim -----------------------------------------------------------------

ppc_sim_l2 <- PPC_fun(obj = fit_l2,
                      DD = DATA_l2,
                      Nsim = 100)


# save summary l2 ------------------------------------------------------------

setwd(paste0(dir, '/Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit_l2, 
                  round = 4,
                  file_name = paste0('si-temp-l2-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('si-temp-l2-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('si-temp-l2-fit-', run_date),
                  add_obj = list(DATA_l2, ppc_sim_l2),
                  add_obj_names = c(paste0('si-temp-l2-data-', run_date),
                                    paste0('ppc-sim-l2-', run_date)),
                  cp_file = c('Model_files/morph-temp-ss2.stan', 
                              '7-si-temp-lag.R'),
                  cp_file_names = c(paste0('morph-temp-ss2-', run_date, '.stan'),
                                    paste0('7-si-temp-lag-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)


# PPC plots ----------------------------------------------------------------

ppc_plot_fun <- function(sim_obj, file_name)
{
  #overall
  pdf(paste0(file_name), height = 5, width = 5)
  plot(density(sim_obj$y), lwd = 2, main = 'PPC', xlab = 'Value')
  for (i in 1:NROW(sim_obj$y_rep))
  {
    lines(density(sim_obj$y_rep[i,]), col = rgb(1,0,0,0.1))
  }
  dev.off()
}

setwd(paste0(dir, '/Results/si-temp-l0-', run_date))
ppc_plot_fun(sim_obj = ppc_sim_l0, file_name = 'PPC_l0.pdf')

setwd(paste0(dir, '/Results/si-temp-l1-', run_date))
ppc_plot_fun(sim_obj = ppc_sim_l1, file_name = 'PPC_l1.pdf')

setwd(paste0(dir, '/Results/si-temp-l2-', run_date))
ppc_plot_fun(sim_obj = ppc_sim_l2, file_name = 'PPC_l2.pdf')


# PPO ---------------------------------------------------------------------

#lag 0
#mu_gamma ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l0,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-mu_gamma-', run_date, '.pdf'))

#mu_theta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l0,
                   params = 'mu_theta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-mu_theta-', run_date, '.pdf'))

#mu_gamma2 ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l0,
                   params = 'mu_gamma2',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-mu_gamma2-', run_date, '.pdf'))

#mu_theta2 ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l0,
                   params = 'mu_theta2',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-mu_theta2-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l0,
                   params = 'sigma_alpha',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-sigma_alpha-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l0,
                   params = 'sigma_beta',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-sigma_beta-', run_date, '.pdf'))

#lambda_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l0,
                   params = 'lambda_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-lambda_proc-', run_date, '.pdf'))

#kappa_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l0,
                   params = 'kappa_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-kappa_proc-', run_date, '.pdf'))

#sigma_gt ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l0,
                   params = 'sigma_gt',
                   ISB = 'FALSE',
                   exact = FALSE,
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-sigma_gt-', run_date, '.pdf'))

#sigma_gt2 ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l0,
                   params = 'sigma_gt2',
                   ISB = 'FALSE',
                   exact = FALSE,
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-sigma_gt2-', run_date, '.pdf'))

#nu ~ gamma(2, 0.1)
PR <- rgamma(10000, 2, 0.1)
MCMCvis::MCMCtrace(fit_l0,
                   params = 'nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l0-ss2-trace-nu-', run_date, '.pdf'))

#lag 1
#mu_gamma ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l1,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-mu_gamma-', run_date, '.pdf'))

#mu_theta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l1,
                   params = 'mu_theta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-mu_theta-', run_date, '.pdf'))

#mu_gamma2 ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l1,
                   params = 'mu_gamma2',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-mu_gamma2-', run_date, '.pdf'))

#mu_theta2 ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l1,
                   params = 'mu_theta2',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-mu_theta2-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l1,
                   params = 'sigma_alpha',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-sigma_alpha-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l1,
                   params = 'sigma_beta',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-sigma_beta-', run_date, '.pdf'))

#lambda_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l1,
                   params = 'lambda_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-lambda_proc-', run_date, '.pdf'))

#kappa_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l1,
                   params = 'kappa_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-kappa_proc-', run_date, '.pdf'))

#sigma_gt ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l1,
                   params = 'sigma_gt',
                   ISB = 'FALSE',
                   exact = FALSE,
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-sigma_gt-', run_date, '.pdf'))

#sigma_gt2 ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l1,
                   params = 'sigma_gt2',
                   ISB = 'FALSE',
                   exact = FALSE,
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-sigma_gt2-', run_date, '.pdf'))

#nu ~ gamma(2, 0.1)
PR <- rgamma(10000, 2, 0.1)
MCMCvis::MCMCtrace(fit_l1,
                   params = 'nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l1-ss2-trace-nu-', run_date, '.pdf'))

#lag 2
#mu_gamma ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l2,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-mu_gamma-', run_date, '.pdf'))

#mu_theta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l2,
                   params = 'mu_theta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-mu_theta-', run_date, '.pdf'))

#mu_gamma2 ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l2,
                   params = 'mu_gamma2',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-mu_gamma2-', run_date, '.pdf'))

#mu_theta2 ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit_l2,
                   params = 'mu_theta2',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-mu_theta2-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l2,
                   params = 'sigma_alpha',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-sigma_alpha-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l2,
                   params = 'sigma_beta',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-sigma_beta-', run_date, '.pdf'))

#lambda_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l2,
                   params = 'lambda_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-lambda_proc-', run_date, '.pdf'))

#kappa_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l2,
                   params = 'kappa_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-kappa_proc-', run_date, '.pdf'))

#sigma_gt ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l2,
                   params = 'sigma_gt',
                   ISB = 'FALSE',
                   exact = FALSE,
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-sigma_gt-', run_date, '.pdf'))

#sigma_gt2 ~ HN(0, 2)
PR_p <- rnorm(10000, 0, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit_l2,
                   params = 'sigma_gt2',
                   ISB = 'FALSE',
                   exact = FALSE,
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-sigma_gt2-', run_date, '.pdf'))

#nu ~ gamma(2, 0.1)
PR <- rgamma(10000, 2, 0.1)
MCMCvis::MCMCtrace(fit_l2,
                   params = 'nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-l2-ss2-trace-nu-', run_date, '.pdf'))

