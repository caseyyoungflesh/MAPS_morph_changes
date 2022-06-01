####################
# si ~ temp space
#
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


# read in data -------------------------------------------------------------

#read in morph fit data
setwd(paste0(dir, '/Results/si-tle-', idx_tle_date))
morph_data <- readRDS(paste0('si-tle-data-', idx_tle_date, '.rds'))

#read in env data
setwd(paste0(dir, '/Data/daymet/processed'))
daymet <- readRDS(paste0('daymet-master-', daymet_process_date, '.rds')) %>%
  mutate(MJJ_tmax = (May_tmax + June_tmax + July_tmax) / 3) #using mean here does not do what might be expected


# process temp data ----------------------------------------------------------------

#data from 5-tle
mdata <- morph_data$pro_data

#merge env with MAPS
mdata2 <- dplyr::select(mdata, sci_name, sp_id, station, cn_id, lat, lng, year, GMTED_elev, size_idx) %>%
  dplyr::left_join(daymet, by = c('station', 'year')) %>%
  dplyr::arrange(sp_id, cn_id)

#sp_id for each cn_id
cn_sp <- unique(mdata2[,c('sp_id', 'cn_id')])


# processing --------------------------------------------------------------

#scale factor for temp
scf_temp <- 10

#mean temp (at each species/station) - scaled across dataset
mn_st_temp <- mdata2 %>%
  dplyr::group_by(sp_id, cn_id) %>%
  dplyr::summarize(mn_MJJ_tmax = mean(MJJ_tmax)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sc_mean_temp = scale(mn_MJJ_tmax, scale = FALSE)[,1] / scf_temp)

Nsp <- length(unique(mdata2$sp_id))

DATA <- list(N = NROW(mdata2),
             Nsp = Nsp, #number of species
             Nsc = length(unique(mdata2$cn_id)), #number of y data points
             y = mdata2$size_idx,
             sp = mdata2$sp_id, #species
             cn_id = mdata2$cn_id,
             cn_sp = cn_sp$sp_id,
             sc_mean_temp = mn_st_temp$sc_mean_temp,
             scf_temp = scf_temp,
             pro_data = mdata2)


# Call Stan model ---------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#model run settings
DELTA <- 0.96
TREE_DEPTH <- 12
STEP_SIZE <- 0.02
CHAINS <- 4
ITER <- 8000

#run model
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/morph-temp-space.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'mu_alpha',
                            'sigma_alpha',
                            'xi',
                            'sigma_xi',
                            'lambda_xi',
                            'kappa_xi',
                            'beta',
                            'mu_beta',
                            'sigma_beta',
                            'nu',
                            'sigma_proc',
                            'lambda_proc',
                            'kappa_proc'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# PPC sim --------------------------------------------------------

#returns list with y, cn_id, sp, y_rep, and mu
PPC_fun <- function(Nsim = 100, 
                    seed = 1)
{
  print('Extracting posteriors')
  alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
  xi_ch <- MCMCvis::MCMCchains(fit, params = 'xi')
  sigma_proc_ch <- MCMCvis::MCMCchains(fit, params = 'sigma_proc')
  nu_ch <- MCMCvis::MCMCchains(fit, params = 'nu')
  
  #random iterations to draw from posteriors
  idx <- sample(1:NROW(alpha_ch), Nsim)
  
  print('Simulating data')
  y_rep <- matrix(NA, nrow = Nsim, ncol = DATA$N)
  mu <- matrix(NA, nrow = Nsim, ncol = DATA$N)
  
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = (DATA$N * Nsim), style = 3)
  counter <- 0
  for (i in 1:length(idx))
  {
    for (j in 1:DATA$N)
    {
      #i <- 1
      #j <- 1
      mu[i,j] <- alpha_ch[idx[i], DATA$sp[j]] + 
        xi_ch[idx[i], DATA$cn_id[j]]
      
      eps <- rt(n = 1, df = nu_ch[idx[i]]) * sigma_proc_ch[idx[i], DATA$sp[j]]
      y_rep[i,j] <- mu[i,j] + eps
      
      #progress bar
      setTxtProgressBar(pb, counter)
      counter <- counter + 1
    }
  }
  
  return(list(y = DATA$y, 
              sci_name = DATA$pro_data$sci_name, 
              sp = DATA$sp, 
              cn_id = DATA$cn_id, 
              lat = DATA$pro_data$lat,
              lng = DATA$pro_data$lng,
              y_rep = y_rep, 
              mu = mu))
}

#run function
ppc_sim <- PPC_fun(Nsim = 100)


# save summary ------------------------------------------------------------

setwd(paste0(dir, '/Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('si-temp-space-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('si-temp-space-', run_date),
                  add_field = 'MJJ tmax',
                  add_field_names = 'Temp data',
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('si-temp-space-fit-', run_date),
                  add_obj = list(DATA, ppc_sim),
                  add_obj_names = c(paste0('si-temp-space-data-', run_date),
                                    paste0('si-temp-space-ppc-sim-', run_date)),
                  cp_file = c('Model_files/morph-temp-space.stan', 
                              '6-si-temp-space.R'),
                  cp_file_names = c(paste0('si-temp-space', run_date, '.stan'),
                                    paste0('6-si-temp-space-', run_date, '.R')))


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

setwd(paste0(dir, '/Results/si-temp-space-', run_date))
ppc_plot_fun(sim_obj = ppc_sim, file_name = 'PPC.pdf')


# PPO ----------------------------------------------------------------

#mu_alpha ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_alpha',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-mu_alpha-', run_date, '.pdf'))


#sigma_alpha ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_alpha',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-sigma_alpha-', run_date, '.pdf'))

#mu_beta ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_beta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-mu_beta-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-sigma_beta-', run_date, '.pdf'))

#lambda_xi ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'lambda_xi',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-lambda_xi-', run_date, '.pdf'))

#kappa_xi ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'kappa_xi',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-kappa_xi-', run_date, '.pdf'))

#lambda_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'lambda_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-lambda_proc-', run_date, '.pdf'))

#kappa_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'kappa_proc',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-kappa_proc-', run_date, '.pdf'))

#nu ~ gamma(2, 0.1)
PR <- rgamma(10000, 2, 0.1)
MCMCvis::MCMCtrace(fit,
                   params = 'nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-temp-space-trace-nu-', run_date, '.pdf'))
