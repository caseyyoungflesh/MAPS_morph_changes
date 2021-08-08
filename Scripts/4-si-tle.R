####################
# Changes in size idx over time/lat/elev
#
####################


# set dir ----------------------------------------------------------------

run_date <- '2021-07-28'
maps_master_date <- '2021-07-28'

dir <- 'XXXX'


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)
library(MCMCvis)


# process input data -------------------------------------------------------------

#read in data
setwd(paste0(dir, '/Data'))
mdata <- readRDS(paste0('MAPS-master-', maps_master_date, '.rds'))

#add cn_id for each species, station
mdata2 <- mdata %>%
  dplyr::group_by(sp_id, station) %>%
  dplyr::mutate(cn_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sp_id, cn_id)

#sp_id for each cn_id
cn_sp <- unique(mdata2[,c('sp_id', 'cn_id')])


# processing --------------------------------------------------------------

#response <- mdata2$wingy_idx
response <- mdata2$size_idx

#scale factor
scf_YR <- 100
scf_LAT <- 100
scf_ELEV <- 10000

#center covariates and scale so param estimates are near unit scale for sampling efficiency
sc_year <- scale(mdata2$year, scale = FALSE)[,1] / scf_YR
lat_cn_id <- unique(mdata2[,c('cn_id', 'lat')])
sc_lat <- scale(lat_cn_id$lat, scale = FALSE)[,1] / scf_LAT
elev_cn_id <- unique(mdata2[,c('cn_id', 'GMTED_elev')])
sc_elev <- scale(elev_cn_id$GMTED_elev, scale = FALSE)[,1] / scf_ELEV

DATA <- list(N = NROW(mdata2),
             Nsp = length(unique(mdata2$sp_id)),
             Nsc = length(unique(mdata2$cn_id)),
             y = response,
             sp = mdata2$sp_id, #species
             cn_id = mdata2$cn_id,
             cn_sp = cn_sp$sp_id,
             year = sc_year,
             lat = sc_lat,
             elev = sc_elev,
             pro_data = mdata2,
             scf_YR = scf_YR,
             scf_LAT = scf_LAT,
             scf_ELEV = scf_ELEV)


# Call Stan model --------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.92
TREE_DEPTH <- 10
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 8000

fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/morph-tle.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'xi',
                            'eta',
                            'gamma',
                            'theta',
                            'mu_alpha',
                            'mu_gamma',
                            'mu_theta',
                            'mu_eta',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma_xi',
                            'sigma_eta',
                            'sigma_gt',
                            'Rho_gt',
                            'sigma_proc',
                            'lambda_xi',
                            'kappa_xi',
                            'lambda_proc',
                            'kappa_proc',
                            'nu'),
                   #'y_rep'), 
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
  beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')
  xi_ch <- MCMCvis::MCMCchains(fit, params = 'xi')
  sigma_proc_ch <- MCMCvis::MCMCchains(fit, params = 'sigma_proc')
  nu_ch <- MCMCvis::MCMCchains(fit, params = 'nu')
  
  #random iterations to draw from posteriors
  idx <- sample(1:NROW(xi_ch), Nsim)
  
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
        beta_ch[idx[i], DATA$cn_id[j]] * DATA$year[j] + 
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
                  file_name = paste0('si-tle-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('si-tle-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('si-tle-fit-', run_date),
                  add_obj = list(DATA, ppc_sim),
                  add_obj_names = c(paste0('si-tle-data-', run_date),
                                    paste0('si-ppc-sim-', run_date)),
                  cp_file = c('Model_files/morph-tle.stan', 
                              '4-si-tle.R'),
                  cp_file_names = c(paste0('morph-tle-', run_date, '.stan'),
                                    paste0('4-si-tle-', run_date, '.R')))


# PPO ---------------------------------------------------------------------

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Results/si-tle-', run_date, '/Figures')),
       dir.create(paste0(dir, 'Results/si-tle-', run_date, '/Figures')),
       FALSE)

setwd(paste0(dir, 'Results/si-tle-', run_date, '/Figures'))

#mu_gamma ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-mu_gamma-', run_date, '.pdf'))

#mu_theta ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_theta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-mu_theta-', run_date, '.pdf'))

#mu_eta ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_eta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-mu_eta-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-sigma_beta-', run_date, '.pdf'))

#sigma_gt[1] ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gt[1]',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-sigma_gt[1]-', run_date, '.pdf'))

#sigma_gt[2] ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gt[2]',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-sigma_gt[2]-', run_date, '.pdf'))

#sigma_eta ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_eta',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-sigma_eta-', run_date, '.pdf'))

#lambda_xi ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'lambda_xi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-lambda_xi-', run_date, '.pdf'))

#kappa_xi ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'kappa_xi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-kappa_xi-', run_date, '.pdf'))

#lambda_proc ~ N(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'lambda_proc',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-lambda_proc-', run_date, '.pdf'))

#kappa_proc ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'kappa_proc',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-kappa_proc-', run_date, '.pdf'))

#mu_alpha ~ N(0, 1)
PR <- rnorm(10000, 0, 1)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_alpha',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-mu_alpha-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_alpha',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-sigma_alpha-', run_date, '.pdf'))

#nu ~ gamma(2, 0.1)
PR <- rgamma(10000, 2, 0.1)
MCMCvis::MCMCtrace(fit,
                   params = 'nu',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('si-tle-trace-nu-', run_date, '.pdf'))


# PPC plots ----------------------------------------------------------------

#overall
pdf('PPC.pdf', height = 5, width = 5)
plot(density(ppc_sim$y), lwd = 2, main = 'PPC', xlab = 'Value')
for (i in 1:NROW(ppc_sim$y_rep))
{
  lines(density(ppc_sim$y_rep[i,]), col = rgb(1,0,0,0.1))
}
dev.off()

#species-specific PPC and r^2
sp_PPC_fun <- function(ppc_obj, sci_name, return)
{
  #get idx for species 
  sp_idx <- which(ppc_sim$sci_name == sci_name)
  
  #data points just for that species
  ty <- ppc_sim$y[sp_idx]
  ty_rep <- ppc_sim$y_rep[,sp_idx]
  tmu <- ppc_sim$mu[,sp_idx]
  
  plot(density(ty), lwd = 2, main = paste0('PPC - ', sci_name), xlab = 'Value')
  for (i in 1:NROW(ty_rep))
  {
    lines(density(ty_rep[i,]), col = rgb(1,0,0,0.1))
  }

}

#run function PPC
pdf('PPC_species.pdf', height = 5, width = 5)
par(mfrow = c(3,3))
usp <- sort(unique(ppc_sim$sci_name))
for (i in 1:length(usp))
{
  sp_PPC_fun(ppc_obj = ppc_sim, 
             sci_name = usp[i], 
             type = 'PPC', 
             return = FALSE)
}
dev.off()
