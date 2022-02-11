// SI ~ (lagged) temporal var in temp

data {
int<lower=0> N;                     // number of obs
int<lower=0> Nsp;                   // number of species
int<lower=0> Nsc;                   // number of species/stations
real y[N];                          // response
int<lower=1, upper=Nsc> cn_id[N];   // species/station ids for each obs
int<lower=1, upper=Nsp> cn_sp[Nsc]; // species ids for each cn_id
int<lower=1, upper=Nsp> sp[N];      // species ids for each y
vector[N] temp;                     // temp
vector[Nsc] mn_st_temp;             // mean station temp
}

parameters {
// real<lower=0> sigma_proc;
vector<lower=0>[Nsp] sigma_proc;
real<lower=0> lambda_proc;
real<lower=0> kappa_proc;
vector[Nsc] alpha_raw;
vector[Nsc] beta_raw;
real mu_gamma;
real mu_theta;
real mu_gamma2;
real mu_theta2;
// real mu_kappa;
// real<lower=0> sigma_kappa;
// vector[Nsp] kappa;
// vector<offset = mu_kappa, multiplier = sigma_kappa>[Nsp] kappa;
real<lower=0> sigma_alpha;
real<lower=0> sigma_beta;
vector<lower=0>[2] sigma_gt;
cholesky_factor_corr[2] L_Rho_gt;   // cholesky factor of corr matrix
matrix[2, Nsp] z_gt;                // z-scores
vector<lower=0>[2] sigma_gt2;
cholesky_factor_corr[2] L_Rho_gt2;   // cholesky factor of corr matrix
matrix[2, Nsp] z_gt2;                // z-scores
real<lower=0> nu;
}

transformed parameters {
vector[N] mu;
vector[Nsc] alpha;
vector[Nsc] mu_alpha;
vector[Nsc] beta;
vector[Nsc] mu_beta;
vector[Nsp] gamma;
vector[Nsp] theta;
matrix[Nsp, 2] gt;                 // gamma, and theta
matrix[2, 2] Rho_gt;               // corr matrix
vector[Nsp] gamma2;
vector[Nsp] theta2;
matrix[Nsp, 2] gt2;                 // gamma, and theta
matrix[2, 2] Rho_gt2;               // corr matrix


// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
gt2 = (diag_pre_multiply(sigma_gt2, L_Rho_gt2) * z_gt2)';
// implies Rho = L_Rho * L_Rho';
Rho_gt2 = multiply_lower_tri_self_transpose(L_Rho_gt2);

gamma2 = mu_gamma2 + gt2[,1];
theta2 = mu_theta2 + gt2[,2];

mu_alpha = gamma2[cn_sp] + theta2[cn_sp] .* mn_st_temp;
alpha = alpha_raw * sigma_alpha + mu_alpha;

// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
gt = (diag_pre_multiply(sigma_gt, L_Rho_gt) * z_gt)';
// implies Rho = L_Rho * L_Rho';
Rho_gt = multiply_lower_tri_self_transpose(L_Rho_gt);

gamma = mu_gamma + gt[,1];
theta = mu_theta + gt[,2];

// gamma is effect of temp on idx at mean mn_st_z
// theta is effect of mn station temp on temp effect
mu_beta = gamma[cn_sp] + theta[cn_sp] .* mn_st_temp;

beta = beta_raw * sigma_beta + mu_beta;
mu = alpha[cn_id] + beta[cn_id] .* temp;
}

model {
mu_gamma ~ normal(0, 2);
mu_theta ~ normal(0, 2);
mu_gamma2 ~ normal(0, 2);
mu_theta2 ~ normal(0, 2);
// mu_kappa ~ normal(0, 2);
// sigma_kappa ~ normal(0, 2);
sigma_alpha ~ normal(0, 2);
sigma_beta ~ normal(0, 2);
lambda_proc ~ std_normal();
kappa_proc ~ std_normal();
sigma_gt ~ normal(0, 2);
sigma_gt2 ~ normal(0, 2);
nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// non-centered
alpha_raw ~ std_normal();
beta_raw ~ std_normal();

// centered
// kappa ~ normal(mu_kappa, sigma_kappa);
sigma_proc ~ normal(lambda_proc, kappa_proc);

to_vector(z_gt) ~ std_normal();
L_Rho_gt ~ lkj_corr_cholesky(1);

to_vector(z_gt2) ~ std_normal();
L_Rho_gt2 ~ lkj_corr_cholesky(1);

y ~ student_t(nu, mu, sigma_proc[sp]); 
}
