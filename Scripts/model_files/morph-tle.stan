// morph ~ time + lat + elev

data {
int<lower=0> N;                     // number of obs
int<lower=0> Nsp;                   // number of species
int<lower=0> Nsc;                   // number of species/stations
real y[N];                          // SI or WI
int<lower=1, upper=Nsp> sp[N];      // species ids for each obs
int<lower=1, upper=Nsc> cn_id[N];   // species/station ids for each obs
int<lower=1, upper=Nsp> cn_sp[Nsc]; // species ids for each cn_id
vector[N] year;                     // time
vector[Nsc] lat;                    // latitude
vector[Nsc] elev;                   // elev
}

parameters {
vector[Nsc] xi_raw;
vector<lower=0>[Nsp] sigma_xi;
vector[Nsc] beta_raw;
real<lower=0> sigma_beta;
real mu_alpha;
real<lower=0> sigma_alpha;
vector<offset = mu_alpha, multiplier = sigma_alpha>[Nsp] alpha;
vector<lower=0>[Nsp] sigma_proc;
// real lambda;
real<lower=0> lambda_xi;
real<lower=0> kappa_xi;
real<lower=0> lambda_proc;
real<lower=0> kappa_proc;
real mu_gamma;
real mu_theta;
vector<lower=0>[2] sigma_gt;
cholesky_factor_corr[2] L_Rho_gt;             // cholesky factor of corr matrix
matrix[2, Nsp] z_gt;                          // z-scores
real mu_eta;
real<lower=0> sigma_eta;
vector<offset = mu_eta, multiplier = sigma_eta>[Nsp] eta;
real<lower=0> nu;
}

transformed parameters {
vector[N] mu;
vector[Nsc] xi;
vector[Nsc] mu_xi;
vector[Nsc] beta;
vector[Nsp] gamma;
vector[Nsp] theta;
matrix[Nsp, 2] gt;                 // gamma, and theta
matrix[2, 2] Rho_gt;               // corr matrix

// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
gt = (diag_pre_multiply(sigma_gt, L_Rho_gt) * z_gt)';
// implies Rho = L_Rho * L_Rho';
Rho_gt = multiply_lower_tri_self_transpose(L_Rho_gt);

gamma = mu_gamma + gt[,1];
theta = mu_theta + gt[,2];

mu_xi = gamma[cn_sp] .* lat + theta[cn_sp] .* elev;

// implies xi_jk ~ normal(mu_xi_jk, sigma_xi)
xi = xi_raw .* sigma_xi[cn_sp] + mu_xi;

// implies beta_jk ~ normal(eta_k, sigma_beta)
beta = beta_raw * sigma_beta + eta[cn_sp];

mu = alpha[sp] + beta[cn_id] .* year + xi[cn_id];
}

model {
mu_gamma ~ normal(0, 5);
mu_theta ~ normal(0, 5);
mu_eta ~ normal(0, 5);
sigma_eta ~ normal(0, 5);
sigma_beta ~ normal(0, 5);
sigma_gt ~ normal(0, 5);        // std dev gamma, theta
lambda_xi ~ std_normal();
kappa_xi ~ std_normal();
lambda_proc ~ std_normal();
kappa_proc ~ std_normal();
mu_alpha ~ std_normal();
sigma_alpha ~ std_normal();
nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// centered
sigma_xi ~ normal(lambda_xi, kappa_xi);
sigma_proc ~ normal(lambda_proc, kappa_proc);

// non-centered
alpha ~ normal(mu_alpha, sigma_alpha);
eta ~ normal(mu_eta, sigma_eta);
xi_raw ~ std_normal();
beta_raw ~ std_normal();

to_vector(z_gt) ~ std_normal();
L_Rho_gt ~ lkj_corr_cholesky(1);

y ~ student_t(nu, mu, sigma_proc[sp]); 
}
