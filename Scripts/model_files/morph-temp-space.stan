// SI ~ spatial var in temp

data {
int<lower=0> N;                       // number of obs
int<lower=0> Nsp;                     // number of species
int<lower=0> Nsc;                     // number of species/stations
real y[N];                            // response
int<lower=1, upper=Nsp> sp[N];        // species ids for each obs
int<lower=1, upper=Nsc> cn_id[N];     // species/station ids for each obs
int<lower=1, upper=Nsp> cn_sp[Nsc];   // species ids for each cn_id
vector[Nsc] sc_mean_temp;             // mean temp at each station
}

parameters {
real mu_alpha;
real<lower=0> sigma_alpha;
vector<offset = mu_alpha, multiplier = sigma_alpha>[Nsp] alpha;
// vector[Nsp] alpha;
real<lower=1> nu;
vector[Nsc] xi_raw;
// vector[Nsc] xi;
vector<lower=0>[Nsp] sigma_xi;
real<lower=0> lambda_xi;
real<lower=0> kappa_xi;
vector<lower=0>[Nsp] sigma_proc;
real<lower=0> lambda_proc;
real<lower=0> kappa_proc;
real mu_beta;
real<lower=0> sigma_beta;
vector<offset = mu_beta, multiplier = sigma_beta>[Nsp] beta;
// vector[Nsp] beta;
}

transformed parameters {
vector[N] mu;
vector[Nsc] xi;
vector[Nsc] mu_xi;

mu_xi = beta[cn_sp] .* sc_mean_temp;

// implies xi_jk ~ normal(mu_xi_jk, sigma_xi)
xi = xi_raw .* sigma_xi[cn_sp] + mu_xi;

mu = alpha[sp] + xi[cn_id];
}

model {
mu_alpha ~ normal(0, 5);
sigma_alpha ~ normal(0, 5);
mu_beta ~ normal(0, 5);
sigma_beta ~ normal(0, 5);
lambda_xi ~ std_normal();
kappa_xi ~ std_normal();
lambda_proc ~ std_normal();
kappa_proc ~ std_normal();
nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// centered
// alpha ~ normal(mu_alpha, sigma_alpha);
// beta ~ normal(mu_beta, sigma_beta);
sigma_proc ~ normal(lambda_proc, kappa_proc);
sigma_xi ~ normal(lambda_xi, kappa_xi);
// xi ~ normal(mu_xi, sigma_xi[cn_sp]);

// non-centered
alpha ~ normal(mu_alpha, sigma_alpha);
beta ~ normal(mu_beta, sigma_beta);
xi_raw ~ std_normal();

y ~ student_t(nu, mu, sigma_proc[sp]); 
}
