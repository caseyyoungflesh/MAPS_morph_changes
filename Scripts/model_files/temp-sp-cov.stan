// SI ~ spatial var in temp

data {
int<lower=0> N;                       // number of obs
real y[N];                            // observed
real<lower=0> sd_y[N];                            // uncertainty in observed
vector[N] mn_temp;                  // mean temp experienced by species
matrix[N, N] R;                   // scaled phylo distance matrix (from phylogeny)
matrix[N, N] I;                   // identity matrix
}

parameters {
vector[N] mu_y;
real<lower=0> sigma;
real alpha;
real beta;
real<lower=0, upper=1> lambda;
}

transformed parameters {
vector[N] mu;
matrix[N, N] R_lambda;
matrix[N, N] S;

// Pagel's lambda
R_lambda = lambda * R + (1 - lambda) * I;
S = R_lambda * sigma;

mu = alpha + beta * mn_temp;
}

model {
alpha ~ normal(0, 5);
beta ~ normal(0, 5);
sigma ~ normal(0, 5);
lambda ~ uniform(0, 1);

mu_y ~ multi_normal(mu, S); 
y ~ normal(mu_y, sd_y);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu_y, sd_y);
}
