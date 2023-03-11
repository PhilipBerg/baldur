functions {
  vector reg_function(vector x, vector p, real I, real I_L, real S, real S_L, int N) {
    vector[N] exp_beta  = .001*exp(p.*(I_L - S_L*x));
              exp_beta +=      exp(I - S*x);
    return exp_beta;
  }
}

data {
  int<lower=0> N;       // Number of observations
  vector<lower=0>[N] y; // observations
  vector[N] x;          // regressor
}

transformed data {
  real v_y = variance(y);
  vector[N] x_star = (x - mean(x))/sd(x);
}
parameters {
  real<lower = 0> alpha;
  real<lower = 0> alpha_mu;
  real I;
  real I_L;
  vector<lower = 0>[2] eta;
  vector<lower = 0, upper = 1>[N] p; // Fuzzy mixture paramater
}

transformed parameters {
  real<lower = 0> S     = eta[1];
  real<lower = 0> S_L   = eta[2]*.1;
}

model {
  alpha         ~ exp_mod_normal(alpha_mu, 1, .1);
  alpha_mu      ~ normal(50, 10);
  I             ~ std_normal();
  I_L           ~ skew_normal(2, 15,  35);
  eta           ~ std_normal();
  p             ~ beta(.5, .5);
  {
    vector[N] exp_beta = reg_function(x_star, p, I, I_L, S, S_L, N);
    y ~ gamma(alpha, alpha ./ exp_beta); // Likelihood
  }
}

generated quantities {
  real nrmse;
  {
    vector[N] se = reg_function(x_star, p, I, I_L, S, S_L, N);
    se -= y;
    se = se^2;
    nrmse = mean(se) / v_y;
  }
  nrmse = sqrt(nrmse);
}
