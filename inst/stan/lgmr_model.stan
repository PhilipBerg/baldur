functions {
  vector reg_function(vector x, vector theta, real I, real I_L, real S, real S_L, int N) {
    vector[N] exp_beta  = .001*exp(theta.*(I_L - S_L*x));
              exp_beta +=      exp(I - S*x);
    return exp_beta;
  }
}

data {
  int<lower=0> N;       // Number of observations
  vector<lower=0>[N] y; // sd
  vector[N] x;          // mean
}

transformed data {
  real v_y = variance(y);                 // for NRMSE
  vector[N] x_star = (x - mean(x))/sd(x); // f(y_bar)
}
parameters {
  real<lower = 0> alpha;    // Shape
  real<lower = 0> alpha_mu; // Shape mean hyperprior
  real I;                   // Intercep common
  real I_L;                 // Intercept Latent
  vector<lower = 0>[2] eta; // For S,S_L
  vector<lower = 0, upper = 1>[N] theta; // Mixture paramater
}

transformed parameters {
  real<lower = 0> S     = eta[1];    // Slope common
  real<lower = 0> S_L   = eta[2]*.1; // Slope Latent
}

model {
  alpha         ~ exp_mod_normal(alpha_mu, 1, .1); //Priors
  alpha_mu      ~ normal(50, 10);
  I             ~ std_normal();
  I_L           ~ skew_normal(2, 15,  35);
  eta           ~ std_normal();
  theta         ~ beta(.5, .5);
  {
    vector[N] exp_beta = reg_function(x_star, theta, I, I_L, S, S_L, N);
    y ~ gamma(alpha, alpha ./ exp_beta); // Likelihood
  }
}

generated quantities {
  real nrmse; // NRMSE calculations
  {
    vector[N] se = reg_function(x_star, theta, I, I_L, S, S_L, N);
    se -= y;
    se = se^2;
    nrmse = mean(se) / v_y;
  }
  nrmse = sqrt(nrmse);
}
