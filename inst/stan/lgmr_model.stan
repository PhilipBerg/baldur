functions {
  // mu
  vector reg_function(
    vector x,
    vector theta,
    real I,
    real I_L,
    real S,
    real S_L,
    int N
    ) {
    vector[N] exp_beta  = .001*exp(theta .* (I_L - S_L*x));
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
  real v_y = variance(y);                      // for NRMSE
  vector[N] x_star = (x - mean(x))/sd(x);      // f(y_bar)
  vector[3] lb = [0, 0, negative_infinity()]'; // Boundaries normal coefs.
}
parameters {
  real<lower = 0> alpha;                 // Shape
  vector<lower = lb>[3] eta;             // For S, S_L, I
  vector<lower = 0, upper = 1>[N] theta; // Mixture parameter
  real I_L; // Intercept Latent
}

transformed parameters {
  real<lower = 0> S   = eta[1]; // Slope common
  real<lower = 0> S_L = eta[2]; // Slope Latent
  real I              = eta[3]; // Intercep common
}

model {
  //Priors
  alpha         ~ cauchy(0, 25);
  eta           ~ std_normal();
  I_L           ~ skew_normal(2, 15,  35);
  theta         ~ beta(1, 1);
  {
    vector[N] exp_beta = reg_function(x_star, theta, I, I_L, S, S_L, N);
    // Likelihood
    y ~ gamma(alpha, alpha ./ exp_beta);
  }
}

generated quantities {
  // NRMSE calculations
  real nrmse;
  {
    vector[N] se = reg_function(x_star, theta, I, I_L, S, S_L, N);
    se -= y;
    se  = square(se);
    nrmse = mean(se) / v_y;
  }
  nrmse = sqrt(nrmse);
}
