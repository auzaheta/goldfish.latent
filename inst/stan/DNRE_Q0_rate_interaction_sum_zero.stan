// Copyright (C) 2025, Alvaro Uzaheta - SNlab-ETH Zurich
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT
// License for more details.
//
// You should have received a copy of the MIT License along with this
// program. If not, see <https://opensource.org/licenses/MIT>.
//
// DNRE_Q0_rate_interaction_sum_zero.stan
// Rate model with P fixed effect covariates, random intercepts with sum-to-zero constraint
// and interactions for a categorical variable coded from 1 to C
// using QR-based sum-to-zero transformation for better identifiability.
//
// Performance improvements:
// - Standardized predictors with standard normal priors for better sampling
// - Sum-to-zero constrained random effects using QR transformation
// - Non-centered parametrization maintained for sampling efficiency
// - Back-transformed parameters in original scale for interpretability

functions {
  // QR-based sum-to-zero transformation from Martin Modrak's blog
  vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    Q_r = Q_r * inv_sqrt(1 - inv(N));
    return Q_r;
  }
  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector[N] x;
    real x_aux = 0;

    for(i in 1:N-1) {
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N-1];
    }
    x[N] = x_aux;
    return x;
  }

  real partial_sum_rate_lpmf(
    array[] int event_subset,
    int start,
    int end,
    array[] int start_rate,
    array[] int end_rate,
    array[] int chose_rate,
    vector timespan,
    array[] int is_dependent,
    real log_crude_rate,
    matrix X_rate_std,
    array[] int sender,
    array[] int interaction,
    array[] vector beta_rate_std,
    array[] real alpha_std,
    vector gamma
 ) {
    real log_lik = 0.0;
    int size_slice;
    int start_event;
    int end_event;
    int interaction_event;
    int chose_event;

    for (t in event_subset) {
      start_event = start_rate[t];
      end_event = end_rate[t];
      interaction_event = interaction[t];
      chose_event = chose_rate[t] - start_event + 1;
      size_slice = end_event - start_event + 1;
      array[size_slice] int event_slice =
        linspaced_int_array(size_slice, start_event, end_event);

      vector[size_slice] xb_rate =
        X_rate_std[event_slice] * beta_rate_std[interaction_event] +
        alpha_std[interaction_event] +
        gamma[sender[event_slice]] + log_crude_rate;

      if (timespan[t] > 0)
        log_lik += (is_dependent[t] ? xb_rate[chose_event] : 0) -
          timespan[t] * exp(log_sum_exp(xb_rate));
    }

    return log_lik;
  }
}

data {
  int<lower=1> N_rate;       // Total number of choices in all events
  int<lower=1> T_rate;       // Number of events
  int<lower=0> P_rate;       // Number of fixed-effect covariates

  matrix[N_rate, P_rate] X_rate_raw; // Raw (unstandardized) data matrix

  // the starting, ending and chosen index observation for each event
  array[T_rate] int<lower=1, upper=N_rate> start_rate;
  array[T_rate] int<lower=1, upper=N_rate> end_rate;
  array[T_rate] int<lower=0, upper=N_rate> chose_rate;

  // random effects vars
  int<lower=1> A;              // Number of actors/groups
  array[N_rate] int<lower=1, upper=A> sender;

  // rate information: time span and right-censored events
  vector<lower=0>[T_rate] timespan;
  array[T_rate] int<lower=0, upper=1> is_dependent;

  // interaction var for event
  int<lower=2> C; // number of categories
  array[T_rate] int<lower=1, upper=C> interaction;

  int<lower=1> grain_size;     // Grain size for map_reduce
}

transformed data {
  // Standardize predictors
  matrix[N_rate, P_rate] X_rate_std;
  vector[P_rate] X_means;
  vector[P_rate] X_sds;

  for (p in 1:P_rate) {
    X_means[p] = mean(X_rate_raw[, p]);
    X_sds[p] = sd(X_rate_raw[, p]);
    X_rate_std[, p] = (X_rate_raw[, p] - X_means[p]) / X_sds[p];
  }

  real log_crude_rate = log(T_rate / (N_rate * mean(timespan)));

  // QR transformation for sum-to-zero constraint
  vector[2 * A] Q_r = Q_sum_to_zero_QR(A);         // QR transformation matrix for sum-to-zero

}

parameters {
  array[C] vector[P_rate] beta_rate_std; // Standardized fixed effects
  array[C] real alpha_std;               // Category-specific intercepts
  real<lower=0> sigma;                   // Variance of the random effect
  vector[A-1] gamma_raw;                 // A-1 uncentered random effects (sum-to-zero)
}

transformed parameters {
  // Back-transform category-specific coefficients to original scale
  array[C] vector[P_rate] beta_rate;
  for (c in 1:C) {
    beta_rate[c] = beta_rate_std[c] ./ X_sds;
  }

  // Back-transform intercepts (include log_crude_rate AND standardization adjustment)
  array[C] real alpha;
  for (c in 1:C) {
    real adjustment = 0.0;
    for (p in 1:P_rate) {
      adjustment += beta_rate_std[c, p] * X_means[p] / X_sds[p];
    }
    alpha[c] = alpha_std[c] + log_crude_rate - adjustment;
  }

  // Sum-to-zero constrained random effects using QR transformation
  vector[A] gamma = sigma * sum_to_zero_QR(gamma_raw, Q_r);
}

model {
  // Priors for standardized coefficients
  for (c in 1:C) {
    target += std_normal_lpdf(beta_rate_std[c]);
    target += normal_lpdf(alpha_std[c] | 0, 4);
  }
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);  // A-1 dimensional prior

  // Event-based parallelized likelihood computation
  array[T_rate] int event_indices = linspaced_int_array(T_rate, 1, T_rate);

  target += reduce_sum(
    partial_sum_rate_lpmf,
    event_indices,
    grain_size,
    start_rate,
    end_rate,
    chose_rate,
    timespan,
    is_dependent,
    log_crude_rate,
    X_rate_std,
    sender,
    interaction,
    beta_rate_std,
    alpha_std,
    gamma
  );
}
