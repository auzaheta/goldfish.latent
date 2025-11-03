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
// DNRE_Q1_choice_int_centered.stan
// Multinomial choice model with P fixed effect covariates, one random effect,
// and category-specific interactions using standardized predictors.
// Centered parametrization with category-specific means for random effects.
//
// Performance improvements:
// - Standardized predictors with standard normal priors for better sampling
// - Category-specific interaction coefficients
// - Centered parametrization with hierarchical means
// - Back-transformed parameters in original scale for interpretability
// - No intercepts (multinomial identifiability)

functions {
  real partial_sum_lpmf(
    array[] int event_subset,
    int start,
    int end,
    array[] int start_choice,
    array[] int end_choice,
    array[] int chose_choice,
    matrix X_choice_std,
    vector Z_std,
    array[] int sender,
    array[] int interaction,
    array[] vector beta_choice_std,
    vector gamma_std) {
    real log_lik = 0.0;
    int size_slice;
    int start_event;
    int end_event;
    int chose_event;
    int interaction_event;

    for (t in event_subset) {
      start_event = start_choice[t];
      end_event = end_choice[t];
      chose_event = chose_choice[t] - start_event + 1;
      interaction_event = interaction[t];
      size_slice = end_event - start_event + 1;
      array[size_slice] int event_slice =
        linspaced_int_array(size_slice, start_event, end_event);

      vector[size_slice] xb_choice = X_choice_std[event_slice] * beta_choice_std[interaction_event] +
        Z_std[event_slice] .* gamma_std[sender[event_slice]];

      log_lik += xb_choice[chose_event] - log_sum_exp(xb_choice);
    }

    return log_lik;
  }
}

data {
  int<lower=1> N_choice;       // Total number of choices in all events
  int<lower=1> T_choice;       // Number of events
  int<lower=0> P_choice;       // Number of fixed-effect covariates (excluding random effect)

  matrix[N_choice, P_choice] X_choice_raw; // Raw (unstandardized) fixed effects data matrix

  // the starting, ending and chosen index observation for each event
  array[T_choice] int<lower=1, upper=N_choice> start_choice;
  array[T_choice] int<lower=1, upper=N_choice> end_choice;
  array[T_choice] int<lower=1, upper=N_choice> chose_choice;

  // random effects vars
  int<lower=1> A;              // Number of actors/groups
  int<lower=1> Q_choice;       // Number of random effects (must be 1)
  array[Q_choice] int<lower = 1, upper = P_choice> V1_choice; // index random effect covariate
  array[N_choice] int<lower=1, upper=A> sender;

  // interaction var for event
  int<lower=2> C; // number of categories
  array[T_choice] int<lower=1, upper=C> interaction;
  array[A] int<lower=1, upper=C> send_int; // sender's interaction category

  int<lower=1> grain_size;     // Grain size for map_reduce
}

transformed data {
  // Standardize fixed effects predictors
  matrix[N_choice, P_choice - 1] X_choice_std;
  vector[P_choice - 1] X_means;
  vector[P_choice - 1] X_sds;

  vector[N_choice] Z_std;
  real Z_mean;
  real Z_sd;

  int p_new = 1;
  for (p in 1:P_choice) {
    if (p != V1_choice[1]) {
      X_means[p_new] = mean(X_choice_raw[, p]);
      X_sds[p_new] = sd(X_choice_raw[, p]);
      X_choice_std[, p_new] = (X_choice_raw[, p] - X_means[p_new]) / X_sds[p_new];
      p_new += 1;
    } else {
      // Extract and standardize random effects predictor (assuming it's within X_choice)
      // Note: This assumes V1_choice[1] refers to a column in the FULL original matrix
      // For this centered model, we'll extract Z from the original raw data
      Z_mean = mean(X_choice_raw[, V1_choice[1]]);
      Z_sd = sd(X_choice_raw[, V1_choice[1]]);
      Z_std = (X_choice_raw[, V1_choice[1]] - Z_mean) / Z_sd;
    }
  }
}

parameters {
  array[C] vector[P_choice - 1] beta_choice_std; // Category-specific standardized fixed effects
  vector[C] mu_std;                      // Category-specific means for random effects
  real<lower=0> sigma;                       // Variance of the random effect
  vector[A] gamma_std;                       // Centered random effects
}

transformed parameters {
  // Back-transform category-specific coefficients to original scale
  array[C] vector[P_choice - 1] beta_choice;
  for (c in 1:C) {
    beta_choice[c] = beta_choice_std[c] ./ X_sds;
  }

  // Back-transform category means to original scale
  vector[C] mu = mu_std / Z_sd;

  // Back-transform random effects to original scale
  vector[A] gamma = (gamma_std - mu_std[send_int]) / Z_sd;
}

model {
  // Priors for standardized coefficients
  for (c in 1:C) {
    target += std_normal_lpdf(beta_choice_std[c]);
  }
  target += std_normal_lpdf(mu_std);  // Standard normal for category means
  target += exponential_lpdf(sigma | 1);

  // Hierarchical prior for random effects
  for (a in 1:A) {
    target += normal_lpdf(gamma_std[a] | mu_std[send_int[a]], sigma);
  }

  // Event-based parallelized likelihood computation
  array[T_choice] int event_indices = linspaced_int_array(T_choice, 1, T_choice);

  target += reduce_sum(
    partial_sum_lpmf,
    event_indices,
    grain_size,
    start_choice,
    end_choice,
    chose_choice,
    X_choice_std,
    Z_std,
    sender,
    interaction,
    beta_choice_std,
    gamma_std
  );
}
