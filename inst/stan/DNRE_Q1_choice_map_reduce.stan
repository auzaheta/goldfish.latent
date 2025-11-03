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
// DNRE_Q1_choice_map_reduce.stan
// Multinomial choice model with P fixed effect covariates and one random effect
// using map_reduce for within-chain parallelization.

functions {
  real partial_sum_lpmf(
    array[] int start_choice,
    int start,
    int end,
    array[] int end_choice,
    matrix X_choice,
    vector Z_temp,
    array[] int chose_choice,
    array[] int sender,
    vector beta_choice,
    vector gamma
 ) {
    real log_lik = 0.0;
    int t_index;
    int size_slice;
    int start_event;
    int end_event;
    int chose_event;

    for (t in 1:(end - start + 1)) {
      t_index = t + start - 1;
      start_event = start_choice[t];
      end_event = end_choice[t_index];
      chose_event = chose_choice[t_index] - start_event + 1;
      size_slice = end_event - start_event + 1;
      array[size_slice] int event_slice =
        linspaced_int_array(size_slice, start_event, end_event);
      
      vector[size_slice] xb_choice = X_choice[event_slice] * beta_choice +
        Z_temp[event_slice] .* gamma[sender[event_slice]];
      
      log_lik += xb_choice[chose_event] - log_sum_exp(xb_choice);
    }

    return log_lik;
  }
}

data {
  int<lower=1> N_choice;       // Total number of choices in all events
  int<lower=1> T_choice;       // Number of events
  int<lower=0> P_choice;       // Number of fixed-effect covariates

  matrix[N_choice, P_choice] X_choice;

  // the starting, ending and chosen index observation for each event
  array[T_choice] int<lower=1, upper=N_choice> start_choice;
  array[T_choice] int<lower=1, upper=N_choice> end_choice;
  array[T_choice] int<lower=1, upper=N_choice> chose_choice;

  // random effects vars
  int<lower=1> A;              // Number of actors/groups
  int<lower=1> Q_choice;       // Number of random effects (must be 1)
  matrix[N_choice, Q_choice] Z_choice;
  array[N_choice] int<lower=1, upper=A> sender;
  // array[A] int<lower=1, upper=T_choice> start_group;

  int<lower=1> grain_size;     // Grain size for map_reduce
}

transformed data {
  // array[A] int<lower=1, upper=T_choice> end_group;
  // for (g in 1:(A - 1)) {
  //   end_group[g] = start_group[g + 1] - 1;
  // }
  // end_group[A] = T_choice;

  vector[N_choice] Z_temp = to_vector(Z_choice);
}

parameters {
  vector[P_choice] beta_choice; // Fixed effects
  real<lower=0> sigma;          // Variance of the random effect
  vector[A] gamma_raw;          // Uncentered random effects
}

transformed parameters {
  vector[A] gamma = sigma * gamma_raw; // Centered random effects
}

model {
  // Priors
  target += std_normal_lpdf(beta_choice);
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);

  // Likelihood
  target += reduce_sum(
    partial_sum_lpmf,
    start_choice,
    grain_size,
    end_choice,
    X_choice,
    Z_temp,
    chose_choice,
    sender,
    beta_choice,
    gamma
  );
}
