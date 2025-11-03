//
// DyNAMCoord1REV05.stan
// DyNAM Choice Coordination with P covariates and one Random Effect - Map-Reduce Parallelization
// it's a multilevel model with RE by groups and category-specific fixed effects
// non centered parametrization, full mixed formulation
//
// Performance improvements over REV01:
// - Pre-computed indices in transformed data to eliminate function calls
// - Reduced memory allocations in loops
// - Optimized indexing patterns for better cache performance
// - Memory buffer reuse across events
// - Standardized predictors with standard normal priors for better sampling
// - Category-specific fixed effects for flexible modeling
// - Map-reduce parallelization for within-chain parallelization
//
// Based on DyNAMCoord1REV04.stan with map-reduce parallelization

functions {
  // transform the rowIndex and colIndex to a linear index
  int index_pos(int rowIndex, int colIndex, int nActors) {
    // index transformation works for square matrices (one-mode)
    int rowLower = rowIndex <= colIndex ? rowIndex : colIndex;
    int colHigher = rowIndex <= colIndex ? colIndex : rowIndex;

    return((nActors * (rowLower - 1)) + colHigher - choose(rowLower + 1, 2));
  }

  real partial_sum_dynam_coord(
    array[] int event_subset,
    int start,
    int end,
    matrix X,
    vector Z,
    array[] int group,
    array[] int category,
    array[,] int X_indices,
    array[,] int dyad_indices,
    array[] int event_nDyads,
    array[] int observed_pos,
    array[] int event_start_flat_row,
    array[] int nA,
    matrix beta_std,
    vector gamma_raw,
    real sigma
  ) {
    real log_lik = 0.0;

    for (t in event_subset) {
      int n_choices = nA[t] - 1;
      int current_nDyads = event_nDyads[t];
      vector[current_nDyads] logLikDyad = rep_vector(0.0, current_nDyads);
      vector[n_choices] lprobs = rep_vector(0.0, n_choices);
      int flat_row = event_start_flat_row[t];
      array[n_choices] int x_row_idx; 
      vector[n_choices] xb;

      for (row in 1:nA[t]) {
        x_row_idx = X_indices[flat_row, 1:n_choices];
        for (i in 1:n_choices)
          xb[i] = X[x_row_idx[i]] * beta_std[category[x_row_idx[i]]]' +
            Z[x_row_idx[i]] * (sigma * gamma_raw[group[x_row_idx[i]]]);

        lprobs = log_softmax(xb);
        logLikDyad[dyad_indices[flat_row, 1:n_choices]] += lprobs;
        flat_row += 1;
      }

      log_lik += logLikDyad[observed_pos[t]] - log_sum_exp(logLikDyad);
    }

    return log_lik;
  }
}

data {
  int N; // number of events * actors * (actors - 1)
  int E; // number of events
  int MxAct; // maximum number of actors in an event

  int P; // number of covariates
  int G; // number of groups
  int C; // number of categories

  matrix[N, P] X_raw; // Raw (unstandardized) data matrix with covariates, include RE fix part
  vector[N] Z_raw; // Raw (unstandardized) Random Effect design vector
  array[N] int<lower = 1, upper = G> group; // Group assignments
  array[N] int<lower = 1, upper = C> category; // Category assignment for each observation

  // number of actors in an event
  array[E] int<lower = 1, upper = MxAct> nA;

  // row col observed event
  array[E] int<lower = 1, upper = MxAct> rowOE;
  array[E] int<lower = 1, upper = MxAct> colOE;

  int<lower=1> grain_size; // Grain size for map_reduce parallelization
}

transformed data {
  // Standardize predictors
  matrix[N, P] X;
  vector[P] X_means;
  vector[P] X_sds;

  for (p in 1:P) {
    X_means[p] = mean(X_raw[, p]);
    X_sds[p] = sd(X_raw[, p]);
    X[, p] = (X_raw[, p] - X_means[p]) / X_sds[p];
  }

  // Standardize Z
  vector[N] Z;
  real Z_mean;
  real Z_sd;

  Z_mean = mean(Z_raw);
  Z_sd = sd(Z_raw);
  Z = (Z_raw - Z_mean) / Z_sd;

  // Pre-compute all indices to eliminate function calls in model evaluation

  // Total number of choice sets across all events
  int total_actors = sum(nA);

  // Pre-computed indices for X matrix access
  array[total_actors, MxAct - 1] int X_indices;

  // Pre-computed dyad indices for likelihood accumulation
  array[total_actors, MxAct - 1] int dyad_indices;

  // Pre-computed observed event positions
  array[E] int observed_pos;

  // Calculate maximum dyads across all events
  int max_dyads = choose(MxAct, 2);

  // Pre-compute exact dyad counts for each event (avoids repeated choose() calls)
  array[E] int event_nDyads;

  // Event boundaries for map-reduce
  array[E] int event_start_flat_row;
  array[E] int event_end_flat_row;

  // Pre-compute all indices
  {
    int flat_idx = 1;
    int flat_row = 1;

    for (t in 1:E) {
      // Pre-compute observed positions
      observed_pos[t] = index_pos(rowOE[t], colOE[t], nA[t]);
      // dyads in event
      event_nDyads[t] = choose(nA[t], 2);

      // Event boundaries for map-reduce
      event_start_flat_row[t] = flat_row;

      for (rowChoice in 1:nA[t]) {
        // Pre-compute all indices for this choice set
        int local_idx = 1;
        for (i in 1:nA[t]) {
          if (i != rowChoice) {
            X_indices[flat_row, local_idx] = flat_idx;
            dyad_indices[flat_row, local_idx] = index_pos(rowChoice, i, nA[t]);
            local_idx += 1;
          }
          flat_idx += 1;
        }
        flat_row += 1;
      }

      // End boundary for this event
      event_end_flat_row[t] = flat_row - 1;
    }
    print("flat_row = ", flat_row);
    print("flat_idx = ", flat_idx);
  }
}

parameters {
  matrix[C, P] beta_std;  // Category-specific standardized coefficients (C × P matrix)
  real<lower=0> sigma; // variance RE
  vector[G] gamma_raw;  // group Random Effect
}

transformed parameters {
  // Back-transform category-specific coefficients to original scale
  matrix[C, P] beta;
  for (c in 1:C) {
    beta[c, ] = beta_std[c, ] ./ X_sds';
  }

  // Back-transform random effects to original scale (corrected: division)
  vector[G] gamma = (sigma * gamma_raw) / Z_sd;
}

model {
  // Enhanced priors for category-specific coefficients
  for (c in 1:C) {
    target += std_normal_lpdf(beta_std[c, ]);  // Standard normal for each category
  }
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);

  // Map-reduce parallelized likelihood computation
  array[E] int event_indices = linspaced_int_array(E, 1, E);

  target += reduce_sum(partial_sum_dynam_coord,
                        event_indices,
                        grain_size,
                        X,
                        Z,
                        group,
                        category,
                        X_indices,
                        dyad_indices,
                        event_nDyads,
                        observed_pos,
                        event_start_flat_row,
                        nA,
                        beta_std,
                        gamma_raw,
                        sigma);
}
