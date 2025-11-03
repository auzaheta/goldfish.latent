//
// DyNAMCoord1REV03.stan
// DyNAM Choice Coordination with P covariates and one Random Effect - Optimized Version with Standardization
// it's a multilevel model with RE by groups
// non centered parametrization, full mixed formulation
//
// Performance improvements over REV01:
// - Pre-computed indices in transformed data to eliminate function calls
// - Reduced memory allocations in loops
// - Optimized indexing patterns for better cache performance
// - Memory buffer reuse across events
// - Standardized predictors with standard normal priors for better sampling
//
// Based on DyNAMCoord1REV02.stan with standardization improvements from V03
// Maintains identical likelihood computation and random effects structure

functions {
  // transform the rowIndex and colIndex to a linear index
  int index_pos(int rowIndex, int colIndex, int nActors) {
    // index transformation works for square matrices (one-mode)
    int rowLower = rowIndex <= colIndex ? rowIndex : colIndex;
    int colHigher = rowIndex <= colIndex ? colIndex : rowIndex;

    return((nActors * (rowLower - 1)) + colHigher - choose(rowLower + 1, 2));
  }
}

data {
  int N; // number of events * actors * (actors - 1)
  int E; // number of events
  int MxAct; // maximum number of actors in an event

  int P; // number of covariates
  int G; // number of groups

  matrix[N, P] X_raw; // Raw (unstandardized) data matrix with covariates, include RE fix part
  vector[N] Z_raw; // Raw (unstandardized) Random Effect design vector
  array[N] int<lower = 1, upper = G> group; // Group assignments

  // number of actors in an event
  array[E] int<lower = 1, upper = MxAct> nA;

  // row col observed event
  array[E] int<lower = 1, upper = MxAct> rowOE;
  array[E] int<lower = 1, upper = MxAct> colOE;
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

  // Starting row index in X for each event
  array[total_actors] int event_start_X;

  // Pre-computed observed event positions
  array[E] int observed_pos;

  // Calculate maximum dyads across all events
  int max_dyads = choose(MxAct, 2);

  // Pre-compute exact dyad counts for each event (avoids repeated choose() calls)
  array[E] int event_nDyads;

  // Pre-compute all indices
  {
    int flat_idx = 1;
    int flat_row = 1;

    for (t in 1:E) {
      // Pre-compute observed positions
      observed_pos[t] = index_pos(rowOE[t], colOE[t], nA[t]);
      // dyads in event
      event_nDyads[t] = choose(nA[t], 2);

      for (rowChoice in 1:nA[t]) {
        // Pre-compute all indices for this choice set
        event_start_X[flat_row] = flat_idx;
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
    }
    print("flat_row = ", flat_row);
    print("flat_idx = ", flat_idx);
  }
}

parameters {
  vector[P] beta_std;  // Standardized coefficients
  real<lower=0> sigma; // variance RE
  vector[G] gamma_raw;  // group Random Effect
}

transformed parameters {
  // Back-transform coefficients to original scale
  vector[P] beta = beta_std ./ X_sds;

  // Back-transform random effects to original scale (corrected: division)
  vector[G] gamma = (sigma * gamma_raw) / Z_sd;
}

model {
  // Enhanced priors
  target += std_normal_lpdf(beta_std);  // Standard normal prior on standardized coefficients
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);

  // Compute utility with random effects using standardized predictors and coefficients
  vector[N] xb = X * beta_std + Z .* (sigma * gamma_raw)[group];

  // Memory-optimized likelihood computation
  vector[max_dyads] logLikDyad = rep_vector(0.0, max_dyads);
  vector[MxAct - 1] lprobs = rep_vector(0.0, MxAct - 1);
  int flat_row = 1;

  // Likelihood computation using pre-computed indices
  for (t in 1:E) {
    // Reset only the portion we actually need
    int current_nDyads = event_nDyads[t];
    array[nA[t] - 1] int slice_event = linspaced_int_array(nA[t] - 1, 1, nA[t] - 1);
    logLikDyad[1:current_nDyads] = rep_vector(0.0, current_nDyads);

    for (row in 1:nA[t]) {
      // probability of i choosing j
      lprobs[slice_event] = log_softmax(xb[X_indices[flat_row, slice_event]]);
      logLikDyad[dyad_indices[flat_row, slice_event]] += lprobs[slice_event];
      flat_row += 1;
    }

    // Add event likelihood using pre-computed observed position
    target += logLikDyad[observed_pos[t]] - log_sum_exp(logLikDyad[1:current_nDyads]);
  }
}