# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`goldfish.latent` is an R package that extends the `goldfish` package with Bayesian relational event models featuring latent variables. It focuses on Dynamic Network Actor Models (DyNAMs) and Relational Event Models (REMs) with random effects and Hidden Markov Models (HMMs), using Stan for Hamiltonian Monte Carlo inference.

## Architecture

### Core Components

- **Data Transformation**: Functions in `R/transform-data.R` and `R/randomEffects-actor.R` convert goldfish data objects into Stan-compatible formats
- **Stan Models**: Located in `inst/stan/` directory containing various model implementations:
  - `DyNAM_*.stan`: Basic DyNAM models for choice and rate sub-models
  - `DNRE_*.stan`: DyNAM with random effects variants
  - `DNHMM_*.stan`: DyNAM with Hidden Markov Model components
  - `DNCHMM_*.stan`: DyNAM with choice-specific HMM
- **Model Generation**: `R/stan_generator.R` contains logic for generating Stan code
- **Utility Functions**: `R/utility.R` provides helper functions for data manipulation and model fitting
- **HMM Processing**: `R/HMM.R` handles Hidden Markov Model-specific functionality

### Key Functions

- `make_data_re()`: Transforms goldfish data into Stan format for random effects models
- `make_data_hmm()`: Prepares data for HMM models
- `make_model_code()`: Generates Stan model code based on data structure
- `compute_log_likelihood()`: Computes marginal/conditional log-likelihood using MCMC samples
- `gather_groups()` and `gather_groups_to_json()`: Data aggregation and export functions

## Development Commands

### R Package Development
```r
# Install development version
remotes::install_github("snlab-ch/goldfish.latent", build_vignettes = TRUE)

# Build package
devtools::build()

# Check package
devtools::check()

# Run tests
devtools::test()
# Or using testthat directly
testthat::test_dir("tests/testthat")

# Generate documentation
devtools::document()

# Show code that doesn't follow the package style
lintr::lint_package()

# Test code coverage
covr::codecov()
covr::report(x = covr::package_coverage(type = "tests"))

# test website
pkgdown::build_site()
```

### Stan Model Development
- Stan models are pre-written in `inst/stan/` and selected/modified by R functions
- No direct Stan compilation needed as models are compiled via `cmdstanr` at runtime
- Models support map-reduce parallelization (files ending with `_map_reduce.stan`)

### Testing
- Test files located in `tests/testthat/`
- Currently focused on random effects functionality (`test-randomEffects.R`)
- Tests use `testthat` framework (version 3.0.0+)

## Dependencies

### Required System Dependencies
- **CmdStan >= 2.29.2**: Stan's command-line interface for model compilation and sampling
- **R >= 4.1.0**: Base R installation

### Key R Package Dependencies
- `goldfish (>= 1.6.10)`: Base package for relational event modeling
- `cmdstanr (>= 0.8.0)`: R interface to CmdStan
- `posterior (>= 1.5.0)`: Working with MCMC output
- `bayesplot (>= 1.10.0)`: Plotting MCMC results
- `loo`: Leave-one-out cross-validation and model comparison

## Model Types Supported

1. **DyNAM with Random Effects**: Actor-level random effects in choice and rate models
2. **DyNAM with HMM**: Hidden Markov Models for capturing latent states
3. **Combined Models**: Models with both random effects and HMM components
4. **Parallel Computation**: Map-reduce implementations for large datasets

## Data Flow

1. Start with `goldfish` data objects (networks, events, actors)
2. Use `make_data_re()` or `make_data_hmm()` to transform for Stan
3. Generate Stan model code with `make_model_code()`
4. Fit model using `cmdstanr` functions
5. Post-process results with package utilities like `compute_log_likelihood()`

## Important Notes

- Models are designed to work with the `goldfish` package's data structures
- All Bayesian inference is performed through Stan/CmdStan
- The package focuses on network event data with temporal dynamics
- Supports both choice (who to interact with) and rate (when to interact) modeling


## Spec-driven development (OpenSpec)

This repo uses OpenSpec (`openspec/`); see the user-level `~/.claude/CLAUDE.md`
for the general workflow. Project disciplines (commit-per-task; run
`devtools::document()` inline on roxygen/export/signature changes; bump
`DESCRIPTION` + `NEWS.md` at milestones) live authoritatively in
`openspec/config.yaml`.

## Personal defaults

R style, performance guidelines, the code-comments policy, and the `pbcopy` copy
command live in the user-level `~/.claude/CLAUDE.md` and apply here automatically.