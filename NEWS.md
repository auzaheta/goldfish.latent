# goldfish.latent (development version)

* Compatibility with the reworked `goldfish::gather_model_data()` output
  (goldfish >= 1.8.x). `make_data_re()` and `make_data_hmm()` again build Stan
  data against the current gather result: the obsolete `"DyNAMRE"` model remap
  is dropped, the 1-based `selected` vector is read as a vector (no longer a
  matrix), and `make_data_hmm()` forwards `control_preprocessing` and reads
  `has_intercept`.
* The R-side data contract is standardised on **snake_case**: return-list fields
  (`data_stan`, `senders_ix`, `names_effects`, `effect_description`) and the
  `sub_model` / `model` S3 attributes now use one spelling at every set- and
  read-site (fixing a silent set-camel/read-snake dispatch bug). Inner
  `data_stan` keys stay Stan-spelled to match the `.stan` `data{}` blocks.
* `make_data_re()` now rejects a `random_effects` list with more than one
  element (single-random-effect inference only, as `compute_log_likelihood()`
  already required).
* `make_data_hmm()` gains a first test suite (`tests/testthat/test-HMM.R`).

# goldfish.latent 0.0.1

* Random-effects for DyNAM-choice model functionality is added:
  * `CreateData()` generates the expected data for `Stan`.
  * `CreateModelCode()` copies the `Stan` code to a temporal folder.
    `cmdstanr` uses the path to the new copied code to compile the model and
    generate sampled from the posterior distribution.
  * `ComputeLogLikelihood()` creates an array with the log-likelihood of either
    the marginal or the conditional version. Output is ready to be used by `loo`
    package.

# goldfish.latent 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.

# goldfish.latent 0.0.2

* `scale` parameter in `CreateData()` allows to standardize variables and
  keep mean and standard deviation information.
* Internal function `RescaleCoef()` uses the standardization information to
  transform parameter values in the original variable scale.
* HMM-DyNAM includes a new Stan variant that allows to change temporal resolution,
  for example, the model considers changes the Hidden Markov happening after each
  day or other time window predefined.

