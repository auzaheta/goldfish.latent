test_that("make_data_hmm rejects invalid k_states", {
  expect_error(
    make_data_hmm(
      rate_effects = depNetwork ~ 1 + indeg,
      choice_effects = NULL,
      k_states = 1,
      data = testData,
      progress = FALSE
    )
  )
})

test_that("make_data_hmm rate-only builds a DNHMM object", {
  res <- make_data_hmm(
    rate_effects = depNetwork ~ 1 + indeg + outdeg,
    choice_effects = NULL,
    k_states = 3,
    data = testData,
    progress = FALSE
  )
  expect_s3_class(res, "goldfish.latent.data")
  expect_equal(attr(res, "model"), "DNHMM")
  expect_equal(attr(res, "sub_model"), "rate")
  # top-level fields are snake_case (match make_data_re's contract)
  expect_named(
    res,
    c("data_stan", "names_effects", "effect_description", "scale_stats")
  )
  # inner data_stan keys stay Stan-spelled (documented domain-3 exception)
  expect_equal(res$data_stan$Trate, 12)
  expect_equal(res$data_stan$Prate, 3)
  expect_equal(res$data_stan$kS, 3)
  expect_length(res$data_stan$startRate, 12)
  expect_length(res$data_stan$endRate, 12)
  expect_equal(ncol(res$data_stan$Xrate), 3)
  expect_true(is.finite(res$data_stan$offsetInt))
})

test_that("make_data_hmm choice-only builds a DNHMM object", {
  res <- make_data_hmm(
    rate_effects = NULL,
    choice_effects = depNetwork ~ recip + trans,
    k_states = 3,
    data = testData,
    progress = FALSE
  )
  expect_s3_class(res, "goldfish.latent.data")
  expect_equal(attr(res, "model"), "DNHMM")
  expect_equal(attr(res, "sub_model"), "choice")
  expect_equal(res$data_stan$Tchoice, 12)
  expect_equal(res$data_stan$Pchoice, 2)
  expect_equal(ncol(res$data_stan$Xchoice), 2)
  expect_equal(res$data_stan$kS, 3)
})

test_that("make_data_hmm both sub-models carries rate and choice blocks", {
  res <- make_data_hmm(
    rate_effects = depNetwork ~ 1 + indeg + outdeg,
    choice_effects = depNetwork ~ recip + trans,
    k_states = 3,
    data = testData,
    progress = FALSE
  )
  expect_s3_class(res, "goldfish.latent.data")
  expect_equal(attr(res, "sub_model"), "both")
  # both Stan-spelled blocks present
  expect_equal(res$data_stan$Trate, 12)
  expect_equal(res$data_stan$Tchoice, 12)
  expect_false(is.null(res$data_stan$Xrate))
  expect_false(is.null(res$data_stan$Xchoice))
})

test_that("make_data_hmm applies a support_constraint to the choice set", {
  res_free <- make_data_hmm(
    rate_effects = NULL,
    choice_effects = depNetwork ~ recip + trans,
    k_states = 2,
    data = testData,
    progress = FALSE
  )
  res_cstr <- make_data_hmm(
    rate_effects = NULL,
    choice_effects = depNetwork ~ recip + trans,
    support_constraint = ~ tie(networkExog),
    k_states = 2,
    data = testData,
    progress = FALSE
  )
  # the constraint restricts the candidate set, so fewer stacked rows
  expect_lt(res_cstr$data_stan$Nchoice, res_free$data_stan$Nchoice)
  # the constraint column is dropped from the estimated effects
  expect_equal(res_cstr$data_stan$Pchoice, 2)
  expect_false(any(grepl("tie", res_cstr$names_effects)))
})
