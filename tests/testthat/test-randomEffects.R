test_that("warnings and stops", {
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      model = "dynam",
      data = testData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      model = "REM",
      sub_model = "choose",
      data = testData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      model = "REM",
      sub_model = "choice_coordination",
      data = testData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1, recip ~ 1),
      fixed_effects = depNetwork ~ trans,
      data = testData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip * trans,
      data = testData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ ego(actorsEx$attr1) * outdeg),
      fixed_effects = depNetwork ~ recip + trans,
      data = testData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      support_constraint = ~ tie(networkExog) + recip(networkExog),
      data = testData
    )
  )
})
test_that("choice empty model RE", {
  res <- make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      data = testData
  )
  expect_type(res, "list")
  expect_length(res, 5)
  expect_equal(
    res$senders_ix,
    data.frame(
      label = sprintf("Actor %d", 1:5),
      index = 1:5
    )
  )
  expect_equal(
    res$names_effects,
    c("recip" = "recip_networkState", "trans" = "trans_networkState",
      "inertia" = "inertia_networkState")
  )
  expect_equal(res$data_stan$T_choice, 12)
  expect_equal(res$data_stan$P_choice, 3)
  expect_equal(res$data_stan$A, nrow(actorsEx))
})
test_that("choice RE with expl effects", {
  res <- make_data_re(
    random_effects = list(inertia ~ outdeg),
    fixed_effects = depNetwork ~ recip + trans,
    data = testData
  )
  expect_type(res, "list")
  expect_length(res, 5)
  expect_equal(
    res$senders_ix,
    data.frame(
      label = sprintf("Actor %d", 1:5),
      index = 1:5
    )
  )
  # the cross-level term is now the native interaction (reMain:predictor); the
  # standalone level-2 operand column is dropped from the design matrix
  expect_equal(
    res$names_effects,
    c("recip" = "recip_networkState", "trans" = "trans_networkState",
      "inertia" = "inertia_networkState",
      'inertia:outdeg(networkState, type = "ego")' =
        "inertia_networkState:outdeg_networkState_ego")
  )
  expect_equal(res$data_stan$T_choice, 12)
  expect_equal(res$data_stan$P_choice, 4)
  expect_equal(res$data_stan$A, nrow(actorsEx))
})
test_that("save code", {
  skip_if_not_installed("cmdstanr")
  expect_error(make_model_code(list()))

  res <- make_data_re(
    random_effects = list(inertia ~ 1),
    fixed_effects = depNetwork ~ recip + trans,
    data = testData
  )
  out_code <- make_model_code(res)
  expect_length(out_code, 1)
  expect_type(out_code, "character")
})
