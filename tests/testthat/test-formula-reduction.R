# Equivalence oracle for the D11 formula reduction (native interaction terms).
# `_oracle_re_choice.rds` was captured from the pre-reduction `model.matrix()`
# path; after the reduction `make_data_re()` must reproduce it for the two clean
# cases (empty_re, cross_lvl). support_constraint is temporarily disabled, so
# the constraint path is not exercised here.

oracle <- readRDS(test_path("_oracle_re_choice.rds"))

# sort design-matrix columns by name (order-insensitive) and drop row names
# (the old model.matrix path carried expanded_df row indices; Stan ignores them)
sorted_cols <- function(m) {
  m <- m[, order(colnames(m)), drop = FALSE]
  rownames(m) <- NULL
  m
}

build_choice <- function(re, sc = NULL) {
  suppressWarnings(
    make_data_re(
      random_effects = re,
      fixed_effects = depNetwork ~ recip + trans,
      support_constraint = sc,
      data = testData
    )
  )
}

test_that("empty RE choice matches the pre-reduction oracle", {
  o <- oracle$empty_re
  res <- build_choice(list(inertia ~ 1))
  ds <- res$data_stan
  expect_equal(ds$P_choice, o$P)
  expect_equal(ds$Q_choice, o$Q)
  expect_equal(ds$N_choice, o$N)
  expect_equal(ds$T_choice, o$T)
  expect_equal(sorted_cols(ds$X_choice), sorted_cols(o$X))
  expect_equal(unname(ds$Z_choice), unname(o$Z))
  expect_equal(ds$chose_choice, o$chose)
  expect_equal(ds$start_choice, o$start)
  expect_equal(ds$end_choice, o$end)
  expect_equal(ds$sender, o$sender)
})

test_that("cross-level RE choice (interaction) matches the pre-reduction oracle", {
  o <- oracle$cross_lvl
  res <- build_choice(list(inertia ~ outdeg))
  ds <- res$data_stan
  expect_equal(ds$P_choice, o$P)
  expect_equal(ds$Q_choice, o$Q)
  expect_equal(ds$N_choice, o$N)
  expect_equal(ds$T_choice, o$T)
  # X carries the cross-level interaction column and drops the standalone
  # level-2 operand; equal to the old model.matrix product up to column order
  expect_equal(sorted_cols(ds$X_choice), sorted_cols(o$X))
  expect_equal(unname(ds$Z_choice), unname(o$Z))
  expect_equal(ds$chose_choice, o$chose)
  expect_equal(ds$start_choice, o$start)
  expect_equal(ds$end_choice, o$end)
})

test_that("support_constraint is rejected while temporarily disabled", {
  expect_error(
    build_choice(list(inertia ~ 1), sc = ~ tie(networkExog)),
    "support_constraint"
  )
})
