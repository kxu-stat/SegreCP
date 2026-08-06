test_that("the installed experiment entry point contains only authors' methods", {
  script <- system.file("scripts", "ours_methods.R", package = "SegreCP")
  expect_true(nzchar(script))
  expect_true(file.exists(script))

  script_environment <- new.env(parent = globalenv())
  sys.source(script, envir = script_environment)
  public_helpers <- ls(script_environment, all.names = FALSE)

  expect_setequal(
    public_helpers,
    c(
      "fit_ours_decomposition",
      "fit_ours_regression",
      "initialize_ours",
      "ours_convergence"
    )
  )
  expect_false(any(grepl(
    "als|orth|ico|rrr|tucker|iht|hopls",
    public_helpers,
    ignore.case = TRUE
  )))
})

test_that("comparison implementations are absent from the installed package", {
  expect_identical(
    system.file("reproducibility", package = "SegreCP"),
    ""
  )

  exports <- getNamespaceExports("SegreCP")
  expect_false(any(grepl(
    "als|orth|ico|rrr|tucker|iht|hopls",
    exports,
    ignore.case = TRUE
  )))
})
