rebuttal_extdata <- function(...) {
  parts <- c(...)
  installed <- do.call(
    system.file,
    c(as.list(c("extdata", parts)), list(package = "SegreCP"))
  )
  candidates <- c(
    installed,
    do.call(testthat::test_path,
            as.list(c("..", "..", "inst", "extdata", parts))),
    do.call(file.path, as.list(c("inst", "extdata", parts))),
    do.call(file.path,
            as.list(c("software", "SegreCP", "inst", "extdata", parts)))
  )
  candidates <- unique(candidates[nzchar(candidates)])
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("Cannot locate the bundled rebuttal extdata: ",
         paste(parts, collapse = "/"), call. = FALSE)
  }
  normalizePath(found[[1L]], mustWork = TRUE)
}

read_rebuttal_csv <- function(...) {
  utils::read.csv(rebuttal_extdata(...), stringsAsFactors = FALSE,
                  check.names = FALSE)
}

expect_frozen_means <- function(observed, expected, tolerance = 5e-12) {
  expect_setequal(names(observed), names(expected))
  observed <- observed[names(expected)]
  expect_true(all(is.finite(observed)))
  expect_lt(max(abs(unname(observed) - unname(expected))), tolerance)
}

test_that("the frozen rebuttal grid is complete and paired", {
  raw <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "final_seven_method_grid_raw.csv"
  )
  checks <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "final_validation_checks.csv"
  )

  methods <- c("CP-ALS", "RRR", "Tucker-RGD", "HOPLS", "TT-IHT",
               "RGD", "RGN")
  expect_equal(nrow(raw), 4L * 3L * 20L * 7L)
  expect_setequal(unique(raw$display_method), methods)
  expect_equal(sort(unique(raw$coherence)), c(0, 0.2, 0.5, 0.75))
  expect_equal(sort(unique(raw$noise_sd)), c(0, 0.5, 1))
  expect_equal(sort(unique(raw$replication)), seq_len(20L))

  paired_counts <- aggregate(
    display_method ~ coherence + noise_sd + replication,
    raw, function(x) length(unique(x))
  )
  paired_seeds <- aggregate(
    seed ~ coherence + noise_sd + replication,
    raw, function(x) length(unique(x))
  )
  expect_true(all(paired_counts$display_method == 7L))
  expect_true(all(paired_seeds$seed == 1L))
  expect_equal(nrow(checks), 10L)
  expect_true(all(checks$passed))

  manifest <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "package_subset_manifest.csv"
  )
  paths <- vapply(
    manifest$file,
    function(file) rebuttal_extdata("rebuttal_tensor_regression", file),
    character(1L)
  )
  observed_hashes <- vapply(
    paths, digest::digest, character(1L), algo = "sha256",
    file = TRUE, serialize = FALSE
  )
  expect_true(all(file.exists(paths)))
  expect_identical(unname(observed_hashes), manifest$sha256)
})

test_that("Panel A checkpoint means reproduce the rebuttal table", {
  checkpoints <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "checkpoint_metrics_raw.csv"
  )
  observed <- aggregate(
    coefficient_error ~ method + iteration, checkpoints, mean
  )
  observed <- setNames(
    observed$coefficient_error,
    paste(observed$method, observed$iteration, sep = "@")
  )
  expected <- c(
    `CP-ALS@3` = 0.000245866895920443,
    `RGD@3` = 0.00244467370308301,
    `RGN@3` = 0.000192464431506970,
    `RRR@3` = 0.0950529703928841,
    `TT-IHT@3` = 0.112198549637987,
    `Tucker-RGD@3` = 0.0478491301149665,
    `CP-ALS@5` = 6.75924338461747e-06,
    `RGD@5` = 0.000171810661613999,
    `RGN@5` = 4.40869509456659e-06,
    `RRR@5` = 0.0397959924148064,
    `TT-IHT@5` = 0.0496940350719988,
    `Tucker-RGD@5` = 0.0295288480594093,
    `CP-ALS@10` = 1.07195302102063e-08,
    `RGD@10` = 2.86734364206143e-07,
    `RGN@10` = 1.15552159343935e-08,
    `RRR@10` = 0.0276398196594550,
    `TT-IHT@10` = 0.00190910205624562,
    `Tucker-RGD@10` = 0.00119914609781764
  )
  expect_frozen_means(observed, expected)
})

test_that("Panel B high-coherence means reproduce the rebuttal table", {
  raw <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "final_seven_method_grid_raw.csv"
  )
  observed <- aggregate(
    clean_test_rmse ~ display_method + noise_sd,
    raw[raw$coherence == 0.75, ], mean
  )
  observed <- setNames(
    observed$clean_test_rmse,
    paste(observed$display_method, observed$noise_sd, sep = "@")
  )
  expected <- c(
    `CP-ALS@0` = 1.16161401903181e-05,
    `HOPLS@0` = 0.00229465428799533,
    `RGD@0` = 1.06743020276140e-05,
    `RGN@0` = 1.25015809831514e-05,
    `RRR@0` = 0.000592022914515644,
    `TT-IHT@0` = 0.000113738620041145,
    `Tucker-RGD@0` = 6.63151154582903e-06,
    `CP-ALS@0.5` = 0.000859706567656363,
    `HOPLS@0.5` = 0.00237598378657836,
    `RGD@0.5` = 0.000851742354301327,
    `RGN@0.5` = 0.000857397214666973,
    `RRR@0.5` = 0.00101353744882348,
    `TT-IHT@0.5` = 0.00118218894017885,
    `Tucker-RGD@0.5` = 0.00119275346055023,
    `CP-ALS@1` = 0.00152245167321292,
    `HOPLS@1` = 0.00270811754430774,
    `RGD@1` = 0.00151147846508366,
    `RGN@1` = 0.00151811990941374,
    `RRR@1` = 0.00225790443926399,
    `TT-IHT@1` = 0.00266509584425674,
    `Tucker-RGD@1` = 0.00229482738971628
  )
  expect_frozen_means(observed, expected)
})

test_that("Panel C uses replication-balanced rebuttal means", {
  raw <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "final_seven_method_grid_raw.csv"
  )
  frozen <- read_rebuttal_csv(
    "rebuttal_tensor_regression", "final_seven_method_grid_summary.csv"
  )
  metrics <- c("coefficient_error", "clean_test_rmse", "test_rmse")
  by_replication <- aggregate(
    raw[metrics],
    list(display_method = raw$display_method,
         replication = raw$replication),
    mean
  )
  rebuilt <- aggregate(
    by_replication[metrics],
    list(display_method = by_replication$display_method),
    mean
  )
  rebuilt <- rebuilt[match(frozen$display_method, rebuilt$display_method), ]
  expect_lt(
    max(abs(as.matrix(rebuilt[metrics]) - as.matrix(frozen[metrics]))),
    5e-14
  )

  expected <- list(
    coefficient_error = c(
      `CP-ALS` = 0.0455094320192771,
      HOPLS = 0.126985181676494,
      RGD = 0.0448471790478060,
      RGN = 0.0449939495233415,
      RRR = 0.0802751614575916,
      `TT-IHT` = 0.0702875554239684,
      `Tucker-RGD` = 0.0685878439130270
    ),
    clean_test_rmse = c(
      `CP-ALS` = 0.000811712762062883,
      HOPLS = 0.00226047111409431,
      RGD = 0.000798588183230607,
      RGN = 0.000801613380475244,
      RRR = 0.00142895026182980,
      `TT-IHT` = 0.00125177170235762,
      `Tucker-RGD` = 0.00122271673151437
    ),
    test_rmse = c(
      `CP-ALS` = 0.0158408480682129,
      HOPLS = 0.0165950926443466,
      RGD = 0.0158403749863720,
      RGN = 0.0158406849769246,
      RRR = 0.0160872316838175,
      `TT-IHT` = 0.0158842450827309,
      `Tucker-RGD` = 0.0159170718984018
    )
  )
  for (metric in names(expected)) {
    observed <- setNames(rebuilt[[metric]], rebuilt$display_method)
    expect_frozen_means(observed, expected[[metric]])
  }
})

test_that("runtime validation uses stable tiers rather than exact clocks", {
  timing <- read_rebuttal_csv(
    "rebuttal_tensor_regression",
    "final_seven_method_serial_timing_summary.csv"
  )
  elapsed <- setNames(timing$mean_fit_sec, timing$display_method)

  expect_true(all(timing$runs == 20L))
  expect_lt(max(elapsed[c("RGD", "RGN")]), elapsed[["CP-ALS"]])
  expect_lt(elapsed[["CP-ALS"]], elapsed[["HOPLS"]])
  expect_lt(elapsed[["HOPLS"]], elapsed[["TT-IHT"]])
  expect_lt(elapsed[["TT-IHT"]],
            min(elapsed[c("RRR", "Tucker-RGD")]))
})

test_that("initializer validation bundles retain every frozen check", {
  decomposition <- read_rebuttal_csv(
    "rebuttal_init", "decomp_ref",
    "validation_checks.csv"
  )
  regression <- read_rebuttal_csv(
    "rebuttal_init", "reg_ref",
    "validation_checks.csv"
  )

  expect_equal(nrow(decomposition), 8L)
  expect_true(all(decomposition$passed))
  expect_equal(nrow(regression), 18L)
  expect_true(all(regression$pass))

  for (section in c("decomp_ref", "reg_ref", "decomp_pkg")) {
    manifest <- read_rebuttal_csv(
      "rebuttal_init", section, "package_subset_manifest.csv"
    )
    paths <- vapply(
      manifest$file,
      function(file) rebuttal_extdata("rebuttal_init", section, file),
      character(1L)
    )
    observed_hashes <- vapply(
      paths, digest::digest, character(1L), algo = "sha256",
      file = TRUE, serialize = FALSE
    )
    expect_true(all(file.exists(paths)))
    expect_identical(unname(observed_hashes), manifest$sha256)
  }
})

test_that("the package rerun satisfies objective and legacy parity contracts", {
  rerun <- read_rebuttal_csv(
    "rebuttal_init", "decomp_pkg",
    "trajectory_raw.csv"
  )
  rerun_checks <- read_rebuttal_csv(
    "rebuttal_init", "decomp_pkg",
    "validation_checks.csv"
  )
  parity <- read_rebuttal_csv(
    "rebuttal_init", "decomp_pkg",
    "legacy_reference_parity.csv"
  )

  expect_true(all(rerun_checks$passed))
  package_paths <- rerun[rerun$method %in% c("CP-ALS", "RGD", "RGN"), ]
  paths <- split(
    package_paths,
    interaction(package_paths$replication, package_paths$initialization,
                package_paths$method, drop = TRUE)
  )
  nonincreasing <- vapply(paths, function(path) {
    path <- path[order(path$iteration), ]
    all(is.finite(path$objective)) && all(diff(path$objective) <= 1e-8)
  }, logical(1L))
  expect_equal(length(nonincreasing), 20L * 2L * 3L)
  expect_true(all(nonincreasing))

  expect_true(parity$reference_available)
  expect_true(parity$method_ranking_preserved)
  expect_true(parity$initializer_contrast_preserved)
  expect_lte(parity$cpca_max_absolute_difference, parity$exact_tolerance)
  expect_lte(parity$cpals_max_absolute_difference, parity$exact_tolerance)
  expect_lte(parity$hosvd_initial_max_absolute_difference,
             parity$exact_tolerance)
  expect_lte(parity$max_absolute_mean_trajectory_difference,
             parity$mean_path_tolerance)
  expect_lte(parity$max_absolute_checkpoint_difference,
             parity$checkpoint_tolerance)
  expect_true(parity$passed)
})

test_that("a quick offline package rerun has non-increasing objectives", {
  set.seed(20260830)
  truth <- orthogonal_cp_tensor(c(4, 4, 4), c(3, 1.5))
  observed <- truth$array + array(
    stats::rnorm(length(truth$array), sd = 0.03), dim(truth$array)
  )
  initialization <- cp_initialize(
    observed, rank = 2, method = "cpca", split = 1
  )$comps

  fits <- list(
    RGD = segre_decompose(
      observed, rank = 2, method = "rgd", init = initialization,
      control = segre_control(
        profile = "practical", max_iter = 3, tolerance = 0,
        step_rule = "backtracking", component_order = "decreasing_weight"
      ),
      keep_history = TRUE
    ),
    RGN = segre_decompose(
      observed, rank = 2, method = "rgn", init = initialization,
      control = segre_control(
        profile = "manuscript", max_iter = 3, tolerance = 0
      ),
      keep_history = TRUE
    )
  )

  for (fit in fits) {
    if (!is.null(fit$objective)) {
      expect_length(fit$objective, fit$iterations + 1L)
      expect_true(all(is.finite(fit$objective)))
      expect_true(all(diff(fit$objective) <= 1e-10))
    }
    expect_length(fit$component_history, fit$iterations + 1L)
  }
})
