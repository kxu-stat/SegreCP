fast_gs_state <- function(observed, components) {
  caches <- lapply(components, function(z)
    SegreCP:::.rank1_vector(z$lambda, z$U))
  list(components = components, caches = caches,
    residual = array(Reduce(`+`, caches) - as.numeric(observed), dim(observed)))
}

fast_gs_reference <- function(state, order_i, alpha) {
  components <- state$components
  for (i in order_i) {
    proposal <- SegreCP:::rgd_step_fused_cpp(state$residual,
      components[[i]]$U, components[[i]]$lambda, alpha, 0L)
    next_component <- SegreCP:::.canonicalize_component(
      list(lambda = proposal$lambda, U = proposal$U))
    SegreCP:::update_rank1_cache_inplace_cpp(state$residual,
      state$caches[[i]], next_component$lambda, next_component$U)
    components[[i]] <- next_component
  }
  state$components <- components
  state$objective <- 0.5 * sum(state$residual * state$residual)
  state
}

fast_gs_component_arrays <- function(components) {
  lapply(components, function(z) rank1_array(z$lambda, z$U))
}

test_that("batched contractions preserve sequential GS trajectories", {
  set.seed(431)
  cases <- list(
    list(dims = c(4L, 5L), rank = 2L),
    list(dims = c(4L, 3L, 5L), rank = 3L),
    list(dims = c(5L, 4L, 3L, 4L), rank = 3L),
    list(dims = c(4L, 3L, 5L, 4L), rank = 2L),
    list(dims = c(4L, 5L, 4L, 5L, 4L), rank = 4L),
    list(dims = c(3L, 4L, 3L, 4L), rank = 1L)
  )
  for (case in cases) {
    dims <- case$dims; rank <- case$rank
    components <- lapply(seq_len(rank), function(i)
      list(lambda = 0.8 + i / rank, U = lapply(dims, unit_vector)))
    observed <- reconstruct_cp(components) + array(rnorm(prod(dims), sd = 0.12), dims)
    original_components <- serialize(components, NULL)
    for (alpha in c(0.1, 0.5, 1)) {
      reference <- fast_gs_state(observed, components)
      candidate <- fast_gs_state(observed, components)
      for (iteration in 1:5) {
        order_i <- if (iteration %% 2L) rev(seq_len(rank)) else
          order(vapply(reference$components, `[[`, numeric(1), "lambda"),
                decreasing = TRUE)
        reference <- fast_gs_reference(reference, order_i, alpha)
        before <- serialize(candidate$components, NULL)
        result <- SegreCP:::rgd_sweep_gs_batch_cpp(candidate$residual,
          candidate$caches, candidate$components, as.integer(order_i), alpha)
        expect_identical(serialize(candidate$components, NULL), before)
        candidate$components <- result$components
        expect_equal(fast_gs_component_arrays(candidate$components),
          fast_gs_component_arrays(reference$components), tolerance = 3e-11)
        expect_equal(candidate$residual, reference$residual, tolerance = 3e-11)
        expect_equal(result$objective, reference$objective, tolerance = 3e-12)
        expect_equal(result$objective,
          0.5 * sum(candidate$residual * candidate$residual), tolerance = 3e-14)
        for (i in seq_len(rank))
          expect_equal(candidate$caches[[i]], as.numeric(rank1_array(
            candidate$components[[i]]$lambda, candidate$components[[i]]$U)),
            tolerance = 2e-13)
        expect_equal(candidate$residual,
          array(Reduce(`+`, candidate$caches) - as.numeric(observed), dims),
          tolerance = 3e-12)
      }
    }
    expect_identical(serialize(components, NULL), original_components)
  }
})

test_that("GS contraction corrections handle zero cross-mode inner products", {
  set.seed(432)
  truth <- orthogonal_cp_tensor(c(4L, 5L, 3L, 4L), c(2, 1.5, 1))
  components <- truth$components
  components[[2L]]$U[[1L]] <- -components[[2L]]$U[[1L]]
  observed <- truth$array + array(rnorm(length(truth$array), sd = 0.03), dim(truth$array))
  for (alpha in c(0.5, 1)) {
    reference <- fast_gs_reference(fast_gs_state(observed, components), c(3L, 1L, 2L), alpha)
    candidate <- fast_gs_state(observed, components)
    result <- SegreCP:::rgd_sweep_gs_batch_cpp(candidate$residual,
      candidate$caches, components, c(3L, 1L, 2L), alpha)
    expect_false(result$sequential_fallback)
    expect_equal(fast_gs_component_arrays(result$components),
      fast_gs_component_arrays(reference$components), tolerance = 3e-12)
    expect_equal(candidate$residual, reference$residual, tolerance = 3e-12)
    expect_true(all(vapply(result$components, function(z) z$lambda > 0, logical(1))))
  }
})

test_that("tiny residuals use all original sequential updates without stopping", {
  set.seed(433)
  dims <- c(5L, 4L, 3L, 4L)
  components <- lapply(c(3, 2, 1), function(weight)
    list(lambda = weight, U = lapply(dims, unit_vector)))
  observed <- reconstruct_cp(components) + array(rnorm(prod(dims), sd = 1e-13), dims)
  for (alpha in c(0.5, 1)) {
    reference <- fast_gs_reference(fast_gs_state(observed, components), c(2L, 3L, 1L), alpha)
    candidate <- fast_gs_state(observed, components)
    before <- serialize(components, NULL)
    result <- SegreCP:::rgd_sweep_gs_batch_cpp(candidate$residual,
      candidate$caches, components, c(2L, 3L, 1L), alpha)
    expect_true(result$sequential_fallback)
    expect_identical(serialize(components, NULL), before)
    expect_equal(result$components, reference$components, tolerance = 0)
    expect_equal(candidate$residual, reference$residual, tolerance = 0)
    expect_equal(candidate$caches, reference$caches, tolerance = 0)
    expect_equal(result$objective, reference$objective, tolerance = 2e-15)
  }
})

test_that("batched GS validates the fixed component permutation", {
  truth <- orthogonal_cp_tensor(c(3L, 4L, 3L), c(2, 1))
  state <- fast_gs_state(truth$array + 0.01, truth$components)
  expect_error(SegreCP:::rgd_sweep_gs_batch_cpp(state$residual,
    state$caches, state$components, c(1L, 1L), 0.5), "permutation")
  expect_error(SegreCP:::rgd_sweep_gs_batch_cpp(state$residual,
    state$caches, state$components, c(2L, 1L), 0), "positive")
  expect_error(SegreCP:::rgd_sweep_gs_batch_cpp(state$residual,
    state$caches, state$components, c(NA_integer_, 1L), 0.5), "permutation")
})

test_that("fused initial state preserves generic rank-one multiplication order", {
  set.seed(434)
  cases <- list(
    list(dims = c(3L, 4L), rank = 1L),
    list(dims = c(1L, 4L, 3L), rank = 2L),
    list(dims = c(4L, 3L, 5L, 4L), rank = 3L),
    list(dims = c(3L, 4L, 3L, 4L), rank = 5L),
    list(dims = c(2L, 3L, 2L, 3L, 2L, 3L), rank = 7L)
  )
  for (case in cases) {
    dims <- case$dims; rank <- case$rank
    components <- lapply(seq_len(rank), function(i)
      list(lambda = 1 + i / rank, U = lapply(dims, unit_vector)))
    clean <- array(Reduce(`+`, lapply(components, function(z)
      SegreCP:::.rank1_vector(z$lambda, z$U))), dims)
    noise <- array(rnorm(prod(dims)), dims)
    for (sigma in c(0, 1e-13, 0.2)) {
      observed <- clean + sigma * noise
      original <- serialize(list(observed, components), NULL)
      expected <- fast_gs_state(observed, components)
      fused <- SegreCP:::initialize_decomposition_state_cpp(observed, components)
      expect_identical(fused$caches, expected$caches)
      expect_identical(fused$residual, expected$residual)
      expect_equal(fused$objective,
        0.5 * sum(expected$residual * expected$residual), tolerance = 3e-14)
      expect_identical(serialize(list(observed, components), NULL), original)
      expect_equal(fused$residual,
        array(Reduce(`+`, fused$caches) - as.numeric(observed), dims), tolerance = 0)
    }
  }
})

test_that("previous objectives preserve generic GS sweeps and early fallback", {
  set.seed(435)
  cases <- list(
    list(dims = c(3L, 4L), rank = 1L),
    list(dims = c(1L, 4L, 3L), rank = 2L),
    list(dims = c(3L, 4L, 3L, 4L), rank = 3L),
    list(dims = c(3L, 4L, 3L, 4L), rank = 5L),
    list(dims = c(2L, 3L, 2L, 3L, 2L, 3L), rank = 7L)
  )
  for (case in cases) {
    dims <- case$dims; rank <- case$rank
    components <- lapply(seq_len(rank), function(i)
      list(lambda = 1 + i / rank, U = lapply(dims, unit_vector)))
    clean <- array(Reduce(`+`, lapply(components, function(z)
      SegreCP:::.rank1_vector(z$lambda, z$U))), dims)
    noise <- array(rnorm(prod(dims)), dims)
    for (sigma in c(0, 1e-13, 0.2)) for (alpha in c(0.2, 0.5, 1)) {
      observed <- clean + sigma * noise
      reference <- fast_gs_state(observed, components)
      scanned <- fast_gs_state(observed, components)
      supplied <- SegreCP:::initialize_decomposition_state_cpp(observed, components)
      supplied$components <- components
      for (iteration in 1:3) {
        order_i <- as.integer(if (iteration %% 2L) rev(seq_len(rank)) else
          order(vapply(reference$components, `[[`, numeric(1), "lambda"),
                decreasing = TRUE))
        reference <- fast_gs_reference(reference, order_i, alpha)
        original <- serialize(supplied$components, NULL)
        a <- SegreCP:::rgd_sweep_gs_batch_cpp(scanned$residual,
          scanned$caches, scanned$components, order_i, alpha)
        b <- SegreCP:::rgd_sweep_gs_batch_cpp(supplied$residual,
          supplied$caches, supplied$components, order_i, alpha,
          supplied$objective)
        expect_identical(serialize(supplied$components, NULL), original)
        expect_equal(a$components, b$components, tolerance = 0)
        expect_equal(scanned$residual, supplied$residual, tolerance = 0)
        expect_equal(scanned$caches, supplied$caches, tolerance = 0)
        expect_identical(a$sequential_fallback, b$sequential_fallback)
        expect_equal(fast_gs_component_arrays(b$components),
          fast_gs_component_arrays(reference$components), tolerance = 3e-11)
        expect_equal(supplied$residual, reference$residual, tolerance = 3e-11)
        expect_equal(b$objective,
          0.5 * sum(supplied$residual * supplied$residual), tolerance = 3e-14)
        expect_equal(supplied$residual,
          array(Reduce(`+`, supplied$caches) - as.numeric(observed), dims),
          tolerance = 3e-12)
        if (sigma <= 1e-13) expect_true(b$sequential_fallback)
        scanned$components <- a$components
        supplied$components <- b$components
        supplied$objective <- b$objective
      }
    }
  }
})

test_that("previous objective and fused initializer validate their inputs", {
  truth <- orthogonal_cp_tensor(c(3L, 4L, 3L), c(2, 1))
  state <- fast_gs_state(truth$array + 0.01, truth$components)
  for (previous in list(-1, c(0, 1), NA_real_, Inf))
    expect_error(SegreCP:::rgd_sweep_gs_batch_cpp(state$residual,
      state$caches, state$components, c(2L, 1L), 0.5, previous),
      "finite nonnegative scalar")
  wrong <- truth$components
  wrong[[1L]]$U[[2L]] <- wrong[[1L]]$U[[2L]][-1L]
  expect_error(SegreCP:::initialize_decomposition_state_cpp(truth$array, wrong),
    "dimensions")
  wrong <- truth$components
  wrong[[1L]]$lambda <- 0
  expect_error(SegreCP:::initialize_decomposition_state_cpp(truth$array, wrong),
    "positive finite")
})

test_that("stored objective rounding cannot change the fallback at its boundary", {
  truth <- orthogonal_cp_tensor(c(3L, 4L, 3L), c(2, 1))
  boundary <- 64 * sqrt(.Machine$double.eps) * 3
  for (offset in c(-4, 0, 4) * .Machine$double.eps) {
    for (rounding in c(-4, 0, 4) * .Machine$double.eps) {
      direct <- fast_gs_state(truth$array, truth$components)
      supplied <- fast_gs_state(truth$array, truth$components)
      direct$residual <- array(0, dim(truth$array))
      supplied$residual <- array(0, dim(truth$array))
      direct$residual[[1L]] <- boundary * (1 + offset)
      supplied$residual[[1L]] <- boundary * (1 + offset)
      previous <- 0.5 * sum(supplied$residual^2) * (1 + rounding)
      a <- SegreCP:::rgd_sweep_gs_batch_cpp(direct$residual, direct$caches,
        direct$components, c(2L, 1L), 0.5)
      b <- SegreCP:::rgd_sweep_gs_batch_cpp(supplied$residual, supplied$caches,
        supplied$components, c(2L, 1L), 0.5, previous)
      expect_identical(a$sequential_fallback, b$sequential_fallback)
      expect_equal(a$components, b$components, tolerance = 0)
      expect_equal(direct$residual, supplied$residual, tolerance = 0)
      expect_equal(direct$caches, supplied$caches, tolerance = 0)
      expect_equal(a$objective, b$objective, tolerance = 0)
    }
  }
})
