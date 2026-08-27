mixed_test_data <- function() {
  data.frame(
    continuous = seq(-1, 1, length.out = 8),
    nominal = factor(rep(c("a", "b"), 4)),
    ordinal = ordered(
      rep(c("low", "middle", "high", "middle"), 2),
      levels = c("low", "middle", "high")
    )
  )
}

test_that("compiled CPU sampling and reconstruction work for both metrics", {
  data <- mixed_test_data()
  samplers <- list(dkps = mcmc.dkps, dkss = mcmc.dkss)
  distance_wrappers <- list(dkps = dkps, dkss = dkss)
  fits <- list()

  for (metric in names(samplers)) {
    set.seed(404)
    fit <- samplers[[metric]](
      data, chains = 2, warmup = 6, draws = 6,
      backend = "cpu", progress = FALSE
    )
    fits[[metric]] <- fit

    expect_true(
      inherits(fit, "kdml_mcmc") &&
        identical(fit$backend, "cpu") &&
        identical(dim(fit$theta), c(6L, 3L, 2L)) &&
        all(is.finite(fit$bandwidth)) &&
        all(is.finite(fit$score)) &&
        all(is.finite(fit$log_target)) &&
        posterior::is_draws_array(fit$draws)
    )

    reconstructed <- distance_wrappers[[metric]](fit)
    expect_true(
      all(is.finite(reconstructed$distances)) &&
        isTRUE(all.equal(
          reconstructed$distances,
          t(reconstructed$distances),
          tolerance = 1e-12
      )) &&
      all(abs(diag(reconstructed$distances)) < 1e-12) &&
      max(reconstructed$distances) <= 1 + 1e-12 &&
      identical(reconstructed$state, fit$map[c("draw", "chain")])
    )
  }

  diagnostics <- suppressWarnings(kdml.diagnostics(fits$dkps))
  summary <- suppressWarnings(summary(fits$dkps, min_ess_per_chain = 1))
  expect_true(
    nrow(diagnostics) == ncol(data) &&
      nrow(summary$acceptance) == 2 * ncol(data)
  )
})

test_that("compiled fixed distances and CUDA inspection keep their contracts", {
  data <- mixed_test_data()
  fixed <- dkps(
    data,
    bw = c(continuous = 0.6, nominal = 0.2, ordinal = 0.3),
    oFUN = "o_liracine",
    stan = FALSE
  )
  expect_true(
    all(is.finite(fixed$distances)) &&
      isTRUE(all.equal(fixed$distances, t(fixed$distances))) &&
      all(abs(diag(fixed$distances)) < 1e-12)
  )

  info <- cuda_info()
  expect_true(
    is.logical(info$available) && length(info$available) == 1L &&
      !is.na(info$available) &&
      identical(
        names(info),
        c(
          "compiled", "available", "compile_version",
          "compile_version_string", "runtime_version", "device_count",
          "devices"
        )
      )
  )
  expect_false("cuda_available" %in% getNamespaceExports("kdml"))
})
