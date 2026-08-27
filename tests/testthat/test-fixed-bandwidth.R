legacy_matrix <- function(fun, data, bandwidth, ...) {
  result <- fun(data, bw = bandwidth, stan = FALSE, ...)
  result[[if (identical(fun, kss)) "similarities" else "distances"]]
}

test_that("legacy type slices handle single and missing blocks", {
  nominal_only <- data.frame(group = factor(c("a", "b", "a", "c")))
  mixed <- data.frame(
    continuous = c(0, 0.5, 1.5, 2),
    ordinal = ordered(
      c("low", "middle", "high", "middle"),
      levels = c("low", "middle", "high")
    )
  )

  for (fun in list(dkps, dkss, kss)) {
    one_column <- legacy_matrix(fun, nominal_only, 0.3)
    expect_true(
      identical(dim(one_column), c(4L, 4L)) && all(is.finite(one_column))
    )

    combined <- legacy_matrix(fun, mixed, c(0.8, 0.3))
    continuous <- legacy_matrix(fun, mixed["continuous"], 0.8)
    ordinal <- legacy_matrix(fun, mixed["ordinal"], 0.3)
    expect_equal(combined, continuous + ordinal, tolerance = 1e-12)
  }
})

test_that("categorical kernels evaluate each feature independently", {
  nominal <- data.frame(
    three_levels = factor(c("a", "b", "c")),
    two_levels = factor(c("x", "y", "x"))
  )
  nominal_bandwidth <- c(0.3, 0.2)
  expected_similarity <- matrix(
    c(1.50, 0.35, 0.95,
      0.35, 1.50, 0.35,
      0.95, 0.35, 1.50),
    nrow = 3,
    byrow = TRUE
  )
  expected_distance <- outer(
    diag(expected_similarity), diag(expected_similarity), "+"
  ) - 2 * expected_similarity

  expect_equal(
    legacy_matrix(
      kss, nominal, nominal_bandwidth,
      uFUN = "u_aitchisonaitken"
    ),
    expected_similarity,
    tolerance = 1e-12
  )
  for (fun in list(dkps, dkss)) {
    expect_equal(
      legacy_matrix(
        fun, nominal, nominal_bandwidth,
        uFUN = "u_aitchisonaitken"
      ),
      expected_distance,
      tolerance = 1e-12
    )
  }

  ordinal <- data.frame(
    first = ordered(c("low", "low"), levels = c("low", "high")),
    second = ordered(c("one", "two"), levels = c("one", "two"))
  )
  expect_equal(
    legacy_matrix(kss, ordinal, c(0.2, 0.4)),
    matrix(c(1.4, 0.92, 0.92, 1.4), nrow = 2),
    tolerance = 1e-12
  )
  expect_equal(
    legacy_matrix(dkss, ordinal, c(0.2, 0.4)),
    matrix(c(0, 0.96, 0.96, 0), nrow = 2),
    tolerance = 1e-12
  )
})
