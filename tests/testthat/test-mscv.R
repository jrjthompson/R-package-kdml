with_stub_optimizer <- function(fun, optimizer) {
  local_environment <- new.env(parent = environment(fun))
  local_environment$optim <- optimizer
  environment(fun) <- local_environment
  fun
}

test_that("MSCV starts within bounds and rejects invalid candidates", {
  data <- data.frame(nominal = factor(rep(c("a", "b", "c"), each = 2)))

  optimizer <- function(par, fn, control) {
    expect_true(length(par) == 1L && par > 0 && par <= 2 / 3)
    expect_identical(fn(-0.1), -Inf)
    expect_identical(fn(0.8), -Inf)
    list(par = par, value = fn(par), convergence = 0)
  }

  for (fun in list(mscv.dkps, mscv.dkss)) {
    set.seed(101)
    result <- with_stub_optimizer(fun, optimizer)(
      data, nstart = 1, ukernel = "u_aitchisonaitken"
    )
    expect_true(is.finite(result$fn_value))
  }
})

test_that("MSCV retains the legacy dotted kernel alias", {
  data <- data.frame(nominal = factor(rep(c("a", "b", "c"), each = 2)))

  for (fun in list(mscv.dkps, mscv.dkss)) {
    set.seed(202)
    canonical <- suppressWarnings(fun(
      data, nstart = 1, ukernel = "u_aitchisonaitken"
    ))
    set.seed(202)
    legacy <- suppressWarnings(fun(
      data, nstart = 1, ukernel = "u.aitchisonaitken"
    ))
    expect_equal(legacy, canonical, tolerance = 1e-10)
  }
})
