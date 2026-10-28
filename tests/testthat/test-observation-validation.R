test_that("invalid censoring codes cannot silently discard observations", {
  model <- dfr_exponential(1)
  for (evaluate in list(loglik(model), score(model), hess_loglik(model))) {
    for (status in list(2, -2, NA_real_, NaN, Inf, "1", factor("1"))) {
      df <- data.frame(t = c(1, 2), delta = c(status, status))
      expect_error(evaluate(df, par = 1), "event indicators")
    }
  }
})

test_that("likelihood methods validate observation times before evaluation", {
  model <- dfr_exponential(1)
  for (evaluate in list(loglik(model), score(model), hess_loglik(model))) {
    for (time in list(-1, NA_real_, NaN, Inf, -Inf, "2", factor("2"))) {
      expect_error(evaluate(data.frame(t = time, delta = 1), par = 1),
                   "observation times")
    }
    expect_error(evaluate(data.frame(other = 1, delta = 1), par = 1),
                 "observation times")
    expect_error(evaluate(list(t = 1, delta = 1), par = 1),
                 "observations must be a data frame")
  }
})

test_that("validation retains custom names, missing status, and zero times", {
  model <- dfr_exponential(2)
  model$ob_col <- "elapsed"
  model$delta_col <- "completed"
  df <- data.frame(elapsed = c(0, 1, 2), completed = c(TRUE, FALSE, TRUE))
  expect_equal(loglik(model)(df), 2 * log(2) - 6)
  expect_equal(score(model)(df), -2)
  expect_equal(hess_loglik(model)(df), matrix(-0.5))
  df$completed <- NULL
  expect_equal(loglik(model)(df), 3 * log(2) - 6)
})

test_that("fit rejects invalid records instead of optimizing partial data", {
  model <- dfr_exponential(1)
  expect_error(fit(model)(data.frame(t = c(1, 2), delta = c(1, NA))),
               "event indicators")
})
