test_that("circa_single_mixed() works", {
  set.seed(42)
  tau_in <- 14
  mixed_data <- function(n) {
    counter <- 1
    for (i in 1:n) {
      x <- make_data(k1 = rnorm(1, 10, 1), alpha1 = 0, phi1 = 0)
      x$id <- counter
      counter <- counter + 1
      if (i == 1) {
        res <- x
      } else {
        res <- rbind(res, x)
      }
    }
    return(res)
  }
  df <- mixed_data(n = 12)

  out <- circa_single_mixed(
    x = df, col_time = "time", col_outcome = "measure",
    col_id = "id", randomeffects = c("k"),
    control = list(decay_params = c("k"))
  )

  expect_true(class(out) == "list")
  expect_true(round(extract_model_coefs(out$fit)["alpha", "estimate"]) == 10)

  df$time <- df$time / 24 * tau_in
  out_tau_adjusted <-
    circa_single_mixed(
      x = df, col_time = "time", col_outcome = "measure",
      col_id = "id", randomeffects = c("k"), period = NA,
      control = list(period_param = T, period_min = tau_in - 4, period_max = tau_in + 4)
    )

  fit_tau <- extract_model_coefs(out_tau_adjusted$fit)["tau", ]
  tau_est <- fit_tau["estimate"]
  tau_ll <- tau_est - 1.96 * fit_tau["std_error"]
  tau_ul <- tau_est + 1.96 * fit_tau["std_error"]
  expect_true(tau_in < tau_ul & tau_in > tau_ll)
})

### make test that weights are used correctly and malformatted weights are detected
test_that("weights work", {
  set.seed(42)
  mixed_data <- function(n) {
    counter <- 1
    for (i in 1:n) {
      x <- make_data(k1 = rnorm(1, 10, 2), alpha1 = 0, phi1 = 0)
      x$id <- counter
      counter <- counter + 1
      if (i == 1) {
        res <- x
      } else {
        res <- rbind(res, x)
      }
    }
    return(res)
  }

  df <- mixed_data(n = 50)

  # no weights used (= all weights are 1), hence fit$apVar should not be populated
  out <- circa_single_mixed(
    x = df, col_time = "time", col_outcome = "measure",
    col_id = "id", randomeffects = c("k")
  )
  expect_true(is(out$fit$apVar, "character"))

  # when weights are not all 1 then fit$apVar should be a matrix
  sw <- runif(n = nrow(df))
  out2 <- circa_single_mixed(
    x = df, col_time = "time", col_outcome = "measure",
    col_id = "id", randomeffects = c("k"), weights = sw
  )
  expect_true(is(out2$fit$apVar, "matrix"))

  # weights must be same length as nrow(x)
  sw2 <- c(sw, 1)
  expect_error(
    circa_single_mixed(
      x = df, col_time = "time", col_outcome = "measure",
      col_id = "id", randomeffects = c("k"), weights = sw2
    )
  )

  # weights must not contain NA
  sw3 <- sw
  sw3[1] <- NA
  expect_error(
    circa_single_mixed(
      x = df, col_time = "time", col_outcome = "measure",
      col_id = "id", randomeffects = c("k"), weights = sw3
    )
  )

  # weights must not be negative
  sw4 <- sw
  sw4[1] <- -1
  expect_error(
    circa_single_mixed(
      x = df, col_time = "time", col_outcome = "measure",
      col_id = "id", randomeffects = c("k"), weights = sw4
    )
  )
})
