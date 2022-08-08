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
