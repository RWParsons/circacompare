test_that("circa_single works", {
  set.seed(42)
  tau <- 15
  data_rhythmic <- make_data(k1 = 0, alpha1 = 0, phi = pi, phi1 = 0, noise_sd = 1)
  out_rhythmic <- circa_single(x = data_rhythmic, col_time = "time", col_outcome = "measure")

  data_rhythmic$time <- data_rhythmic$time / 24 * tau
  out_rhythmic_free_tau <-
    circa_single(
      x = data_rhythmic, col_time = "time", col_outcome = "measure",
      period = NA,
      control = list(
        main_params = c("k", "alpha", "phi", "tau"),
        period_min = tau - 5,
        period_max = tau + 5
      )
    )


  # out_rhythmic_free_tau$plot + ggplot2::geom_vline(xintercept=out_rhythmic_free_tau$summary[5, 'value'])

  data_arrhythmic <- make_data(alpha = 0)
  data_arrhythmic <- data_arrhythmic[data_arrhythmic$group == "g1", ]
  out_arrhythmic <- circa_single(x = data_arrhythmic, col_time = "time", col_outcome = "measure")

  expect_true(class(out_rhythmic) == "list") # no errors when running circa_single()
  expect_true(out_rhythmic$summary[1, 2] < 0.01) # amplitude_p for rhythmic data is small
  expect_true(out_arrhythmic$summary[1, 2] > 0.05) # amplitude_p for arrhythmic data is large.

  fit_tau <- extract_model_coefs(out_rhythmic_free_tau$fit)["tau", ]
  tau_est <- fit_tau["estimate"]
  tau_ll <- tau_est - 1.96 * fit_tau["std_error"]
  tau_ul <- tau_est + 1.96 * fit_tau["std_error"]
  expect_true(tau < tau_ul & tau > tau_ll) # period estimate is approx well estimated to be close to tau (ln 5)


  # assess whether decay on amplitude per-hour is modelled well when period is parameterized
  alpha_decay_in <- 0.01
  tau_in <- 16
  df <- make_data(k = 5, k1 = 0, alpha = 20, alpha1 = 0, phi = 2, phi1 = 0, hours = 96, noise_sd = 1)
  df$time <- df$time / 24 * tau_in
  df$measure <- df$measure * exp(-alpha_decay_in * (df$time))
  out_alpha_decay <- circa_single(
    x = df,
    col_time = "time",
    col_outcome = "measure",
    period = NA,
    control = list(
      main_params = c("k", "alpha", "phi", "tau"),
      decay_params = c("alpha"),
      period_min = 12,
      period_max = 20
    )
  )
  out_alpha_decay

  fit_alpha_decay <- extract_model_coefs(out_alpha_decay$fit)["alpha_decay", ]
  alpha_decay_est <- fit_alpha_decay["estimate"]
  alpha_decay_ll <- alpha_decay_est - 1.96 * fit_alpha_decay["std_error"]
  alpha_decay_ul <- alpha_decay_est + 1.96 * fit_alpha_decay["std_error"]
  expect_true(alpha_decay_in < alpha_decay_ul & alpha_decay_in > alpha_decay_ll)
})


### make test to capture output and test that running with/without suppress_all works to suppress messages to console
test_that("suppress_all works", {
  y <- structure(list(
    time = c(
      1L, 1L, 1L, 1L, 5L, 5L, 9L, 9L, 13L,
      13L, 17L, 17L, 17L, 21L, 21L, 21L
    ),
    value = c(
      6.46491702175632,
      6.37210528510888, 6.75505623236344, 6.4457897862926, 6.63766950190431,
      6.48725138475295, 6.40819847507183, 6.42253808100338, 6.37486222182972,
      6.51868394085349, 6.41506838906571, 6.40449437273951, 6.47273627195726,
      6.76905314588271, 6.59233676207294, 6.44481187866212
    )
  ),
  class = "data.frame", row.names = c(NA, -16L)
  )

  set.seed(1)
  output <- capture.output(
    {
      circa_single(x = y, col_time = "time", col_outcome = "value", return_figure = FALSE)
    },
    type = "message"
  )
  expect_true(length(output) > 1)

  output <- capture.output(
    {
      circa_single(x = y, col_time = "time", col_outcome = "value", return_figure = FALSE, suppress_all = TRUE)
    },
    type = "message"
  )
  expect_true(length(output) == 0)
})
