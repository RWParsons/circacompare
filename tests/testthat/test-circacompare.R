test_that("circacompare() fits a good model to generated data", {
  tau_in <- 15
  phi1_in <- 12

  withr::with_seed(42, {
    df <- make_data(phi1 = (phi1_in / 24) * (2 * pi), noise_sd = 2)
    out <- circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure")

    df$time <- df$time / 24 * tau_in
    out_tau_adjusted <- circacompare(
      x = df, col_time = "time", col_group = "group", col_outcome = "measure",
      period = NA,
      control = list(
        main_params = c("k", "alpha", "tau", "phi"),
        grouped_params = c("k", "alpha", "phi"),
        period_min = tau_in - 4, period_max = tau_in + 4
      )
    )
  })

  both_groups_rhythmic <- as.logical(out$summary[1, "value"] < 0.05 & out$summary[2, "value"] < 0.05)
  phase_shift_estimated_within_2hours <- abs(abs(out$summary[13, "value"]) - phi1_in) < 2

  expect_true(both_groups_rhythmic)
  expect_true(phase_shift_estimated_within_2hours)

  fit_tau <- extract_model_coefs(out_tau_adjusted$fit)["tau", ]
  tau_est <- fit_tau["estimate"]
  tau_ll <- tau_est - 1.96 * fit_tau["std_error"]
  tau_ul <- tau_est + 1.96 * fit_tau["std_error"]
  expect_true(tau_in < tau_ul & tau_in > tau_ll) # period estimate is approx well estimated to be close to tau (ln 5)


  # create some longer time period data and keep all parameters the same except amplitude
  # create some decay in one group for amplitude and test whether it's well estimated by the model.

  withr::with_seed(42, {
    df <- make_data(k1 = 0, alpha1 = 10, phi1 = 0, seed = 42, hours = 96, noise_sd = 2)
    df$time <- df$time / 24 * tau_in
    alpha_decay1_in <- 0.05
    # note that decay is on a scale of time in radians, not time in hours.
    df$measure[df$group == "g2"] <- df$measure[df$group == "g2"] * exp(-alpha_decay1_in * df$time[df$group == "g2"])

    out_alpha_decay <-
      circacompare(
        x = df, "time", "group", "measure", period = NA,
        control = list(
          main_params = c("k", "alpha", "phi", "tau"),
          decay_params = c("alpha"),
          grouped_params = c("alpha", "alpha_decay"),
          period_min = tau_in - 4, period_max = tau_in + 4
        )
      )
  })

  fit_alpha_decay1 <- extract_model_coefs(out_alpha_decay$fit)["alpha_decay1", ]
  alpha_decay1_est <- fit_alpha_decay1["estimate"]
  alpha_decay1_ll <- alpha_decay1_est - 1.96 * fit_alpha_decay1["std_error"]
  alpha_decay1_ul <- alpha_decay1_est + 1.96 * fit_alpha_decay1["std_error"]
  expect_true(alpha_decay1_in < alpha_decay1_ul & alpha_decay1_in > alpha_decay1_ll)
})

### make test that weights are used correctly and malformatted weights are detected
test_that("weights work", {
  # all weights should be 1
  df <- make_data(phi1 = 6)
  out <- circacompare(
    x = df, col_time = "time", col_outcome = "measure", col_group = "group"
  )
  expect_true(all(out$fit$weights == 1))

  # all weights should not be 1
  sw <- runif(n = nrow(df))
  out2 <- circacompare(
    x = df, col_time = "time", col_outcome = "measure", col_group = "group",
    weights = sw
  )
  expect_false(all(out2$fit$weights == 1))

  # weights must be same length as nrow(x)
  sw2 <- c(sw, 1)
  expect_error(
    circacompare(
      x = df, col_time = "time", col_outcome = "measure", col_group = "group",
      weights = sw2
    )
  )

  sw3 <- sw
  sw3[1] <- NA
  expect_error(
    circacompare(
      x = df, col_time = "time", col_outcome = "measure", col_group = "group",
      weights = sw3
    )
  )
})
