test_that("circacompare() fits a good model to generated data", {
  set.seed(42)
  tau_in <- 15
  phi1_in <- 12
  df <- make_data(phi1=phi1_in)
  out <- circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure")

  df$time <- df$time/24*tau_in
  out_tau_adjusted <- circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure",
                                   period=NA,
                                   control = list(main_params=c("k", "alpha", "tau", "phi"),
                                                  grouped_params=c("k",  "alpha", "phi")))
  out_tau_adjusted$plot

  both_groups_rhythmic <- as.logical(out$table[1, 'value'])
  phase_shift_estimated_within_2hours <- abs(abs(out$table[14,'value']) - phi1_in) < 2

  expect_true(both_groups_rhythmic)
  expect_true(phase_shift_estimated_within_2hours)


  fit_tau <- summary(out_tau_adjusted$fit)$coef['tau', ]
  tau_est <- fit_tau['Estimate']
  tau_ll <- tau_est - 1.96*fit_tau['Std. Error']
  tau_ul <- tau_est + 1.96*fit_tau['Std. Error']
  expect_true(tau_in < tau_ul & tau_in > tau_ll) # period estimate is approx well estimated to be close to tau (ln 5)
})
