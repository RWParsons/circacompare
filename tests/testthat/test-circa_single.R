test_that("circa_single works", {
  set.seed(42)
  tau <- 15
  data_rhythmic <- make_data(phi1=0)
  out_rhythmic <- circa_single(x = data_rhythmic, col_time = "time", col_outcome = "measure")

  data_rhythmic$time <- data_rhythmic$time/24*tau
  out_rhythmic_free_tau <-
    circa_single(x = data_rhythmic, col_time = "time", col_outcome = "measure",
                 period=NA,
                 control=list(
                   main_params=c("k", "alpha", "phi", "tau"),
                   period_min=tau-5,
                   period_max=tau+5
                 ))
  out_rhythmic_free_tau$plot + ggplot2::geom_vline(xintercept=out_rhythmic_free_tau$summary$peak_time_hours)
  out_rhythmic_free_tau$model
  data_arrhythmic <- make_data(phi1=12)
  out_arrhythmic <- circa_single(x = data_arrhythmic, col_time = "time", col_outcome = "measure")

  expect_true(class(out_rhythmic) == "list") # no errors when running circa_single()
  expect_true(out_rhythmic$summary[1,2] < 0.01)   # amplitude_p for rhythmic data is small

  expect_true(out_arrhythmic$summary[1,2] > 0.05) # amplitude_p for arrhythmic data is large.


  fit_tau <- summary(out_rhythmic_free_tau$fit)$coef['tau', ]
  tau_est <- fit_tau['Estimate']
  tau_ll <- tau_est - 1.96*fit_tau['Std. Error']
  tau_ul <- tau_est + 1.96*fit_tau['Std. Error']
  expect_true(tau_est < tau_ul & tau_est > tau_ll) # period estimate is approx well estimated to be close to tau (ln 5)

})

