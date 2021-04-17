test_that("circa_single works", {
  set.seed(42)
  tau <- 15
  data_rhythmic <- make_data(k1=0, alpha1=0, phi=pi, phi1=0, noise_sd=1)
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


  out_rhythmic_free_tau$plot + ggplot2::geom_vline(xintercept=out_rhythmic_free_tau$summary[5, 'value'])
  out_rhythmic_free_tau$model
  data_arrhythmic <- make_data(alpha=0)
  data_arrhythmic <- data_arrhythmic[data_arrhythmic$group == "g1",]
  out_arrhythmic <- circa_single(x = data_arrhythmic, col_time = "time", col_outcome = "measure")

  expect_true(class(out_rhythmic) == "list") # no errors when running circa_single()
  expect_true(out_rhythmic$summary[1,2] < 0.01)   # amplitude_p for rhythmic data is small
  expect_true(out_arrhythmic$summary[1,2] > 0.05) # amplitude_p for arrhythmic data is large.

  fit_tau <- extract_model_coefs(out_rhythmic_free_tau$fit)['tau',]
  tau_est <- fit_tau['estimate']
  tau_ll <- tau_est - 1.96*fit_tau['std_error']
  tau_ul <- tau_est + 1.96*fit_tau['std_error']
  expect_true(tau < tau_ul & tau > tau_ll) # period estimate is approx well estimated to be close to tau (ln 5)

})

