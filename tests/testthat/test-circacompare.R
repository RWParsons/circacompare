test_that("circacompare() fits a good model to generated data", {
  df <- make_data(phi1=12)
  out <- circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure")

  both_groups_rhythmic <- as.logical(out[[2]][1,2])
  phase_shift_estimated_within_2hours <- abs(abs(out[[2]][14,2]) - 12) < 2

  expect_true(both_groups_rhythmic)
  expect_true(phase_shift_estimated_within_2hours)
})
