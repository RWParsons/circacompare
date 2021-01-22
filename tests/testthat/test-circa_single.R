test_that("circa_single works", {
  data_rhythmic <- make_data(phi1=0)
  data_arrhythmic <- make_data(phi1=12)

  out_rhythmic <- circa_single(x = data_rhythmic, col_time = "time", col_outcome = "measure")
  out_arrhythmic <- circa_single(x = data_arrhythmic, col_time = "time", col_outcome = "measure")

  expect_true(class(out_rhythmic) == "list")
  expect_true(out_rhythmic[[2]][3] < 0.1)
  expect_true(out_arrhythmic[[2]][3] > 0.05)
})
