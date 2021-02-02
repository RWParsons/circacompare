test_that("circacompare_mixed() works", {
  set.seed(99)
  mixed_data <- function(n){
    counter <- 1
    for(i in 1:n){
      x <- make_data(k1=0, alpha1=0, phi1=rnorm(1, 6, 1))
      x$id <- counter
      counter <- counter + 1
      if(i==1){res <- x}else{res <- rbind(res, x)}
    }
    return(res)
  }
  df <- mixed_data(n=50)

  out <- circacompare_mixed(x = df, col_time = "time", col_group = "group",
                            col_outcome = "measure", col_id = "id",
                            randomeffects = c("phi", "phi1"))

  phi1_est <- summary(out[[3]])$tTable[6,1]
  phi1_est_close_to_generator <- (abs(phi1_est) - 0.5*pi) < 1
  expect_true(phi1_est_close_to_generator)

})
