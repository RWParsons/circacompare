test_that("circacompare_mixed() works", {
  set.seed(99)
  tau_in <- 8
  phi1_in <- 6
  mixed_data <- function(n){
    counter <- 1
    for(i in 1:n){
      x <- make_data(k1=0, alpha1=0, phi1=rnorm(1, phi1_in, 1))
      x$id <- counter
      counter <- counter + 1
      if(i==1){res <- x}else{res <- rbind(res, x)}
    }
    return(res)
  }

  x1 <- mixed_data(n=10)
  x2 <- mixed_data(n=10)
  x2$time <- x2$time + 24
  df <- rbind(x1, x2)

  out <- circacompare_mixed(
    x = df,
    col_time = "time",
    col_group = "group",
    col_outcome = "measure",
    col_id = "id",
    randomeffects = c("phi1"), verbose=F
  )


  phi1_est <- extract_model_coefs(out$fit)['phi1', 'estimate']
  phi1_est_close_to_generator <- (abs(phi1_est) - (phi1_in/(2*pi))) < 1
  expect_true(phi1_est_close_to_generator)

  out$plot + ggplot2::geom_smooth(ggplot2::aes(group=interaction(as.character(id), group)), span=0.2, se=F)
})
