test_that("circa_single_mixed() works", {
  set.seed(42)
  mixed_data <- function(n){
    counter <- 1
    for(i in 1:n){
      x <- make_data(k1=rnorm(1, 10, 2), alpha1=0, phi1=0)
      x$id <- counter
      counter <- counter + 1
      if(i==1){res <- x}else{res <- rbind(res, x)}
    }
    return(res)
  }
  df <- mixed_data(n=50)

  out <- circa_single_mixed(x = df, col_time = "time", col_outcome = "measure",
                            col_id = "id", randomeffects = c("k"))
  expect_true(exists("out"))
})
