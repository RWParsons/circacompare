test_that("circacompare_mixed() works", {
  phi1_in <- 3.15
  mixed_data <- function(n) {
    counter <- 1
    for (i in 1:n) {
      x <- make_data(k1 = 0, alpha1 = 0, phi1 = rnorm(1, phi1_in, 0.5), hours = 72, noise_sd = 1)
      x$id <- counter
      counter <- counter + 1
      if (i == 1) {
        res <- x
      } else {
        res <- rbind(res, x)
      }
    }
    return(res)
  }

  withr::with_seed(99, {
    df <- mixed_data(20)

    out <- circacompare_mixed(
      x = df,
      col_time = "time",
      col_group = "group",
      col_outcome = "measure",
      col_id = "id",
      control = list(grouped_params = c("phi"), random_params = c("phi1"))
    )
  })

  phi1_fit <- extract_model_coefs(out$fit)["phi1", ]
  phi1_est <- phi1_fit["estimate"]
  phi1_ll <- phi1_est - 1.96 * phi1_fit["std_error"]
  phi1_ul <- phi1_est + 1.96 * phi1_fit["std_error"]
  expect_true(phi1_in < phi1_ul & phi1_in > phi1_ll)
})

### make test that weights are used correctly and malformatted weights are detected
test_that("weights work", {
  phi1_in <- 3.15
  mixed_data <- function(n) {
    counter <- 1
    for (i in 1:n) {
      x <- make_data(k1 = 0, alpha1 = 0, phi1 = rnorm(1, phi1_in, 0.5), hours = 72, noise_sd = 1)
      x$id <- counter
      counter <- counter + 1
      if (i == 1) {
        res <- x
      } else {
        res <- rbind(res, x)
      }
    }
    return(res)
  }

  withr::with_seed(99, {
    df <- mixed_data(20)

    # no weights used (= all weights are 1), hence fit$apVar should not be populated
    out <- circacompare_mixed(
      x = df,
      col_time = "time",
      col_group = "group",
      col_outcome = "measure",
      col_id = "id",
      control = list(grouped_params = c("phi"), random_params = c("phi1"))
    )

    # when weights are not all 1 then fit$apVar should be a matrix
    sw <- runif(n = nrow(df))
    out2 <- circacompare_mixed(
      x = df,
      col_time = "time",
      col_group = "group",
      col_outcome = "measure",
      col_id = "id",
      control = list(grouped_params = c("phi"), random_params = c("phi1")),
      weights = sw
    )

    sw2 <- c(sw, 1)
    expect_error(
      circacompare_mixed(
        x = df,
        col_time = "time",
        col_group = "group",
        col_outcome = "measure",
        col_id = "id",
        control = list(grouped_params = c("phi"), random_params = c("phi1")),
        weights = sw2
      )
    )

    # weights must not contain NA
    sw3 <- sw
    sw3[1] <- NA
    expect_error(
      circacompare_mixed(
        x = df,
        col_time = "time",
        col_group = "group",
        col_outcome = "measure",
        col_id = "id",
        control = list(grouped_params = c("phi"), random_params = c("phi1")),
        weights = sw3
      )
    )

    # weights must not be negative
    sw4 <- sw
    sw4[1] <- -1
    expect_error(
      circacompare_mixed(
        x = df,
        col_time = "time",
        col_group = "group",
        col_outcome = "measure",
        col_id = "id",
        control = list(grouped_params = c("phi"), random_params = c("phi1")),
        weights = sw4
      )
    )
  })

  expect_true(is(out$fit$apVar, "character"))
  expect_true(is(out2$fit$apVar, "matrix"))
})
