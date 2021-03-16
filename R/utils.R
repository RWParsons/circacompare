extract_model_coefs <- function(model){
  if("nls" %in% class(model)){
    results <- summary(model)$coef[, c(1, 4)]
    dimnames(results)[[2]] <- c("estimate", "p_value")
  }
  if("nlme" %in% class(model)){
    results <- summary(model)$tTable[, c(1, 5)]
    dimnames(results)[[2]] <- c("estimate", "p_value")
  }
  return(results)
}


model_each_group <- function(data, type, args=list(timeout_n=10000, alpha_threshold=0.05)){
  success <- FALSE
  n <- 0
  while(!success){
    starting_params <- c(
      k=mean(data$measure, na.rm = TRUE)*2*stats::runif(1),
      alpha=(max(data$measure, na.rm = TRUE) - min(data$measure, na.rm = TRUE)) * stats::runif(1),
      phi=stats::runif(1)*6.15 - 3.15
    )

    if(type=="nlme"){
      ranefs <- intersect(c("k", "alpha", "phi"), args$randomeffects)
      if(length(ranefs)!=0){
        ranefs_formula <- stats::formula(paste(paste0(ranefs, collapse="+"), "~ 1 | id"))
        }else{
          ranefs_formula <- NULL
      }
      fit <- try({nlme::nlme(model = measure~k+alpha*cos(time_r-phi),
                             random = ranefs_formula,
                             fixed = k+alpha+phi~1,
                             data = data,
                             start = starting_params,
                             method = args$nlme_method,
                             control = args$nlme_control,
                             verbose = args$verbose)},
                 silent = ifelse(args$verbose, FALSE, TRUE)
      )
    }else if(type=="nls"){
      fit <- try({stats::nls(measure~k + alpha*cos(time_r-phi),
                             data = data,
                             start = lapply(split(starting_params, names(starting_params)), unname))},
                 silent = TRUE)
    }

    if("try-error" %in% class(fit)){
      n <- n + 1
    }else{
      coefs <- extract_model_coefs(fit)
      success <- ifelse(coefs['alpha', 'estimate'] > 0 & coefs['phi', 'estimate'] > 0 & coefs['phi', 'estimate'] < 2*pi ,TRUE,FALSE)
      n <- n + 1

    }
    if(n >= args$timeout_n){
      return(list(
        model=NA,
        timeout=TRUE,
        rhythmic=NA,
        k_estimate=NA,
        alpha_estimate=NA,
        phi_estimate=NA
      ))
    }
  }

  return(list(
    model=fit,
    timeout=FALSE,
    rhythmic=coefs['alpha', 'p_value'] < args$alpha_threshold,
    k_estimate=coefs['k', 'estimate'],
    alpha_estimate=coefs['alpha', 'estimate'],
    alpha_p=coefs['alpha', 'p_value'],
    phi_estimate=coefs['phi', 'estimate']
  ))
}


random_start_phi1 <- function(p){
  p <- stats::rnorm(1, p, 2)
  if(abs(p) > pi){
    if(p < 0){
      p <- p + 2*pi
    }else{
      p <- p - 2*pi
    }
  }
  return(p)
}
