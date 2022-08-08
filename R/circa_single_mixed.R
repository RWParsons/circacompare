#' @title circa_single_mixed
#' @name circa_single_mixed
#'
#' @description \code{circa_single_mixed} is similar to \code{circa_single} but allows for some simple, user-specified random-effects on the rhythmic parameters of choice.
#'
#' @param x data.frame.  This is the data.frame which contains the rhythmic data in a tidy format.
#' @param col_time The name of the column within the data.frame, x, which contains time in hours at which the data were collected.
#' @param col_outcome The name of the column within the data.frame, x, which contains outcome measure of interest.
#' @param col_id The name of the column within the data.frame, \code{x}, which contains the identifying values for the random effect, such as \code{subject_id}.
#' @param randomeffects which rhythmic parameters to allow random effects. The default is to include all rhythmic parameters.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, \code{24}.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is to \code{0.05}.
#' @param nlme_control A list of control values for the estimation algorithm to replace the default values returned by the function nlme::nlmeControl. Defaults to an empty list.
#' @param nlme_method A character string. If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "ML".
#' @param verbose An optional logical value. If \code{TRUE} information on the evolution of the iterative algorithm is printed. Default is \code{FALSE}.
#' @param timeout_n The upper limit for the model fitting attempts. Default is \code{10000}.
#' @param return_figure Whether or not to return a ggplot graph of the rhythm and cosine model.
#' @param control \code{list}. Used to control the parameterization of the model.
#'
#' @return list
#' @export
#'
#' @examples
#' set.seed(42)
#' mixed_data <- function(n){
#'   counter <- 1
#'   for(i in 1:n){
#'       x <- make_data(k1=rnorm(1, 10, 2), alpha1=0, phi1=0)
#'       x$id <- counter
#'       counter <- counter + 1
#'       if(i==1){res <- x}else{res <- rbind(res, x)}
#'   }
#'   return(res)
#' }
#' df <- mixed_data(n=50)
#' out <- circa_single_mixed(x = df, col_time = "time", col_outcome = "measure",
#'                          col_id = "id", randomeffects = c("k"))
circa_single_mixed <- function (x,
                                col_time,
                                col_outcome,
                                col_id,
                                randomeffects = c("k", "alpha", "phi"),
                                period = 24,
                                alpha_threshold = 0.05,
                                nlme_control = list(),
                                nlme_method = "ML",
                                verbose = FALSE,
                                timeout_n = 10000,
                                return_figure = TRUE,
                                control = list()){
  if(!requireNamespace("ggplot2", quietly = TRUE) & return_figure){
    return(message("Please install 'ggplot2'"))
  }

  controlVals <- circa_single_mixed_control()
  controlVals[names(control)] <- control

  if(controlVals$period_param){
    controlVals$main_params <- c(controlVals$main_params, "tau")
  }
  if("tau" %in% controlVals$main_params){
    controlVals$period_param=TRUE
  }

  randomeffects <- unique(c(randomeffects, controlVals$random_params))

  x <- x[c(col_time, col_outcome, col_id)]
  colnames(x) <- c("time", "measure", "id")

  if(length(controlVals$decay_params)>0){
    p <- append(controlVals$main_params, paste0(controlVals$decay_params, "_decay"))
  }else{
    p <- controlVals$main_params
  }

  if(length(setdiff(randomeffects, p)) != 0){
    return(message('"randomeffects" should only include the names of parameters\nthat represent rhythmic characteristics in the model.\nThey should be a subset of, or equal to c("k", "alpha", "phi")'))
  }

  if(length(randomeffects) == 0){
    return(message("If you do not want to include any random effects, than you ought to use 'circa_single' rather than 'circa_single_mixed'"))
  }

  if (!class(x$time) %in% c("numeric", "integer")) {
    return(message(paste("The time variable which you gave was a '",
                         class(x$time), "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the time variable in your dataframe to be of one of these classes",
                         sep = "")))
  }

  if (!class(x$measure) %in% c("numeric", "integer")) {
    return(message(paste("The measure variable which you gave was a '",
                         class(x$measure), "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the measure variable in your dataframe to be of one of these classes",
                         sep = "")))
  }

  if(!class(period) %in% c("numeric", "integer") & !controlVals$period_param){
    return(message(paste0("The period argument must be a number representing the period of the rhythm in hours\n",
                          "If you would like the period to be estimated as part of the model, use:\ncontrol=list(period_param=TRUE)")))
  }

  if(controlVals$period_param & !is.na(period)){
    message(paste0("control$period_param is TRUE\n'period=", period, "' is being ignored.\nSet 'period=NA' to avoid this message"))
  }

  x$time_r <- x$time*2*pi
  if(!controlVals$period_param){
    x$period <- period
  }else{
    if(is.null(controlVals$period_min) | is.null(controlVals$period_min)){
      message(paste0("If you want the model to estimate the period using a parameter,",
                     "you may get faster convergence if you provide an approximate range using 'period_min' and 'period_max' in control()",
                     "\nCurrently assuming period is between: period_min=", controlVals$period_min,
                     "and period_max=", controlVals$period_max))
    }
  }

  success <- FALSE
  n <- 0
  form <- create_formula(main_params = controlVals$main_params, decay_params=controlVals$decay_params)$formula
  fixedeffects_formula <- stats::formula(paste(paste0(p, collapse="+"), "~ 1"))
  randomeffects_formula <- stats::formula(paste(paste0(randomeffects, collapse="+"), "~ 1 | id"))

  while(!success){
    fit.nlme <- try({nlme::nlme(model = form,
                                random = randomeffects_formula,
                                fixed = fixedeffects_formula,
                                data = x,
                                start = unlist(start_list(outcome=x$measure, controlVals=controlVals)),
                                control = nlme_control,
                                method = nlme_method,
                                verbose = verbose)},
                    silent = ifelse(verbose, FALSE, TRUE)
    )

    if ("try-error" %in% class(fit.nlme)) {
      n <- n + 1
    }
    else {
      nlme_coefs <- extract_model_coefs(fit.nlme)
      V <- nlme_coefs[, 'estimate']
      success <- assess_model_estimates(param_estimates=V, controlVals=controlVals)
      n <- n + 1
    }
    if (n >= timeout_n) {
      return(message("Failed to converge data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  data_rhythmic <- nlme_coefs['alpha', 'p_value'] < alpha_threshold

  if(!controlVals$period_param){V['tau'] <- period}

  eq_expression <- create_formula(main_params = controlVals$main_params, decay_params=controlVals$decay_params)$f_equation
  eval(parse(text=eq_expression))

  if(return_figure){
    if(data_rhythmic){
      fig_out <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
        ggplot2::stat_function(fun = eq, size = 1) +
        ggplot2::geom_point() +
        ggplot2::xlim(min(floor(x$time/period) * period),
                      max(ceiling(x$time/period) * period)) +
        ggplot2::labs(subtitle = "Data is rhythmic", x = "time (hours)")
    }else{
      fig_out <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
        ggplot2::geom_point() +
        ggplot2::xlim(min(floor(x$time/period) * period),
                      max(ceiling(x$time/period) * period)) +
        ggplot2::labs(subtitle = "Data is arrhythmic", x = "time (hours)")
    }
  }

  output_parms <- circa_summary(model=fit.nlme, period=period, control=controlVals)

  if(return_figure){
    return(list(fit=fit.nlme, summary=output_parms, plot=fig_out))
  }else{
    return(list(fit=fit.nlme, summary=output_parms))
  }
}


circa_single_mixed_control <- function(period_param=F, period_min=20, period_max=28,
                                       main_params=c("k", "alpha", "phi"), decay_params=c(),
                                       random_params=c("k")){
  list(period_param=period_param, period_min=period_min, period_max=period_max,
       main_params=main_params, decay_params=decay_params, random_params=random_params)
}
