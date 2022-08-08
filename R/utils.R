utils::globalVariables(c('time', 'measure', 'group', 'eq', 'eq_1', 'eq_2'))


extract_model_coefs <- function(model){
  if("nls" %in% class(model)){
    results <- summary(model)$coef[, c(1, 2, 4)]
  }
  if("nlme" %in% class(model)){
    results <- summary(model)$tTable[, c(1, 2, 5)]
  }
  dimnames(results)[[2]] <- c("estimate", "std_error", "p_value")
  return(results)
}


model_each_group <- function(data, type, form=stats::as.formula("measure~k+alpha*cos(time_r-phi)"),
                             controlVals=list(),
                             args=list(timeout_n=10000, alpha_threshold=0.05)){
  success <- FALSE
  n <- 0
  while(!success){
    starting_params <- start_list(outcome=data$measure, controlVals=controlVals)
    # use the nls function below if the only random effects are on group parameters
    if(type=="nlme" & length(args$randomeffects)>0){
      ranefs <- intersect(controlVals$non_grouped_params, args$randomeffects)
      ranefs_formula <- stats::formula(paste(paste0(ranefs, collapse="+"), "~ 1 | id"))
      fixefs_formula <- stats::formula(paste(paste0(controlVals$non_grouped_params, collapse="+"), "~ 1"))

      fit <- try({nlme::nlme(model = form,
                             random = ranefs_formula,
                             fixed = fixefs_formula,
                             data = data,
                             start = unlist(starting_params),
                             method = args$nlme_method,
                             control = args$nlme_control,
                             verbose = args$verbose)},
                 silent = ifelse(args$verbose, FALSE, TRUE)
      )
    }else{
      fit <- try({stats::nls(formula=form,
                             data = data,
                             start = starting_params)
        },
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


create_formula <- function(main_params=c("k", "alpha", "phi"), decay_params=c(), grouped_params=c()){
  if(length(grouped_params[!grouped_params %in% c(main_params, paste0(decay_params, "_decay"))])>0){
    stop("All grouped parameters must be within the main or decay parameters.")
  }
  if(all(c("phi", "tau") %in% grouped_params)){
    stop("Either phase or period can have grouped effects, but not both.")
  }

  build_component <- function(pattern, eq=F){
    t <- "time_r"

    main <- main_params[grep(pattern, main_params)]
    decay <- decay_params[grep(pattern, decay_params)]
    grouped <- grouped_params[grep(pattern, grouped_params)]
    if(length(main)<length(decay)|(length(main)+length(decay))<length(grouped)){
      stop(paste0("You can't model the decay of a parameter which doesn't exist in the model!\nAdd `",
                  gsub("^.", "", pattern), "` to the `main_params`"))
    }
    if(eq & length(main)!=0){
      main <- paste0("V['", main, "']")
    }

    group_match <- grouped
    group_main <- group_match[!grepl("decay", group_match)]
    if(length(group_main)!=0){
      if(eq){
        group_main <- paste0("V['", group_main, "1']")
      }else{
        group_main <- paste0(group_main, "1")
      }
      main <- paste0("(",main, "+", group_main, ")")
    }

    if(length(decay)!=0){
      decay <- paste0(decay, "_decay")
      group_decay <- group_match[grepl("decay", group_match)]

      if(length(group_decay)!=0){
        if(eq){
          decay <- gsub(group_decay, paste0("V['",group_decay, "']", '+', "V['", group_decay, "1']"), decay)
        }else{
          decay <- gsub(group_decay, paste0(group_decay, '+', group_decay, "1"), decay)
        }
      }else{
        if(eq){
          decay <- paste0("V['", decay, "']")
        }
      }
      main <- paste0(main, "*exp(-(", decay, ")*time)")
    }
    component <- main

    if(pattern=="^tau"){
      if(length(component)==0){
        component <- ifelse(eq, "(1/V['tau'])*", "(1/period)*")
      }else{
        component <- paste0("(1/(", component, "))*")
      }
    }
    return(component)
  }

  res_formula <- paste0("measure~", build_component("^k"), "+(", build_component("^alpha"),
                        ")*cos(", build_component("^tau"), "time_r-(", build_component("^phi"), "))")
  res_formula <- gsub("(?<=[a-z])1", "1*x_group", res_formula, perl=T)
  res_formula <- stats::as.formula(res_formula)
  res_equation <- paste0(build_component("^k", eq=T), "+(", build_component("^alpha", eq=T),
                         ")*cos(", build_component("^tau", eq=T), "time_r-(", build_component("^phi", eq=T), "))")


  if(length(grouped_params)>0){
    res_equation <- list(
      g1=paste0("eq_1 <- function(time) {time_r<-time*2*pi;return(", gsub("\\+V\\['[a-z]*_?[a-z]*1'\\]", "", res_equation), ")}"),
      g2=paste0("eq_2 <- function(time) {time_r<-time*2*pi;return(", res_equation, ")}")
    )
  }else{
    res_equation <- paste0("eq <- function(time) {time_r<-time*2*pi;return(", res_equation, ")}")
  }
  return(list(formula=res_formula, f_equation=res_equation))
}


start_list <- function(outcome, controlVals){
  res <-
    list(
      k=mean(outcome, na.rm = TRUE) * 2 * stats::runif(1),
      alpha=(max(outcome, na.rm = TRUE) - min(outcome, na.rm = TRUE)) * stats::runif(1),
      phi=stats::runif(1) * 6.15 - 3.15
    )
  if("tau" %in% controlVals$main_params){
    res <-
      append(
        res,
        list(tau=stats::runif(1, min=controlVals$period_min, max=controlVals$period_max))
      )
  }
  if(length(controlVals$decay_params)!=0){
    vec <- stats::rnorm(length(controlVals$decay_params), mean=0, sd=0.2)
    names(vec) <- paste0(controlVals$decay_params, "_decay")
    res <-
      append(
        res,
        as.list(vec)
      )
  }
  return(res)
}


start_list_grouped <- function(g1, g2, grouped_params=c("k", "alpha", "phi")){
  V_g1 <- extract_model_coefs(g1)[, 'estimate']
  se_g1 <- extract_model_coefs(g1)[, 'std_error']
  V_g2 <- extract_model_coefs(g2)[, 'estimate']
  se_g2 <- extract_model_coefs(g2)[, 'std_error']

  params <- names(V_g1)

  vec <- rep(0, length(params))
  names(vec) <- params

  for(gp in grouped_params){
    vec[gp] <- V_g1[gp]*2*stats::runif(1)
    if(gp=="phi"){
      vec[paste0(gp, "1")] <- random_start_phi1(p=(V_g2[gp] - V_g1[gp]))
    }else{
      vec[paste0(gp, "1")] <- (V_g2[gp] - V_g1[gp])*2*stats::runif(1)
    }
  }

  for(shared_param in params[!params %in% grouped_params]){
    vec[shared_param] <- stats::runif(n=1, V_g1[shared_param], V_g1[shared_param])
  }

  order <- c(
    names(vec)[grepl("^k", names(vec))],
    names(vec)[grepl("^alpha", names(vec))],
    names(vec)[grepl("^tau", names(vec))],
    names(vec)[grepl("^phi", names(vec))]
  )
  lst <- as.list(vec)

  lst <- lst[order]
  return(lst)
}


assess_model_estimates <- function(param_estimates, controlVals){
  res <- TRUE

  if("alpha" %in% names(param_estimates)){
    res <- ifelse(param_estimates['alpha']<0, FALSE, res)
  }

  if("alpha1" %in% names(param_estimates)){
    res <- ifelse(param_estimates['alpha'] + param_estimates['alpha1']<0, FALSE, res)
  }

  if("phi" %in% names(param_estimates)){
    if("phi1" %in% names(param_estimates)){
      res <- ifelse(param_estimates['phi1']<(-pi) | param_estimates['phi1']>pi,
                    FALSE, res)
    }else{
      res <- ifelse(param_estimates['phi']<0 | param_estimates['phi']>2*pi,
                    FALSE, res)
    }
  }

  fx_check_period <- function(period) {
    res <- TRUE
    if(!is.null(controlVals$period_min)) {
      res <- ifelse(period < controlVals$period_min, FALSE, res)
    }
    if(!is.null(controlVals$period_max)) {
      res <- ifelse(period > controlVals$period_max, FALSE, res)
    }

    res <- ifelse(period < 0, FALSE, res)
    return(res)
  }

  if("tau" %in% names(param_estimates)){
    res <- ifelse(!fx_check_period(param_estimates['tau']), FALSE, res)
  }

  if("tau1" %in% names(param_estimates)){
    res <- ifelse(!fx_check_period(param_estimates['tau'] + param_estimates['tau1']), FALSE, res)
  }

  return(res)
}


circa_summary <- function(model, period, control,
                          g1=NULL, g2=NULL, g1_text=NULL, g2_text=NULL){
  grouped <- ifelse(length(control$grouped_params)>0, TRUE, FALSE)
  row_adder <- function(data, parameter, value){
    new_row <- data.frame(parameter=parameter, value=value)
    row.names(new_row) <- NULL
    return(rbind(data, new_row))
  }

  coefs <- extract_model_coefs(model)
  V <- coefs[, 'estimate']

  if(!"tau" %in% names(V)){
    V["tau"] <- period
  }

  if("phi" %in% names(V)){
    if(V['phi'] > 2*pi){
      while(V['phi'] > 2*pi){
        V['phi'] <- V['phi'] - 2*pi
      }
    }
    if(V['phi'] < -2*pi){
      while(V['phi'] < -2*pi){
        V['phi'] <- V['phi'] + 2*pi
      }
    }
  }

  peak_time_diff <- V['phi1']*V['tau']/(2*pi)

  res <- data.frame(parameter=c(), value=c())
  if(!grouped){
    res <- row_adder(res, "rhythmic_p", coefs['alpha', 'p_value'])
    res <- row_adder(res, "mesor", V['k'])
    if("k_decay" %in% names(V)){res <- row_adder(res, "mesor_decay", V['k_decay'])}
    res <- row_adder(res, "amplitude", V['alpha'])
    if("alpha_decay" %in% names(V)){res <- row_adder(res, "alpha_decay", V['alpha_decay'])}
    res <- row_adder(res, "phase_radians", V['phi'])
    if("phi_decay" %in% names(V)){res <- row_adder(res, "phase_radians_decay", V['phi_decay'])}
    res <- row_adder(res, "peak_time_hours", (V['phi']/(2*pi)) * V['tau'])
    res <- row_adder(res, "period", V['tau'])
    if("tau_decay" %in% names(V)){res <- row_adder(res, "period_decay", V['tau_decay'])}
    return(res)
  }

  res <- row_adder(res, paste0("Presence of rhythmicity (p-value) for ", g1_text), g1$alpha_p)
  res <- row_adder(res, paste0("Presence of rhythmicity (p-value) for ", g2_text), g2$alpha_p)

  if("k1" %in% names(V)){
    res <- row_adder(res, paste0(g1_text, " mesor estimate"), V['k'])
    res <- row_adder(res, paste0(g2_text, " mesor estimate"), V['k']+V['k1'])
    res <- row_adder(res, "Mesor difference estimate", V['k1'])
    res <- row_adder(res, "P-value for mesor difference", coefs['k1', 'p_value'])
  }else{
    res <- row_adder(res, "Shared mesor estimate", V['k'])
  }
  if("k_decay" %in% names(V)){
    if("k_decay1" %in% names(V)){
      res <- row_adder(res, paste0(g1_text, " mesor decay estimate"), V['k_decay'])
      res <- row_adder(res, paste0(g2_text, " mesor decay estimate"), V['k_decay']+V['k_decay1'])
      res <- row_adder(res, "Mesor decay difference estimate", V['k_decay1'])
      res <- row_adder(res, "P-value for mesor decay difference", coefs['k_decay1', 'p_value'])
    }else{
      res <- row_adder(res, "Shared mesor decay estimate", V['k_decay'])
    }
  }

  if("alpha1" %in% names(V)){
    res <- row_adder(res, paste0(g1_text, " amplitude estimate"), V['alpha'])
    res <- row_adder(res, paste0(g2_text, " amplitude estimate"), V['alpha']+V['alpha1'])
    res <- row_adder(res, "Amplitude difference estimate", V['alpha1'])
    res <- row_adder(res, "P-value for amplitude difference", coefs['alpha1', 'p_value'])
  }else{
    res <- row_adder(res, "Shared amplitude estimate", V['alpha'])
  }
  if("alpha_decay" %in% names(V)){
    if("alpha_decay1" %in% names(V)){
      res <- row_adder(res, paste0(g1_text, " amplitude decay estimate"), V['alpha_decay'])
      res <- row_adder(res, paste0(g2_text, " amplitude decay estimate"), V['alpha_decay']+V['alpha_decay1'])
      res <- row_adder(res, "Amplitude decay difference estimate", V['alpha_decay1'])
      res <- row_adder(res, "P-value for amplitude decay difference", coefs['alpha_decay1', 'p_value'])
    }else{
      res <- row_adder(res, "Shared amplitude decay estimate", V['alpha_decay'])
    }
  }

  peak_time <- V['phi']*V['tau']/(2*pi)
  while(peak_time > V['tau'] | peak_time < 0){
    if(peak_time > V['tau']){
      peak_time <- peak_time - V['tau']
    }
    if(peak_time<0){
      peak_time <- peak_time + V['tau']
    }
  }
  if("phi1" %in% names(V)){
    peak_time_diff <- V['phi1']*V['tau']/(2*pi)
    g2_peak_time <- (V['phi']+V['phi1'])*V['tau']/(2*pi)
    while(g2_peak_time>V['tau']| g2_peak_time < 0){
      if(g2_peak_time>V['tau']){
        g2_peak_time <- g2_peak_time - V['tau']
      }
      if(g2_peak_time<0){
        g2_peak_time <- g2_peak_time + V['tau']
      }
    }
    res <- row_adder(res, paste0(g1_text, " peak time hours"), peak_time)
    res <- row_adder(res, paste0(g2_text, " peak time hours"), g2_peak_time)
    res <- row_adder(res, "Phase difference estimate", peak_time_diff)
    res <- row_adder(res, "P-value for difference in phase", coefs['phi1', 'p_value'])
  }else{
    res <- row_adder(res, "Shared  peak time hours", peak_time)
  }
  if("phi_decay" %in% names(V)){
    if("phi_decay1" %in% names(V)){
      res <- row_adder(res, paste0(g1_text, " phase decay estimate"), V['phi_decay'])
      res <- row_adder(res, paste0(g2_text, " phase decay estimate"), V['phi_decay']+V['phi_decay1'])
      res <- row_adder(res, "Phase decay difference estimate", V['phi_decay1'])
      res <- row_adder(res, "P-value for phase decay difference", coefs['phi_decay1', 'p_value'])
    }else{
      res <- row_adder(res, "Shared phase decay estimate", V['phi_decay'])
    }
  }

  if("tau1" %in% names(V)){
    res <- row_adder(res, paste0(g1_text, " period estimate"), V['tau'])
    res <- row_adder(res, paste0(g2_text, " period estimate"), V['tau']+V['tau1'])
    res <- row_adder(res, "Period difference estimate (hours)", V['tau1'])
    res <- row_adder(res, "P-value for difference in period", coefs[V['tau1'], 'p_value'])
  }else{
    res <- row_adder(res, "Shared period estimate", V['tau'])
  }
  if("tau_decay" %in% names(V)){
    if("tau_decay1" %in% names(V)){
      res <- row_adder(res, paste0(g1_text, " period decay estimate"), V['tau_decay'])
      res <- row_adder(res, paste0(g2_text, " period decay estimate"), V['tau_decay']+V['tau_decay1'])
    }else{
      res <- row_adder(res, "Shared period decay estimate", V['tau_decay'])
    }
  }

  return(res)
}


