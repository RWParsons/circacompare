utils::globalVariables(c('time', 'measure', 'group', 'eq'))

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


create_formula <- function(main_params=c("k", "alpha", "phi"), decay_params=c(), grouped_params=c()){
  if(length(grouped_params[!grouped_params %in% main_params])>0){
    stop("All grouped parameters must be within the main parameters.")
  }
  if(all(c("phi", "tau") %in% grouped_params)){
    stop("Either phase or period can have grouped effects, but not both.")
  }
  build_component <- function(pattern, eq=F){
    t <- ifelse(eq, "time", "time_r")

    main <- main_params[grep(pattern, main_params)]
    decay <- decay_params[grep(pattern, decay_params)]
    grouped <- grouped_params[grep(pattern, grouped_params)]
    if(length(main)<length(decay)|length(main)<length(grouped)){
      stop(paste0("You can't model the decay of a parameter which doesn't exist in the model!\nAdd `",
                  gsub("^.", "", pattern), "` to the `main_params`"))
    }
    if(eq & length(main)!=0){
      main <- paste0("V['", main, "']")
    }
    if(length(decay)!=0){
      decay <- paste0(decay, "_decay")
      if(eq){
        decay <- paste0("V['", decay, "']")
      }
      main <- paste0(main, "*exp(-(", decay, ")*", t, ")")
    }
    component <- main
    group_match <- grouped
    group_decay <- group_match[grepl("decay", group_match)]

    if(length(group_decay)!=0){
      component <- gsub(group_decay, paste0(group_decay, '+', group_decay, "1"), component)
    }
    group_main <- group_match[!grepl("decay", group_match)]
    if(length(group_main)!=0){
      if(eq){
        group_main <- paste0("V['", group_main, "1']")
      }else{
        group_main <- paste0(group_main, "1")
      }
      component <- paste0(component, "+", group_main)
    }

    if(pattern=="^tau"){
      if(length(component)==0){
        if(eq){
          component <- "(2*pi/period)*"
        }else{
          component <- ""
        }
      }else{
        if(eq){
          component <- paste0("(24/(", component, "))*")
        }else{
          component <- paste0("24/(", component, ")*")
        }
      }
    }
    return(component)
  }
  res_formula <- paste0("measure~", build_component("^k"), "+(", build_component("^alpha"),
                        ")*cos(", build_component("^tau"), "(time_r-(", build_component("^phi"), ")))")
  res_formula <- gsub("1", "1*x_group", res_formula)
  res_formula <- stats::as.formula(res_formula)
  res_equation <- paste0(build_component("^k", eq=T), "+(", build_component("^alpha", eq=T),
                         ")*cos(", build_component("^tau", eq=T), "(time-(", build_component("^phi", eq=T), ")))")

  if(length(grouped_params)>0){
    res_equation <- list(
      g1=paste0("eq1 <- function(time) {", gsub("\\+V\\['[a-z]*1'\\]", "", res_equation), "}"),
      g2=paste0("eq2 <- function(time) {", res_equation, "}")
    )
  }else{
    res_equation <- paste0("eq <- function(time) {time<-(time/24)*2*pi;return(", res_equation, ")}")
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

