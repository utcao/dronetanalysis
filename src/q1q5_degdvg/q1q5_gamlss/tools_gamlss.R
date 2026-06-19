fit_one <- function(index, dat, meta,
                  mu_formula, sigma_formula, offsets, sigma_start=0.1, mu_sum_vars = NULL) {
  # Build the model.matrix call for mu. If `mu_sum_vars` is provided,
  # apply sum-to-zero coding (contr.sum) to those covariates only.
  if(!is.null(mu_sum_vars) && length(mu_sum_vars) > 0){
    contrasts_text <- paste0("contrasts.arg = list(",
                             paste(paste0(mu_sum_vars, " = contr.sum"), collapse = ", "),
                             ")")
    mm_call <- glue("model.matrix({mu_formula}, meta, {contrasts_text})")
  } else {
    mm_call <- glue("model.matrix({mu_formula}, meta)")
  }
  formula_text <- mm_call

  if(is.null(sigma_formula)){
    sigma_design <- as.formula("~1")
  }else{
    sigma_design_text <- glue("~ model.matrix({sigma_formula}, meta)")
    sigma_design <- as.formula(sigma_design_text)
    formula_text_sig <- glue("model.matrix({sigma_formula}, meta)")
  }

  mu_design <- as.formula(glue("dat[{index},] ~ {mm_call} + offset(offsets)"))

  model <- tryCatch({
    gamlss(mu_design,                   # drop intercept – GAMLSS adds it
        sigma.formula = sigma_design,
        family = NBI(),                        # negative-binomial with NBI par.
        sigma.start = sigma_start,
        control = gamlss.control(n.cyc = 100, trace = FALSE))
      }, error = function(e) {
          futile.logger::flog.error(
              glue("Error in fitting model for gene {index} {rownames(dat)[{index}]}: {e$message}"))
          return(NULL)
  })
  if(!is.null(model)){
      if( ! is.null(sigma_formula)){
      names(model$sigma.coefficients) <- map_chr( names(model$sigma.coefficients),
                                      str_remove, pattern = fixed(formula_text_sig))
  }
  names(model$mu.coefficients) <- map_chr( names(model$mu.coefficients),
                                      str_remove, pattern = fixed(formula_text))
  return(model)
  }
}

fit_sum <- function(model) {
  if(is.null(model)){
    model_sum <- NULL
  }else{
    model_sum <- invisible(summary(model))
  }
  return(model_sum)
}