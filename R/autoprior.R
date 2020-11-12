#' Generate brms priors given a formula and data
#'
#' @param formula a formula object, or something that can be coerced to a formula object
#' @param data data frame containing the individual level data
#' @return returns an object of class brmsprior
#'
#' @examples
#'
#' @export
autoprior <- function(formula, data) {
    stopifnot(inherits(formula, "formula"))
    depvar <- all.vars(formula)[1]
    cats <- levels(data[,depvar])
    stopifnot(length(cats) > 0)
    
    dpars <- paste0("mu", cats[-1])
    ### We don't need dpars if we just have two outcomes
    if (length(cats) == 2) {
        dpars <- ""
    }
    
    vars <- all.vars(formula)[-1]

    fct_vars <- brms::get_prior(formula, data)$group
    fct_vars <- unique(fct_vars[fct_vars != ""])

    cont_vars <- brms::get_prior(formula, data)$coef
    cont_vars <- unique(cont_vars[cont_vars != ""])
    cont_vars <- unique(cont_vars[cont_vars != "Intercept"])
    
### "Continuous" variables
    for (d in dpars) {
        if (exists("priors")) {
            priors <- c(priors,
                        brms::set_prior("normal(0, 2.5)",
                                        class = "Intercept",
                                        dpar = d))
        } else {
            priors <- brms::set_prior("normal(0, 2.5)",
                                      class = "Intercept",
                                      dpar = d)
        }
        
        
### Kick off by over-writing the list
        for (v in cont_vars) {
            if (exists("priors")) { 
                priors <- c(priors,
                            brms::set_prior("normal(0, 1)",
                                            class = "b",
                                            coef = v,
                                            dpar = d))
            } else {
                priors <- brms::set_prior("normal(0, 1)",
                                          class = "b",
                                          coef = v,
                                          dpar = d)
            }
        }
        for (v in fct_vars) {
            if (exists("priors")) { 
                priors <- c(priors,
                            brms::set_prior("normal(0, 2.5)",
                                            class = "sd",
                                            group = v,
                                            dpar = d))
            } else {                 
                priors <- brms::set_prior("normal(0, 2.5)",
                                          class = "sd",
                                          group = v,
                                          dpar = d)
            }
        }
    }
    priors
}
