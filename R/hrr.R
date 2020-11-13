#' Estimate hierarchical related regression model
#'
#' @param formula a formula object, or something that can be coerced to a formula object
#' @param data data frame containing the individual level data
#' @param ps data frame containing the post-stratification data; must contain variable `count`
#' @param result data frame containing the results; must contain column names equal to levels of the dependent variable in `formula`
#' @param areavar a character variable containing the name of the variable giving the area
#' @param testing if testing is set to TRUE, brms will use the prior only, no data
#' @param adjust Can aggregate predictions from individual evidence be adjusted?
#' @param dirichlet Use Dirichlet-Multinomial rather than multinomial
#' @param mrp_only Ignore the aggregate element
#' @param weights a character variable containing the name of the weighting variable
#' @param adapt_delta adapt_delta parameter
#' @param max_treedepth max_treedepth parameter
#' @param ... additional parameters passed to cmdstanr
#' @return returns TRUE or reports an error
#'
#' @examples
#'
#' @export
hrr <- function(formula, data, ps, result, areavar,
                testing = FALSE,
                adjust = FALSE,
                dirichlet = FALSE,
                mrp_only = FALSE,
                weights = NULL,
                adapt_delta = 0.9, max_treedepth = 11, ...) {

    ### Input class checking
    if (!inherits(formula, "formula")) {
        formula <- stats::as.formula(formula)
    }
    if (!inherits(data, "data.frame")) {
        stop("data must be a data frame")
    }
    if (!inherits(ps, "data.frame")) {
        stop("data must be a data frame")
    }
    if (!inherits(result, "data.frame")) {
        stop("data must be a data frame")
    }
    if (!inherits(areavar, "character")) {
        stop("areavar must be a character")
    }

### Coerce data frames to data.frame
### (i.e., lose tibbles)
    if (inherits(data, "tbl_df")) data <- as.data.frame(data)
    if (inherits(ps, "tbl_df")) ps <- as.data.frame(ps)
    if (inherits(result, "tbl_df")) result <- as.data.frame(result)
    
### Does everything contain the area var?
    if (!is.element(areavar, all.vars(formula))) {
        stop(paste0("Variable ",
                    areavar,
                    " must be present in formula"))
    }
    if (!is.element(areavar, names(data))) {
        stop(paste0("Variable ",
             areavar,
             " not present in data"))
    }
    if (!is.element(areavar, names(ps))) {
        stop(paste0("Variable ",
                    areavar,
                    " not present in result"))
    }
    if (!is.element(areavar, names(result))) {
        stop(paste0("Variable ",
                    areavar,
                    " not present in result"))
    }

    
### Check the dep. var and its levels
    depvar <- all.vars(formula)[1]
    if (!is.element(depvar, names(data))) {
        stop(paste0("Dependent variable ",
                    depvar,
                    " not present in data"))
    }

### If it's not a factor, coerce it as such
    if (!inherits(data[,depvar], "factor")) {
        data[, depvar] <- factor(data[,depvar])
    }

    cats <- levels(data[,depvar])
    if (length(cats) == 2) {
        family <- "bernoulli"
    } else {
        family <- "categorical"
    }
    
### Check whether the categories have any spaces in them
    if (any(grepl(" ", cats))) {
        stop("Levels of dependent variable cannot have whitespace")
    }

### Check the results file has these columns
    missing_columns <- setdiff(cats, names(result))
    if (length(missing_columns) > 0) {
        stop(paste0("result is missing variables ",
                    paste0(missing_columns, collapse = ", ")))
    }
### Arrange accordingly
    result <- result[,c(areavar, cats)]

### Check whether the ps frame has counts in it
    if (!is.element("count", names(ps))) {
        stop("ps data frame lacks variable count")
    }
    if (any(ps$count < 0)) {
        stop("ps data frame has negative counts")
    }

### Coerce to integer
    ps$count <- as.integer(ps$count)
    
### Check whether the results frame has missing values or negative counts in it
    if (any(is.na(result[,cats]))) {
        stop("results data frame has missing values")
    }
    
    if (any(result[,cats] < 0)) {
        stop("results data frame has negative counts")
    }

### Check the weighting variable is present and all non-negative
    if (!is.null(weights)) {
        if (!is.element(weights, names(data))) {
            stop("weights variable not present in data")
        }
        if (!is.numeric(data[,weights])) {
            stop("weights variable not numeric")
        }
        if (any(data[,weights] < 0)) {
            stop("weights variable has negative entries")
        }
    }
    
### Check whether these data frames are okay to use
    test <- compare_dfs(stats::update(formula, 1 ~ .), data, ps)

    ### Check priors
    my_dots <- dots(...)

    if (!is.element("prior", names(my_dots))) {
        warning("No priors specified. Generating some tight default priors",
                immediate. = TRUE)
        my_dots[["prior"]] <- autoprior(formula, data)
    }

    ### Check threads argument
    if (!is.element("threads_per_chain", names(my_dots))) {
        warning("threads_per_chain not specified. Setting to two threads per chain",
                immediate. = TRUE)
        my_dots[["threads_per_chain"]] <- 2
    }

    
### Generate some preliminary code
    prelim_code <- brms::make_stancode(formula,
                                 data = data,
                                 family = family,
                                 threads = brms::threading(my_dots[["threads_per_chain"]]))

    
    addons <- hrr_code_func(formula = formula,
                            code = prelim_code,
                            data = data,
                            ps = ps,
                            results = result,
                            cats = cats,
                            areavar = areavar,
                            depvar = depvar,
                            adjust = adjust,
                            dirichlet = dirichlet,
                            mrp_only = mrp_only,
                            weights = weights)

### Deal with weights
    new_formula <- add_weights_to_formula(formula, weights)
    
### Got to remove some stuff from my dots
    p <- my_dots[["prior"]]
    my_dots[["prior"]] <- NULL
    thread_count <- my_dots[["threads_per_chain"]]
    my_dots[["threads_per_chain"]] <- NULL

    ### Set grainsize
    grainsize <- floor(nrow(data) / (thread_count * 2))
    message(paste0("Setting grainsize to ", grainsize))
### Permute data to help multithreading (see manual, "To ensure that
    ### chunks (whose size is defined by grainsize) require roughly
    ### the same amount of computing time, we recommend storing
### observations in random order in the data.")
    data <- data[sample(1:nrow(data), size = nrow(data), replace = FALSE),]

    if (testing) {
        code <- brms::make_stancode(formula = new_formula,
                                    data = data,
                                    family = family,
                                    stanvars = addons,
                                    prior = p,
                                    threads = brms::threading(thread_count,
                                                              grainsize = grainsize))
        data <- brms::make_standata(formula = new_formula,
                                    data = data,
                                    family = family,
                                    stanvars = addons,
                                    prior = p,
                                    threads = brms::threading(thread_count,
                                                              grainsize = grainsize))
        retval <- list(code = code,
                       data = data,
                       prior = p,
                       addons = addons)
    } else { 
        brms_args <- c(list(formula = new_formula,
                            data = data,
                            family = family,
                            prior = p,
                            control = list(adapt_delta = adapt_delta,
                                           max_treedepth = max_treedepth),
                            threads = brms::threading(thread_count,
                                                      grainsize = grainsize),
                            backend = "cmdstanr",
                            stanvars = addons),
                       my_dots)
        
        
        if ("algorithm" %in% names(my_dots)) {
            if (my_dots[["algorithm"]] != "sampling") {
                brms_args[["threads"]] <- NULL
                brms_args[["control"]] <- NULL
            }
        }

### What pars do we need?
        ## r_{nvars}_mu{cats}_1
        ## aggmu
        ## adj (if adjust)
        ## prec (if dirichlet)
        ## ps_J_{nvars}_counts
        ## b_mu{cats}_Intercept
        
### we can get this from the second entry in the tparameters block
        tpars_block <- lapply(addons, function(x)x$block)
        tpars_block <- which(unlist(tpars_block) == "tparameters")[2]
        tpars_block <- addons[[tpars_block]]$scode

        rpars <- stringr::str_extract_all(tpars_block,
                                          pattern = "(r_[^,]*)")[[1]]
        bpars <- stringr::str_extract_all(tpars_block,
                                          pattern = "(b_[^,]*)")[[1]]
        intpars <- paste0(bpars, "_Intercept")

        nvars <- sub("r_", "", rpars)
        nvars <- sub("_.*", "", nvars)
        nvars <- as.numeric(nvars)
        count_vars <- glue::glue("ps_J_{nvars}_counts")

        my_save_pars <- c(rpars, bpars, intpars,
                       count_vars,
                       "aggmu")
        if (adjust) my_save_pars <- c(my_save_pars, "adj")
        if (dirichlet) my_save_pars <- c(my_save_pars, "prec")

        brms_args[["save_pars"]] <- save_pars(group = FALSE,
                                              manual = my_save_pars)
        
### Execute the call
        retval <- do.call(brms::brm,
                          args = brms_args)

### Add on hrr specific elements
        retval$areavar <- areavar
        class(retval) <- c("hrrfit", class(retval))
        
    }
    
    return(retval)

}


hrr_code_func <- function(formula, code, data, ps, results,
                          cats, areavar, depvar, adjust, dirichlet, mrp_only, weights) {
### Purpose: pull together all the stanvars
### Input: code and data
### Output: stanvars

    retval <- pll_to_pred_prob(code, adjust)

    if (dirichlet) {
        retval <- retval + dm_func()
    }
    ## add_ps_data_code(code) +
    retval <- retval +
        add_ps_tdata_code(levels(factor(data[,depvar]))) +
        hrr_data_func(formula, data = data, ps = ps,
                      results = results, cats = cats,
                      areavar = areavar, depvar = depvar)
    
    if (adjust|dirichlet) {
        retval <- retval + 
            add_pars_code(code, adjust, dirichlet)
    }
    
    retval <- retval +
        add_tpars_code(code, adjust, dirichlet)

    if (!mrp_only) {
        retval <- retval +
            add_model_code(code, adjust, dirichlet, data, results)
    }
    
    retval <- retval + 
        add_genquant_code(code, formula, data, adjust, dirichlet)
    
    retval
}

dm_func <- function() {
    scode <- "
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
  
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y))) 
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
"
    brms::stanvar(scode = scode,
                  block = "functions")
                  
}

add_ps_tdata_code <- function(cats) {
### Purpose: add data transformations to the transformed data block
### Input: list of output levels
### Output: stanvar

    cats <- cats[-1]
    if (any(grepl(" ", cats))) {
        stop("Don't yet know how to handle outcome levels with whitespace")
    }

    if(length(cats) > 1) {
        family <- "categorical"
    } else {
        family <- "bernoulli"
    }
    
    if (family == "categorical") { 
        pscode_tdata_top <- glue::glue('
  matrix[ps_N, K_mu<<cats>> - 1] ps_Xc_mu<<cats>>;  // 
', .open = "<<", .close = ">>")
    
        pscode_tdata_bottom <- glue::glue('
  for (i in 2:K_mu<<cats>>) {
    ps_Xc_mu<<cats>>[, i - 1] = ps_X_mu<<cats>>[, i] - means_X_mu<<cats>>[i - 1];
  } 
', .open = "<<", .close = ">>")

    } else {
        pscode_tdata_top <- '
  matrix[ps_N, K - 1] ps_Xc;  // 
'
        
        pscode_tdata_bottom <- '
  for (i in 2:K) {
    ps_Xc[, i - 1] = ps_X[, i] - means_X[i - 1];
  } 
'
        
    }

    brms::stanvar(pscode_tdata_top,
            scode = pscode_tdata_top,
            block = "tdata",
            position = "start") +
    brms::stanvar(pscode_tdata_top,
            scode = pscode_tdata_bottom,
            block = "tdata",
            position = "end")

}

pll_to_pred_prob <- function(code, adjust) {
### Input: preliminary code
### output: stanvar
    ### This might not work with weights
    ## require(stringr)
    ## require(brms)
    ## Get the function code
    if (grepl("categorical_logit_lpmf", code)) {
        family <- "categorical"    
    } else {
        family <- "bernoulli"
    }
        
    code <- stringr::str_extract(code,
                                 pattern = stringr::regex("real partial_log_lik.*return ptarget;\\s+\\}",
                                                          dotall = TRUE,
                                                          multiline = TRUE))
    
### Amend the function declaration: name and return type
    if (adjust) {
        if (family == "categorical") { 
            code <- stringr::str_replace(code,
                                         stringr::fixed("real partial_log_lik(int[] seq, "),
                                         "vector pred_prob(vector adj, ")
        } else {
            code <- stringr::str_replace(code,
                                         stringr::fixed("real partial_log_lik(int[] seq, "),
                                         "real pred_prob(real adj, ")
        }
        
    } else {
        if (family == "categorical") { 
            
            code <- stringr::str_replace(code,
                                         stringr::fixed("real partial_log_lik(int[] seq, "),
                                         "vector pred_prob(")
        } else {
            code <- stringr::str_replace(code,
                                         stringr::fixed("real partial_log_lik(int[] seq, "),
                                         "real pred_prob(")
        }
        
    }
        
### Amend the function declaration: additional argument
    if (family == "categorical") { 
        code <- stringr::str_replace(code,
                                     "vector pred_prob(.*?)\\) \\{",
                                     "vector pred_prob\\1, int[] psweights) {")
    } else {
        code <- stringr::str_replace(code,
                                     "real pred_prob(.*?)\\) \\{",
                                     "real pred_prob\\1, int[] psweights) {")
    }
    
### Replace the current declaration of the return variable
    if (family == "categorical") { 
        code <- stringr::str_replace(code,
                                     "real ptarget = 0;",
                                     "vector[ncat] pp = rep_vector(0, ncat);")
    } else {
        code <- stringr::str_replace(code,
                                     "real ptarget = 0;",
                                     "real pp = 0;")
    }
    

### Amend the calculation
    oldcode <- code
    if (family == "categorical") { 
        if (adjust) {
            code <- stringr::str_replace(code,
                                         stringr::fixed("ptarget += categorical_logit_lpmf(Y[nn] | mu[n]);"),
                                         "pp += softmax(adj + mu[n]) * psweights[nn];")
        } else {
            code <- stringr::str_replace(code,
                                         stringr::fixed("ptarget += categorical_logit_lpmf(Y[nn] | mu[n]);"),
                                         "pp += softmax(mu[n]) * psweights[nn];")
        }
    } else {
        if (adjust) {
            code <- stringr::str_replace(code,
                                         stringr::fixed("ptarget += bernoulli_logit_lpmf(Y[nn] | mu[n]);"),
                                         "for (n in 1:N) {
      int nn = n + start - 1;
      pp += inv_logit(adj + mu[n] + Xc[start:end] * b) * psweights[nn];
}")
        } else {
            code <- stringr::str_replace(code,
                                         stringr::fixed("ptarget += bernoulli_logit_glm_lpmf(Y[start:end] | Xc[start:end], mu, b);"),
                                         "for (n in 1:N) {
      int nn = n + start - 1;
      real mu0 = Xc[nn] * b;
      real unweighted = inv_logit(mu[n] + mu0);
      pp += unweighted * psweights[nn] / sum(psweights[start:end]);
}")
        }
    }
    
    if (code == oldcode) {
        stop("Terribly sorry, but some code we expected to find hasn't been found, and so can't be replaced")
    }
    
### Amend the return value
    if (family == "categorical") { 
        code <- stringr::str_replace(code,
                                     stringr::fixed("return ptarget;"),
                                     "return (pp / sum(pp));")
    } else {
        code <- stringr::str_replace(code,
                                     stringr::fixed("return ptarget;"),
                                     "return pp;")
    }

    brms::stanvar(scode = code,
            block = "functions")
} 

add_pars_code <- function(code, adjust, dirichlet) {
    if (grepl("categorical_logit_lpmf", code)) {
        family <- "categorical"    
    } else {
        family <- "bernoulli"
    }
    code <- ""
    if (adjust) {
        code <- "  vector [ncat-1] adj0; "

    }
    if (dirichlet) {
        if (family == "categorical") { 
        code <- paste0(code, "
  real <lower=0>prec;
")
        } else {
        code <- paste0(code, "
  real <lower=0>kappa;
")
        }
    }
    
    code <- brms::stanvar(scode = code,
                          position = "start",
                          block = "parameters")
    return(code)
}

add_tpars_code <- function(code, adjust, dirichlet) {
### Purpose: add transformed parameters declaration
### Input: stan code
### Output: stanvar object
    if (grepl("categorical_logit_lpmf", code)) {
        family <- "categorical"    
    } else {
        family <- "bernoulli"
    }

    if (family == "categorical") {
        top <- brms::stanvar(scode = "matrix[nAreas, ncat] aggmu;",
                   block = "tparameters",
                   position = "start")
        if (adjust) { 
            top <- top +
                brms::stanvar(scode = "vector [ncat] adj; ",
                              block = "tparameters",
                              position = "start")
        }
    } else {
        top <- brms::stanvar(scode = "vector[nAreas] aggmu;",
                             block = "tparameters",
                             position = "start")
        if (adjust) {
            top <- top +
                brms::stanvar(scode = "real adj; ",
                              block = "tparameters",
                              position = "start")
        }
        
    }
    

    old_func <- stringr::str_extract(code,
                            pattern = "reduce_sum\\(partial_log_lik.*")

### Amend the beginning of the function call
    if (adjust) {
        new_func <- stringr::str_replace(old_func,
                                         stringr::fixed("reduce_sum(partial_log_lik, seq, grainsize, "),
                                         "pred_prob(adj, areastart[i], areastop[i], ")
    } else { 
        new_func <- stringr::str_replace(old_func,
                                         stringr::fixed("reduce_sum(partial_log_lik, seq, grainsize, "),
                                         "pred_prob(areastart[i], areastop[i], ")
    }
    
### Amend the variable calls in the middle

    new_func <- stringr::str_replace_all(new_func,
                                stringr::fixed("Xc_"),
                                "ps_Xc_")

    new_func <- stringr::str_replace_all(new_func,
                                stringr::fixed("Xc,"),
                                "ps_Xc,")
    
    new_func <- stringr::str_replace_all(new_func,
                                stringr::fixed("J_"),
                                "ps_J_")
    
    new_func <- stringr::str_replace_all(new_func,
                                stringr::fixed("Z_"),
                                "ps_Z_")

### Amend the end by adding psweights arg);
### only transpose if we need to
    if (family == "categorical") { 
        new_func <- stringr::str_replace(new_func,
                                         stringr::fixed(");"),
                                         ", ps_counts)';")
    } else {
        new_func <- stringr::str_replace(new_func,
                        stringr::fixed(");"),
                        ", ps_counts);")
    }
    
    bottom <- paste0("
    for (i in 1:nAreas) {
        aggmu[i] = ",
    new_func,
    "
        }")

    bottom <- brms::stanvar(scode = bottom,
                      block = "tparameters",
                      position = "end")
    if (adjust) {
        if (family == "categorical") { 
        bottom <- 
            brms::stanvar(scode = "adj = append_row(0, adj0); ",
                          block = "tparameters",
                          position = "end") +
            bottom
        } else {
            bottom <- 
                brms::stanvar(scode = "adj = adj0; ",
                              block = "tparameters",
                              position = "end") +
                bottom
        }
        
    }
    
    if (family == "bernoulli" & dirichlet) {
        top <- brms::stanvar(scode = "real prec;",
                             block = "tparameters",
                             position = "start") +
            top
        bottom <- brms::stanvar(scode = "prec = 1 / kappa;",
                                block = "tparameters",
                                position = "end") +
            bottom
        
    }
    

    top + bottom
    
}

add_model_code <- function(code, adjust, dirichlet, data, results) {
### Purpose: add a sampling statement
### Input: preliminary code
### Output: stanvar
    if (grepl("categorical_logit_lpmf", code)) {
        family <- "categorical"    
    } else {
        family <- "bernoulli"
    }

    if (family == "categorical") { 
        if (dirichlet) {
            scode <- "
  if (!prior_only) {
   for (i in 1:nAreas) { 
      aggy[i] ~ dirichlet_multinomial(prec * to_vector(aggmu[i]));
   }
  }
"
        } else { 
            scode <- "
  if (!prior_only) {
   for (i in 1:nAreas) { 
      aggy[i] ~ multinomial(to_vector(aggmu[i]));
   }
  }
"
        }
    } else {
            if (dirichlet) {
            scode <- "
  if (!prior_only) {
   for (i in 1:nAreas) {
      real share;
      // multiplying by 1.0 nec to promote to real from int
      share = (1.0 * aggy[i, 2]) / (1.0 * (aggy[i,1] + aggy[i, 2]));
      share ~ beta_proportion(aggmu[i], prec);
   }
  }
"
        } else { 
            scode <- "
  if (!prior_only) {
   for (i in 1:nAreas) {
      int totarea = (aggy[i, 1] + aggy[i, 2]);
      aggy[i, 2] ~ binomial(totarea, aggmu[i]);
   }
  }
"
        }
    }
    
    
    

    if (adjust) {
        scode <- paste0(scode, "
   target += normal_lpdf(adj0 | 0, 1);
", collapse = "\n")                        
    }
            
    if (dirichlet) {
        if (family == "categorical") { 
        ## The median on the normal scale should equal nrow(dat) / nAreas
        ## and we should assign zero probability to a pseudo count greater than the observed total in any area
        log_m <- log(nrow(data) / nrow(results))
        ## Upper bound
        ub <- mean(rowSums(results[,-1]))
        ## Possible values for SD on the log scale
        inseq <- seq(1, 4, length.out = 100)
        ## Pseudo counts out
        pseudo_counts <- stats::qlnorm(p = 0.99,
                                       meanlog = log_m,
                                       sdlog = inseq)
        log_sd <- max(inseq[which(pseudo_counts < ub)])
        scode <- paste0(scode, "
// have a sensible default prior for the precision parameter here
   target += lognormal_lpdf(prec | ",
log_m,
", ",
log_sd,
");",
collapse = "\n")                        
        } else { ## binomial
            shares <- results[,3] / rowSums(results[,-1])
            share_sd <- sd(shares)
                    scode <- paste0(scode, "
// have a sensible default prior for the precision parameter here
   target += normal_lpdf(kappa | 0, ",
share_sd,
");",
collapse = "\n")                        
        }
    }
    
    
    retval <- brms::stanvar(scode = scode,
                            block = "model")

}

hrr_data_func <- function(formula, data, ps, results, cats, areavar, depvar) {
### Purpose: construct all the additional data as a stanvar
### Input: whole bunch of things
### Output: stanvar

    if (length(cats) > 2) {
        family <- "categorical"
    } else {
        family <- "bernoulli"
    }
    
    ### Do I need to change the threading argument?
    prelim_data <- brms::make_standata(formula,
                          data,
                          family = family,
                          threads = brms::threading(4))

### Amend the ps data so that it includes the depvar
    ps[,depvar] <- sample(data[,depvar],
                          size = nrow(ps),
                          replace = TRUE)

### Make sure that the ps data is ordered by the geog. var
    ps <- ps[order(ps[,areavar]),]
### And so too are the results
    results <- results[order(results[,areavar]),]
    
    ps_data <- brms::make_standata(formula,
                                   ps,
                                   family = family,
                                   threads = brms::threading(4))
    

### Retain, from ps_data, anything which begins with an X_, a Z_, or or a J_
    new_ps_data <- ps_data[c(grep("X_|Z_|J_", names(ps_data)),
                             which(names(ps_data) == "X"))]
### Rename these variables
    names(new_ps_data) <- paste0("ps_", names(new_ps_data))
    new_ps_data$ps_N <- ps_data$N
    new_ps_data$nAreas <- length(unique(ps[, areavar]))

    ps$rowno <- 1:nrow(ps)
    
    areapos_start <- tapply(ps$rowno, ps[,areavar], min)
    areapos_stop <- tapply(ps$rowno, ps[,areavar], max)

    new_ps_data$areastart <- areapos_start
    new_ps_data$areastop <- areapos_stop
    new_ps_data$ps_counts <- ps$count

### create a stanvar from this?
    i <- 1
    data_block <- brms::stanvar(new_ps_data[[i]],
                                name = names(new_ps_data)[i],
                                block = "data")
    
    for (i in 2:length(new_ps_data)) {
        data_block <- c(data_block,
                        brms::stanvar(new_ps_data[[i]],
                                      name = names(new_ps_data)[i],
                                      block = "data"))
    }

    ### Need to add ncat if it's binomia;
    if (family == "bernoulli") {
        ncat <- 2
        data_block <- c(data_block,
                        brms::stanvar(ncat,
                                      scode = "int ncat;",
                                      block = "data"))

    }
    
### Combine with the other data
    
    aggy <- as.matrix(results[,cats])
    
    data_block <- c(data_block,
                    brms::stanvar(aggy,
                                  scode = "int aggy[nAreas, ncat];",
                                  block = "data"))


    
### Get aggregate results out
### Make sure they follow the same order as the ps data
    return(data_block)
}

dots <- function(...) {
  eval(substitute(alist(...)))
}


add_genquant_code <- function(code,
                              formula,
                              data,
                              adjust,
                              dirichlet) {
    if (grepl("categorical_logit_lpmf", code)) {
        family <- "categorical"    
    } else {
        family <- "bernoulli"
    }
    ## Declarations
    ## Big ps count objects
    new_code <- "
  int psw_counts[ps_N, ncat];

"
### Get the number of groups considered
    tmpcode <- brms::make_stancode(formula,
                                   data = data,
                                   family = "categorical")
    grp_indicators <- grep("int<lower=1> J_",
                           stringr::str_split(tmpcode, pattern = "\n")[[1]],
                           value = TRUE)
    nvars <- unlist(stringr::str_extract_all(grp_indicators, "J_([0-9]+)"))
    nvars <- as.numeric(sub("J_", "", nvars))
    
    ## Now objects for all the cats.
    ## Need to get the right number here
    new_code <- stringr::str_c(new_code,
                           stringr::str_c(glue::glue("
      int ps_J_<<nvars>>_counts[N_<<nvars>>, ncat];
", .open = "<<", .close = ">>"), collapse = "\n"))

### Initialize the psw+counts
    new_code <- stringr::str_c(new_code,
                   "
  for (p in 1:ps_N) {
    for (k in 1:ncat) {
      psw_counts[p, k] = 0;
    }
  }
", collapse = "\n")
    
    ### Initialize each
    new_code <- stringr::str_c(new_code,
                   stringr::str_c(glue::glue("
  for (p in 1:ps_N) {
    for (k in 1:ncat) {
      ps_J_<<nvars>>_counts[ps_J_<<nvars>>[p],k] = 0;
    }
  }

", .open = "<<", .close = ">>"), collapse = "\n"), collapse = "\n")

    tpars_code <- add_tpars_code(code, adjust, dirichlet)
### Extract that part of the tpars_code which involves pred_probs
    pos <- lapply(tpars_code, function(x)grepl("pred_prob", x$scode))
    pos <- unlist(pos)
    pos <- which(pos)
    
    tpars_code <- tpars_code[[pos]]$scode
    
    tpars_code <- stringr::str_replace(tpars_code,
                                       "nAreas",
                                       "ps_N")

    if (family == "categorical") { 
        tpars_code <- stringr::str_replace(tpars_code,
                                           stringr::fixed("aggmu[i] = "),
                                           "row_vector [ncat] tmp = ")
    } else {
        tpars_code <- stringr::str_replace(tpars_code,
                                           stringr::fixed("aggmu[i] = "),
                                           "real tmp = ")
    }
    
    
    tpars_code <- stringr::str_replace(tpars_code,
                                       stringr::fixed("areastart[i], areastop[i]"),
                                       "i, i")

    if (family == "categorical") { 
    tpars_code <- stringr::str_replace(tpars_code,
                                           stringr::fixed("}"),
                                           "
psw_counts[i] = multinomial_rng(to_vector(tmp), ps_counts[i]);
}
")
    } else {
            tpars_code <- stringr::str_replace(tpars_code,
                                           stringr::fixed("}"),
                                           "
psw_counts[i, 2] = binomial_rng(ps_counts[i], tmp);
psw_counts[i, 1] = ps_counts[i] - psw_counts[i, 2];
}
")
        }
   

### Add this on
    new_code <- stringr::str_c(new_code,
                           tpars_code,
                           collapse = "\n")
    
### Then for each

    new_code <- stringr::str_c(new_code,
                           stringr::str_c(glue::glue("
  for (p in 1:ps_N) {
    for (k in 1:ncat) {
      ps_J_<<nvars>>_counts[ps_J_<<nvars>>[p],k] = ps_J_<<nvars>>_counts[ps_J_<<nvars>>[p],k] +
                               psw_counts[p, k];
    }
  }

", .open = "<<", .close = ">>"), collapse = "\n"), collapse = "\n")

### Return as a stanvar
    brms::stanvar(scode = new_code,
                  block = "genquant")

}


add_weights_to_formula <- function(formula, weights) {

    if (is.null(weights)) {
        return(formula)
    } else {
        f <- as.character(formula)
        f[2] <- paste0(f[2], "|weights(", weights, ")")
        return(brms:::bf(paste(f[-1], collapse = "~")))
    }
}

    
