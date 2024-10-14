make_stan_code <- function(f, data, ps, aux, res, adjust, overdispersed, threading, mrp_only) {

    function_code <- make_function_code(f, data, ps, aux, adjust)
    data_code <- make_data_code(f, data, ps, aux)
    tdata_code <- make_tdata_code(f, data, ps, aux)
    params_code <- make_params_code(f, data, ps, aux, adjust, overdispersed, mrp_only)
    tparams_code <- make_tparams_code(f, data, ps, aux, adjust, overdispersed, mrp_only)
    model_code <- make_model_code(f, data, ps, aux, res, adjust, overdispersed, threading, mrp_only)
    genquant_code <- make_genquant_code(f, data, ps, aux, adjust)

    ## Concatenate all the code blocks
    paste(function_code,
          data_code,
          tdata_code,
          params_code,
          tparams_code,
          model_code,
          genquant_code,
          sep = "\n"
          )
    
}

make_function_code <- function(f, data, ps, aux, adjust) {
    require(Formula)
    f <- Formula(f)
    
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)
    ncatvars <- length(ind_predictors)
    
### For each of the levels of the categorical outcome variable
    dv <- terms(f, lhs = TRUE, rhs = c(FALSE, FALSE))
    dv <- all.vars(dv)

    dv_levels <- levels(factor(data[,dv]))

### Need to write a function to generate predicted probabilities
### and a function for reduce_sum (different output signatures)
    code <- "functions {\n"

    code <- paste0(code,
                   "  array[] int sequence(int start, int end) { \n  array[end - start + 1] int seq;\n  for (n in 1:num_elements(seq)) {\n seq[n] = n + start - 1; \n   }\n    return seq; \n  } \n")
    
    code <- paste0(code,
                   "\n\n")

    code <- paste0(code,
                   "
  real dirichlet_multinomial_lpmf(array[] int y, vector alpha) {
    real alpha_plus = sum(alpha);
  
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y))) 
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

"
)
    
### (1) Predicted probabilities
### (a) function beginning
    code <- paste0(code,
                   " vector pred_prob(")

    if (adjust) {
        code <- paste0(code,
                       "vector adj, ")
    }
    
    code <- paste0(code,
                   "int start, int end, int ncat, array[] int Y, matrix Xc, \n")

### (b) Alternative specific parameters
    for(d in dv_levels[-1]) {
        addon <- paste0(" vector b_mu", d, ", ")
        code <- paste0(code,
                       addon)
        addon <- paste0(" real Intercept_mu", d, ", ")
        code <- paste0(code,
                       addon)
        for (k in 1:(ncatvars+1)) {
            addon <- paste0("vector r_", k, "_", d, ", ")
            code <- paste0(code,
                           addon)
        }
        
    }

### (c) Data
    for (k in 1:(ncatvars+1)) {
        addon <- paste0("array[] int J_", k, ", ")
        code <- paste0(code, addon)
    }

### Close data block with weights
    code <- paste0(code,
                   " array[] int weights) {")

    code <- paste0(code,
                   "
    vector[ncat] pp = rep_vector(0, ncat);
    int N = end - start + 1;
")
    
    for (d in dv_levels[-1]) {
        addon <- paste0("    vector[N] mu",
                        d,
                        " = Intercept_mu",
                        d,
                        " + Xc[start:end] * b_mu",
                        d,
                        ";\n")
        code <- paste0(code,
                       addon)
    }

    code <- paste0(code,
                   "    // linear predictor matrix
    array[N] vector[ncat] mu;\n")

    code <- paste0(code,
                   "      for (n in 1:N) {\n")
    
    code <- paste0(code,
                   "      int nn = n + start - 1;\n")


    for (d in dv_levels[-1]) {
        addon <- paste0("mu", d, "[n] += ")
        code <- paste0(code, addon)
        for (k in 1:(ncatvars+1)) {
            code <- paste0(code,
                           paste0("r_", k, "_", d, "[J_", k, "[nn]]"))
            if (k < (ncatvars + 1)) {
                code <- paste0(code,
                               " + ")
            } else {
                code <- paste0(code,
                               "; \n")
            }
        }
    }

    code <- paste0(code,
                   "\n}\n")

    code <- paste0(code,
                   paste0("
        for (n in 1:N) {
      mu[n] = transpose([0, "))
    for (d in dv_levels[-1]) {
        addon <- paste0("mu", d, "[n]")
        if (d == dv_levels[length(dv_levels)]) {
            addon <- paste0(addon, "]);")
        } else {
            addon <- paste0(addon, ", ")
        }
        code <- paste0(code, addon)
    }
        
    code <- paste0(code,
                   "\n}\n")

    code <- paste0(code,
                   "
    for (n in 1:N) {
      int nn = n + start - 1;
      pp += softmax(")
    if (adjust) {
        code <- paste0(code,
                       "adj + ")
    }
    code <- paste0(code,
                   "mu[n]) * weights[nn];
    }
    return (pp / sum(pp));
}
")

    ### (2) Reduce sum
    code <- paste0(code,
                   "  real partial_log_lik(array[] int seq, int start, int end, int ncat, array[] int Y, matrix Xc,  ")

    ### Copy code from above
### (b) Alternative specific parameters
    for(d in dv_levels[-1]) {
        addon <- paste0(" vector b_mu", d, ", ")
        code <- paste0(code,
                       addon)
        addon <- paste0(" real Intercept_mu", d, ", ")
        code <- paste0(code,
                       addon)
        for (k in 1:(ncatvars+1)) {
            addon <- paste0("vector r_", k, "_", d, ", ")
            code <- paste0(code,
                           addon)
        }
        
    }

### (c) Data
    for (k in 1:(ncatvars+1)) {
        addon <- paste0("array[] int J_", k)
        if (k < (ncatvars+1)) {
            addon <- paste0(addon, ", ")
        }
        code <- paste0(code, addon)
    }

### Close data block with weights
    code <- paste0(code,
                   ") {")

    code <- paste0(code,
                   "
    real ptarget = 0;
    int N = end - start + 1;
")
    
    for (d in dv_levels[-1]) {
        addon <- paste0("    vector[N] mu",
                        d,
                        " = Intercept_mu",
                        d,
                        " + Xc[start:end] * b_mu",
                        d,
                        ";\n")
        code <- paste0(code,
                       addon)
    }

    code <- paste0(code,
                   "    // linear predictor matrix
    array[N] vector[ncat] mu;\n")

    code <- paste0(code,
                   "      for (n in 1:N) {\n")
    
    code <- paste0(code,
                   "      int nn = n + start - 1;\n")


    for (d in dv_levels[-1]) {
        addon <- paste0("mu", d, "[n] += ")
        code <- paste0(code, addon)
        for (k in 1:(ncatvars+1)) {
            code <- paste0(code,
                           paste0("r_", k, "_", d, "[J_", k, "[nn]]"))

            if (k < (ncatvars+1)) {
                code <- paste0(code,
                               " + ")
            } else {
                code <- paste0(code,
                               "; \n")
            }
        }
    }

    code <- paste0(code,
                   "\n}\n")

    code <- paste0(code,
                   paste0("
        for (n in 1:N) {
      mu[n] = transpose([0, "))
    for (d in dv_levels[-1]) {
        addon <- paste0("mu", d, "[n]")
        if (d == dv_levels[length(dv_levels)]) {
            addon <- paste0(addon, "]);")
        } else {
            addon <- paste0(addon, ", ")
        }
        code <- paste0(code, addon)
    }
        
    code <- paste0(code,
                   "\n}\n")

    
    code <- paste0(code,
                   "
    for (n in 1:N) {
      int nn = n + start - 1;
      ptarget += categorical_logit_lpmf(Y[nn] | mu[n]);
    }
    return ptarget;
}
")

    ### Finish code block
    code <- paste0(code,
                   "\n}\n\n")

    return(code)
}

make_data_code <- function(f, data, ps, aux) {

    require(Formula)
    f <- Formula(f)

    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)

    code <- "data {\n"

    code <- paste0(code,
                   " int<lower=1> N; // total number of observations\n")

    code <- paste0(code,
                   " int nAreas; \n")
    
### Dependent variable
    depvar_type <- get_depvar_type(f, data)

    if (depvar_type == "cont") {
                                        # do nothing
        stop("Continuous variables not yet supported")
    } else if (depvar_type == "bin") {
        ## do nothing
        ## stop("Binary variables not yet supported")

        ## Treat if as though it were a categorical variable
        code <- paste0(code,
                       " int<lower=2> ncat;  // number of categories\n")
        code <- paste0(code,
                       " int<lower=1, upper=ncat> Y[N];  // response variable\n")
    } else {
        code <- paste0(code,
                       " int<lower=2> ncat;  // number of categories\n")
        code <- paste0(code,
                       " int<lower=1, upper=ncat> Y[N];  // response variable\n")
    }

### auxiliary predictors
    code <- paste0(code,
                   " int<lower=1> K_X;  // number of population-level effects\n")

    code <- paste0(code,
                   " matrix[N, K_X] X;  // population-level design matrix\n")

    code <- paste0(code,
                   "  int grainsize;  // grainsize for threading\n")

### For each categorical predictor
    ncatvars <- length(ind_predictors)
    for(i in 1:(ncatvars+1)) {
        addon <- paste0("  int<lower=1> N_", i,
                        ";  // number of grouping levels\n",
                        " int<lower=1, upper = N_", i,
                        "> J_", i,
                        "[N];  // grouping indicator per observation\n")
        code <- paste0(code,
                       addon)
    }

    code <- paste0(code,
                   " int prior_only;  // should the likelihood be ignored?\n")

### Now post-stratification stuff
    code <- paste0(code,
                   " int<lower=1> ps_N; \n")

    code <- paste0(code,
                   " matrix[ps_N, K_X] ps_X;  // \n")    

    for(i in 1:(ncatvars+1)) {
        addon <- paste0(" int<lower=1, upper=N_", i, "> ps_J_", i,
                        "[ps_N];  // grouping indicator per observation\n")
        code <- paste0(code,
                       addon)
    }

    code <- paste0(code,
                   "\n int<lower=1,upper=ps_N> areastart[nAreas];\n int<lower=2,upper=ps_N> areastop[nAreas];\n")

    code <- paste0(code,
                   " int<lower=1, upper=nAreas> ps_area[ps_N];\n")
    
    code <- paste0(code,
                   " int ps_counts[ps_N];\n")

    if (depvar_type == "cat" || depvar_type == "bin") { 
        code <- paste0(code,
                       " int aggy[nAreas, ncat];\n")
    } else {
        stop("Only categorical variables supported")
    }

    ### Finish code block
    code <- paste0(code,
                   "\n}\n\n")
    return(code)
}


get_depvar_type <- function(f, data) {
    require(Formula)
    f <- Formula(f)

    ## Does it contain all the variables in the formula?
    dv <- terms(f, lhs = TRUE, rhs = c(FALSE, FALSE))
    dv <- all.vars(dv)

    dv_class <- class(data[,dv])

    luniq_dv <- length(unique(data[,dv]))
    if (dv_class == "numeric") {
        return("cont")
    } else if (luniq_dv == 2) {
        return("bin")
    } else {
        return("cat")
    }
    
}

make_tdata_code <- function(f, data, ps, aux) {
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)

    ncatvars <- length(ind_predictors)
    code <- "transformed data {
 matrix[N, K_X] Xc;  // centered version of X
 matrix[ps_N, K_X] ps_Xc;  // centered version of ps_X
 vector[K_X] means_X;  // column means of X
 int seq[N] = sequence(1, N);
"
    code <- paste0(code,
"
for (i in 1:K_X) {
    means_X[i] = mean(X[, i]);
    Xc[, i] = X[, i] - means_X[i];
    ps_Xc[, i] = ps_X[, i] - means_X[i];
  }
}
")
    return(code)
}

make_params_code <- function(f, data, ps, aux, adjust, overdispersed, mrp_only) {
    require(Formula)
    f <- Formula(f)
        
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)

    
### For each of the levels of the categorical outcome variable
    dv <- terms(f, lhs = TRUE, rhs = c(FALSE, FALSE))
    dv <- all.vars(dv)

    dv_levels <- levels(factor(data[,dv]))

    code <- "parameters {\n"
    for (d in dv_levels[-1]) {
        addon <- paste0(" vector[K_X] b_mu", d,
                        ";  // population-level effects\n",
                        " real Intercept_mu", d,
                        ";  // temporary intercept for centered predictors\n")
        code <- paste0(code,
                       addon)
    }

### We want parameters for each option/variable combination
    ncatvars <- length(ind_predictors)
    for (i in 1:(ncatvars+1)) {
        for (d in dv_levels[-1]) {
            addon <- paste0(" vector[N_", i, "] z_",
                            i,
                            "_",
                            d,
                            ";  // standardized group-level effects\n")
            code <- paste0(code,
                           addon)

                        addon <- paste0(" real<lower=0> sd_",
                            i,
                            "_",
                            d,
                            ";  // SD of group intercepts\n")
            code <- paste0(code,
                           addon)
        }
    }

    if (adjust) { 
        code <- paste0(code, " vector [ncat-1] adj0; \n")
    }
    if (overdispersed) { 
        code <- paste0(code,
                       "\n real <lower=0>invprec; \n")
    }
    
    
    code <- paste0(code,
                   "\n}\n\n")
    
    
    return(code)
}

make_tparams_code <- function(f, data, ps, aux, adjust, overdispersed, mrp_only) {

    require(Formula)
    f <- Formula(f)
        
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)

    
### For each of the levels of the categorical outcome variable
    dv <- terms(f, lhs = TRUE, rhs = c(FALSE, FALSE))
    dv <- all.vars(dv)

    dv_levels <- levels(factor(data[,dv]))

    code <- "transformed parameters {\n"

### We want parameters for each option/variable combination
    ncatvars <- length(ind_predictors)
    for (i in 1:(ncatvars+1)) {
        for (d in dv_levels[-1]) {
            addon <- paste0(" vector[N_", i, "] r_",
                            i,
                            "_",
                            d,
                            ";  // real (not standardized) group-level effects\n")
            code <- paste0(code,
                           addon)

        }
    }

    if (!mrp_only) { 
        code <- paste0(code,
                       " matrix[nAreas, ncat] aggmu;\n")
    }
    

    if (adjust) { 
        code <- paste0(code, " vector [ncat] adj; \n")
    }
    if (overdispersed) {
        code <- paste0(code, " real<lower=0> prec; \n")
    }
    
    for (i in 1:(ncatvars+1)) {
        for (d in dv_levels[-1]) {
            addon <- paste0(" r_",
                            i,
                            "_",
                            d,
                            " = (sd_", i, "_", d, " * z_",
                            i,
                            "_",
                            d,
                            ");\n")
            code <- paste0(code,
                           addon)
            
        }
    }


    if (adjust) { 
        code <- paste0(code,
                       "adj = append_row(0, adj0); \n ")
    }
    if (overdispersed) {
        code <- paste0(code,
                       "prec = 1.0 / invprec; \n")
    }
   
### Big chunk to get aggmu
    if (!mrp_only) {
        code <- paste0(code,
                       " for (i in 1:nAreas) {\n")
        code <- paste0(code,
                       "  aggmu[i] = pred_prob(")
        if (adjust) {
            code <- paste0(code,
                           "adj, ")
        }
        
        code <- paste0(code,
                       "areastart[i], areastop[i], ncat, Y, ")

        code <- paste0(code,
                       "ps_Xc, ")

        for(d in dv_levels[-1]) {
            addon <- paste0(" b_mu", d, ", ")
            code <- paste0(code,
                           addon)
            addon <- paste0(" Intercept_mu", d, ", ")
            code <- paste0(code,
                           addon)
            for (k in 1:(ncatvars+1)) {
                addon <- paste0(" r_", k, "_", d, ", ")
                code <- paste0(code,
                               addon)
            }
            
        }

        for (i in 1:(ncatvars+1)) {
            addon <- paste0("ps_J_", i, ", ")
            code <- paste0(code, addon)
        }

        code <- paste0(code,
                       " ps_counts)'; \n }")

    }
    ### End block for aggmu
    code <- paste0(code,
                   "\n}\n\n")

    return(code)    
}

make_model_code <- function(f, data, ps, aux, res, adjust, overdispersed, threading, mrp_only) {

    require(Formula)
    f <- Formula(f)
        
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)

    ncatvars <- length(ind_predictors)
    
### For each of the levels of the categorical outcome variable
    dv <- terms(f, lhs = TRUE, rhs = c(FALSE, FALSE))
    dv <- all.vars(dv)

    dv_levels <- levels(factor(data[,dv]))

    code <- "model { \n"

    if (!mrp_only) { 
        if (overdispersed) { 
            code <- paste0(code,
                           "if (!prior_only) {\n",
                           " for (i in 1:nAreas) {\n",
                           "  aggy[i] ~ dirichlet_multinomial(",
                           "prec * to_vector(aggmu[i]));",
                           "\n }\n}\n\n")
        } else {
            code <- paste0(code,
                           "if (!prior_only) {\n",
                           " for (i in 1:nAreas) {\n",
                           "  aggy[i] ~ multinomial(",
                           "to_vector(aggmu[i]));",
                           "\n }\n}\n\n")
        }

    }
    
    
    code <- paste0(code,
                   "if (!prior_only) {\n")

    if (threading) { 
        code <- paste0(code,
                       " target += reduce_sum(partial_log_lik, seq, grainsize, ncat, Y, Xc, ")
    } else {
        code <- paste0(code,
                       " target += partial_log_lik(sequence(1, N), 1, N, ncat, Y, Xc, ")
    }
    

    for(d in dv_levels[-1]) {
        addon <- paste0(" b_mu", d, ", ")
        code <- paste0(code,
                       addon)
        addon <- paste0(" Intercept_mu", d, ", ")
        code <- paste0(code,
                       addon)
        for (k in 1:(ncatvars+1)) {
            addon <- paste0(" r_", k, "_", d, ", ")
            code <- paste0(code,
                           addon)
        }
        
    }

    code <- paste0(code, "\n")
    
    for (i in 1:(ncatvars+1)) {
        addon <- paste0("J_", i)
        if (i != (ncatvars+1)) {
            addon <- paste0(addon, ", ")
        }
        
        code <- paste0(code, addon)
    }

### Weights?
    
    code <- paste0(code,
                   ");\n}\n")

### Prior specifications
### for each beta, normal(0, 1)
    for (d in dv_levels[-1]) {
        code <- paste0(code,
                       "for (i in 1:K_X) {\n")
        code <- paste0(code,
                       paste0(" target += normal_lpdf(b_mu", d, "[i] | 0, 1);\n"))
        code <- paste0(code,
                       "}\n\n")
    }

### Prior specifications for intercepts
    for (d in dv_levels[-1]) {
        prior_mean <- log(sum(res[,d]) / sum(res[,dv_levels[1]]))
        code <- paste0(code,
                       paste0(" target += normal_lpdf(Intercept_mu",
                              d,
                              " | ",
                              prior_mean,
                              ", 2.5);\n"))
    }

### Standard normals for the standardized effects
    for (i in 1:(ncatvars+1)) {
        for (d in dv_levels[-1]) {
            code <- paste0(code,
                           paste0("target += std_normal_lpdf(z_", i, "_", d, ");\n"))
        }
    }

    for (i in 1:(ncatvars+1)) {
        for (d in dv_levels[-1]) {
            code <- paste0(code,
                           paste0("target += normal_lpdf(sd_", i, "_", d, "|  0, 2.5) - 1 * normal_lccdf(0 | 0, 2.5);\n"))
        }
    }

### Prior specification for precision parameter
    if (overdispersed) {
        sd_in_agg_data <- sd(unlist(res[,dv_levels]))
        code <- paste0(code,
                       paste0("target += normal_lpdf(invprec | 0, ", sd_in_agg_data, ");\n"))
        
    }
    
### Finish it off
    code <- paste0(code,
                   "} \n\n")
    return(code)

    
}

make_genquant_code <- function(f, data, ps, aux, adjust) {

    require(Formula)
    f <- Formula(f)
        
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)

    ncatvars <- length(ind_predictors)

    ### For each of the levels of the categorical outcome variable
    dv <- terms(f, lhs = TRUE, rhs = c(FALSE, FALSE))
    dv <- all.vars(dv)

    dv_levels <- levels(data[,dv])
    
    code <- "generated quantities { \n"

### Initialize psw_counts
    code <- paste0(code,
                   " int psw_counts[ps_N, ncat];\n")

### Initialize category-specific counts
    for (i in 1:(ncatvars+1)) {
        code <- paste0(code,
                       paste0(" int ps_J_", i, "_counts[N_", i, ", ncat];\n"))
    }

    code <- paste0(code,
                   paste0(" int ps_area_counts[nAreas, ncat];\n"))
    code <- paste0(code,
                   "\n\n\n")

### Set them to zero
    for (i in 1:(ncatvars+1)) {
        code <- paste0(code,
                       "for (i in 1:N_", i, ") {\n")
        code <- paste0(code,
                       " for (k in 1:ncat) {\n")
        code <- paste0(code,
                       "  ps_J_", i, "_counts[i, k] = 0;\n")
        code <- paste0(code,
                       "\n }\n}\n")
    }

    code <- paste0(code,
                   "for (i in 1:nAreas) {\n")
    code <- paste0(code,
                   " for (k in 1:ncat) {\n")
    code <- paste0(code,
                   "  ps_area_counts[i, k] = 0;\n")
    code <- paste0(code,
                   "\n }\n}\n")
    
### Add everything on to psw_counts
    code <- paste0(code,
                   "for (i in 1:ps_N) {\n")
    
    code <- paste0(code,
                   " row_vector [ncat] tmp = pred_prob(")
    if (adjust) {
        code <- paste0(code,
                       "adj, ")
    }

    code <- paste0(code,
                   "i, i, ncat, Y, ")

    code <- paste0(code,
                   "ps_Xc, ")

    for (d in dv_levels[-1]) {
        code <- paste0(code, "\n")
        addon <- paste0("  b_mu", d,
                        ", Intercept_mu", d, ", ")
        code <- paste0(code, addon)
        for (i in 1:(ncatvars+1)) {
            addon <- paste0("r_", i, "_", d, ", ")
            code <- paste0(code, addon)
        }
    }

    code <- paste0(code, "\n  ")
    for (i in 1:(ncatvars+1)) {
        addon <- paste0("ps_J_", i)
        if (i != (ncatvars+1)) {
            addon <- paste0(addon, ", ")
        }
        
        code <- paste0(code, addon)
    }
    code <- paste0(code, ", ps_counts)';\n\n")

    code <- paste0(code,
                   " psw_counts[i] = multinomial_rng(to_vector(tmp), ps_counts[i]);")
    code <- paste0(code,
                   "\n}\n")

### Start on adding on things
    for (i in 1:(ncatvars+1)) {
        code <- paste0(code,
                       "for (p in 1:ps_N) {\n  for (k in 1:ncat) {\n")

        code <- paste0(code,
                       "   ps_J_", i, "_counts[ps_J_", i, "[p],k] = ps_J_", i, "_counts[ps_J_", i, "[p],k] + psw_counts[p, k];\n")
        code <- paste0(code,
                       " }\n}\n")
    }

### Do the area counts
### Start on adding on things
    code <- paste0(code,
                   "for (p in 1:ps_N) {\n")
    code <- paste0(code,
                   " for (k in 1:ncat) {\n")
    
    code <- paste0(code,
                   "   ps_area_counts[ps_area[p],k] = ps_area_counts[ps_area[p],k] + psw_counts[p, k];\n")
    code <- paste0(code,
                   " }\n}\n")
    
    code <- paste0(code,
                   "}\n\n")
    return(code)
}
