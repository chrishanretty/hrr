#' Adjust an MRP model to match known aggregate outcomes
#'
#' \code{adjust_mrp} overwrites the random area intercepts contained
#' in an existing model.
#'
#' @param obj an object created by calling `hrr` with argument `mrp_only` set to TRUE
#' @param nthin number of iterations to keep
#'
#' @export
adjust_mrp <- function(obj, nthin = 25) {
### Thin everything down
### Get the predicted probabilities
### get adjustment factors
### Adjust the coefficients accordingly (intercept, area random effects)
### Recalculate group support and area support

### Thin everything down
    message("Thinning fit...")
    obj$fit <- sub_sample(obj$fit, nthin, keep_warmup = FALSE)
    
### Get the the latent scale
    mu <- get_mu(obj)

    ### oldsup <- recreate_support(obj)
### Get adjustment factors
### About half an hour for 100 iters
    message("Finding adjustment parameters...\n")
    counts <- fudge(mu, obj)

    counts <- data.frame(name = names(counts),
                         count = counts,
                         row.names = NULL)

    counts$var_idx <- sub("ps_J_([0-9]+).*", "\\1", counts$name)
### This will not work for area counts
    counts$var_idx[grepl("area_counts", counts$name)] <- length(obj$catlu) + 1
    
    counts$iter <- sub(".*_counts\\[([0-9]+),.*", "\\1", counts$name)
    counts$var_level_idx <- sub(".*_counts\\[[0-9]+,([0-9]+),.*", "\\1", counts$name)

    counts$party <- sub(".*,([0-9]+)\\]$", "\\1", counts$name)

    counts$name <- NULL
    counts$var_idx <- as.numeric(counts$var_idx)
    counts$var_level_idx <- as.numeric(counts$var_level_idx)
    counts$iter <- as.numeric(counts$iter)
    counts$party <- as.numeric(counts$party)

### Label these things
### Add area to catlu
    catlu <- obj$catlu
    catlu[[length(obj$catlu) + 1]] <- list(var = "area",
                                           levels = obj$areas)

    ### Name the variables
    counts$var <- sapply(catlu, function(x)x[["var"]])[counts$var_idx]
    counts$var_level <- ""

    ### Name the variable levels
    for (i in 1:length(catlu)) {
        counts$var_level[which(counts$var_idx == i)] <-
            catlu[[i]]$levels[counts$var_level_idx[which(counts$var_idx == i)]]
    }

### Name the parties
    counts$y.value <- colnames(obj$data$aggy)[counts$party]

### Remove the indices
    counts$var_idx <- NULL
    counts$var_level_idx <- NULL

    ### counts.bak <- counts
### Aggregate by iter, var, and var_level
    ### Do this to create counts which we can merge back in
    tots <- aggregate(counts$count,
                      by = list(var = counts$var,
                                iter = counts$iter,
                                var_level = counts$var_level),
                      sum)

    counts <- merge(counts, tots, all = TRUE)
    counts$prop <- counts$count / counts$x
    
### Aggregate by var_name, var_level, y.value
### summarizing with mean, sd, `2.5%`, `97.5%
    out <- aggregate(prop ~ y.value + var + var_level,
                     data = counts,
                     FUN = function(x) c(meanval = mean(x),
                                         sdval = sd(x),
                                         lo = quantile(x, 0.025),
                                         med = quantile(x, 0.5),
                                         hi = quantile(x, 1 - 0.025)))

    out <- cbind(out[,1:3], as.data.frame(out$prop))
    names(out) <- c(names(out)[1:3],
                    "mean", "sd", "2.5%", "50%", "97.5%")

### Make sure the columns are correctly ordered
    
    area_smry <- out[out$var == "area",]
    area_smry <- area_smry[,c("var_level", "y.value", "mean", "sd",
                              "2.5%", "50%", "97.5%")]
    names(area_smry)[1] <- "area"
    
    grp_smry <- out[out$var != "area",]
    grp_smry <- grp_smry[,c("var", "var_level", "y.value", "mean", "sd",
                              "2.5%", "50%", "97.5%")]
    names(grp_smry)[1] <- "var_name"
    
    return(list(area_smry = area_smry,
                grp_smry = grp_smry))
}


sub_sample <- function(stanfit, n, keep_warmup = FALSE) {
  sim <- stanfit@sim
  samp <- sim$samples
  W <- sim$warmup
  I <- sim$iter
  sel <- c(if (keep_warmup) 1:W, sample((W + 1):I, size = n))
  subsamp <- lapply(samp, function(chain_samp) {
    lapply(chain_samp, function(x) x[sel])
  })
  stanfit@sim$samples <- subsamp
### Adjust the n_save argument
  nchains <- length(stanfit@sim$n_save)
  n_per_chain <- n / nchains
  
  ## if (round(n_per_chain) != n_per_chain) {
  ##     stop("Number of iterations to keep must be perfectly divisible by the number of chains")
  ## }
  stanfit@sim$n_save <- rep(n,
                            nchains)
  stanfit@sim$iter <- n
  stanfit@sim$warmup <- 0
  stanfit@sim$warmup2 <- rep(0, nchains)
  
### Overwrite the permutation
  stanfit@sim$permutation <- lapply(stanfit@sim$permutation,
                                    function(x) {
                                        1:n
                                    })
  
  stanfit
}

get_inits <- function(obj) {
### Get area by party
    area_smry <- obj$area_smry
    area_smry <- area_smry[, c("area", "y.value", "mean")]
    area_smry <- reshape(area_smry,
                         direction = "wide",
                         idvar = "area",
                         timevar = "y.value")
    names(area_smry)[-1] <- sub("mean.", "",
                                names(area_smry)[-1],
                                fixed = TRUE)
    area_smry$area <- NULL
    
### Compare to agg results
    aggy <- obj$data$aggy
    aggy <- aggy / rowSums(aggy)

### Convert both to log odds ratios
    area_smry <- log(area_smry[,-1] / area_smry[,1])
    aggy <- log(aggy[,-1] / aggy[,1])

### Get differences on these scales
    delta <- aggy - area_smry

### Sweep out the global adjustment
    if (is.null(dim(delta))) {
### We have a binary outcome
        global_adj <- array(mean(delta), dim = 1)
        area_adj <- matrix(delta - mean(delta), ncol = 1)
    } else {
        global_adj <- colMeans(delta)
        area_adj <- delta - colMeans(delta)
    }
    return(list(global_adj = global_adj,
                area_adj = area_adj))
}

fudge <- function(mu, obj, debug = FALSE) {
    
    param_names <- names(obj$fit)
    nIter <- nrow(mu)

    nParties <- ifelse(is.na(dim(mu)[3]),
                       2,
                       1 + dim(mu)[3])

    if (nParties == 2) {
        newmu <- array(0, dim = c(nrow(mu), ncol(mu), 1))
        newmu[,,1] <- mu
        mu <- newmu
        rm(newmu)
    }
    
    stan_data <- list(nIter = nIter,
                  nAreas = length(unique(obj$data$ps_area)),
                  nParties = nParties,
                  ps_N = obj$data$ps_N,
                  ps_w8 = obj$data$ps_counts,
                  area_idx = obj$data$ps_area,
                  area_start = aggregate(1:obj$data$ps_N,
                                  list(a = obj$data$ps_area),
                                  min)[,"x"],
                  area_stop = aggregate(1:obj$data$ps_N,
                                  list(a = obj$data$ps_area),
                                  max)[,"x"],
                  aggy = obj$data$aggy,
                  debug = as.numeric(debug))
    
    stan_data$mu <- mu
    
    data_code <- "
data {
     int<lower=1>ps_N;
     int<lower=0>ps_w8[ps_N];
     int<lower=1>nAreas;
     int<lower=2>nParties;
     int<lower=1>nIter;
     int<lower=0, upper = 1> debug;
     int<lower=1, upper = nAreas> area_idx[ps_N];
     int<lower=1, upper = ps_N> area_start[nAreas];
     int<lower=1, upper = ps_N> area_stop[nAreas];
     int<lower=0> aggy[nAreas, nParties];
     real mu[nIter, ps_N, nParties - 1];

"

### Add on the number of categories for each entry
    for (i in 1:length(obj$catlu)) {
        addon_code <- paste0(" int<lower=1> N_", i, "; \n")
        data_code <- paste(data_code, addon_code)
        stan_data[[paste0("N_", i)]] <- obj$data[[paste0("N_", i)]]
    }

### Add on the indices from the post-stratification frame
    for (i in 1:length(obj$catlu)) {
        addon_code <- paste0(" int<lower=1, upper = N_", i, "> ps_J_", i, "[ps_N]; \n")
        data_code <- paste(data_code, addon_code)
        stan_data[[paste0("ps_J_", i)]] <- obj$data[[paste0("ps_J_", i)]]
    }
### Add on the closing bracket
    data_code <- paste0(data_code, "\n}\n")
    
stan_code <- "
parameters {
     real global_adj[nParties - 1];
     real area_adj[nAreas, nParties - 1];
     real adj[nIter, nAreas, nParties - 1];
}
transformed parameters {
    real<lower=0> theta[nIter, nAreas, nParties]; // counts and proportions 
    real newmu[nIter, ps_N, nParties];

    // Create the adjusted mu (newmu)
    for (i in 1:nIter) {
       for (j in 1:ps_N) {

          newmu[i, j, 1] = 0; // zero for the reference parties
          for (k in 2:nParties) { 
             newmu[i, j, k] = mu[i, j, k-1] +
                global_adj[k-1] +
                area_adj[area_idx[j], k-1] +
                adj[i, area_idx[j], k-1];
          }

          newmu[i, j, 1:nParties] = to_array_1d(softmax(to_vector(newmu[i, j, 1:nParties])));
       }
    }

    // Aggregate for different slices
    for (i in 1:nIter) { 
       for (j in 1:nAreas) {
          int start = area_start[j];
          int stop = area_stop[j];
          // initialize
	  for (k in 1:nParties) {
   	     theta[i, j, k] = 0.0;
   	     for (m in start:stop) {
	        theta[i, j, k] += ps_w8[m] .* newmu[i, m, k];
             }
          }
       }
    }

    // Transform theta to simplex
    for (i in 1:nIter){
     for (j in 1:nAreas) {
       real thetasum = sum(theta[i,j,1:nParties]);
       for (k in 1:nParties) { 
         theta[i,j,k] = theta[i,j,k] / thetasum;
       }
     }
    }
       
       
}
model {
   to_vector(global_adj) ~ normal(0, 2.5);
   for (i in 1:nIter) {
      for (j in 1:nAreas) {
         to_vector(adj[i, j, ]) ~ std_normal();
         aggy[j, ] ~ multinomial(to_vector(theta[i, j, ]));
      }
   }
   for (j in 1:nAreas) {
      to_vector(area_adj[j,]) ~ std_normal();
   }
}
"
### Add on to the stan code
    ## We'll cannibalize the code from the object
    genquant_code <- "
generated quantities {
int<lower=0>pr[nIter, ps_N, nParties];
"
    ### Set-up count objects
    for (i in 1:length(obj$catlu)) {
        addon_code <- paste0(" int<lower=0> ps_J_", i, "_counts[nIter, N_", i, ", nParties]; \n")
        genquant_code <- paste0(genquant_code, addon_code)
    }
    addon_code <- " int<lower=0> ps_area_counts[nIter, nAreas, nParties]; \n"
    genquant_code <- paste0(genquant_code, addon_code)

### Start initializing all these
    for (i in 1:length(obj$catlu)) {
        addon_code <- paste0("for (i in 1:nIter) {
   for (j in 1:N_", i, ") {
      for (k in 1:nParties) {
         ps_J_", i, "_counts[i, j, k] = 0;
      }
   }
}\n\n")
        genquant_code <- paste0(genquant_code, addon_code)
    }
    
### Initialize area
    addon_code <- "
for (i in 1:nIter) {
   for (j in 1:nAreas) {
      for (k in 1:nParties) {
         ps_area_counts[i, j, k] = 0;
      }
   }
}

"

    genquant_code <- paste0(genquant_code, addon_code)

### Create predicted counts
    addon_code <- "
    for (i in 1:nIter) {
        for (p in 1:ps_N) {
            pr[i,p, 1:nParties] = multinomial_rng(to_vector(newmu[i, p, ]), ps_w8[p]);
        }
    }
    "
    genquant_code <- paste0(genquant_code, addon_code)
    
### Start adding predicted counts on
    for (i in 1:length(obj$catlu)) {
        addon_code <- paste0("
for (i in 1:nIter) {
   for (p in 1:ps_N) {
      for (k in 1:nParties) { 
         ps_J_", i, "_counts[i, ps_J_", i, "[p],k] = ps_J_", i, "_counts[i, ps_J_", i, "[p],k] + pr[i, p, k];
      }
   }
}

")
        genquant_code <- paste0(genquant_code, addon_code)
    }

### Do area counts
      addon_code <- paste0("
for (i in 1:nIter) {
   for (p in 1:ps_N) {
      for (k in 1:nParties) { 
         ps_area_counts[i, area_idx[p],k] = ps_area_counts[i, area_idx[p],k] + pr[i, p, k];
      }
   }
}

")

    genquant_code <- paste0(genquant_code, addon_code)

### Closing bracket!
    genquant_code <- paste0(genquant_code, "\n\n}\n\n")
    stan_code <- paste0(data_code,
                        stan_code,
                        genquant_code)
    


    sm <- rstan::stan_model(model_code = stan_code)

### Build the parameter list
    ## Oh no, optimizing doesn't support pars!
    ## pars <- c(paste0("ps_J_", 1:length(obj$catlu), "_counts"),
    ##           "ps_area_counts")
    so <- rstan::optimizing(sm,
                            data = stan_data,
                            verbose = FALSE,
                            init = get_inits(obj))

    if (debug) {
        theta <- so$par[grep("theta", names(so$par))]
        theta <- array(theta,
                       dim = c(stan_data$nIter,
                               stan_data$nAreas,
                               stan_data$nParties))
### Compare this to actual
        aggy <- stan_data$aggy / rowSums(stan_data$aggy)
        rmse_holder <- list()
        for (i in 1:stan_data$nIter) {
            rmse_holder[[i]] <- sqrt(mean((theta[i, , ] - as.matrix(aggy)) ^ 2))
        }
        pr <- so$par[grep("^pr", names(so$par))]        
    }
    

    if (so$return_code > 0) {
        stop("L-BFGS optimization failed: try thinning more")
        message("LBFGS optimization failed, trying Newton")
        so <- rstan::optimizing(sm,
                                algorithm = "Newton",
                                data = stan_data,
                                init = 0)
    }

### Get the different supports
    counts <- grep("_counts\\[", names(so$par))
    counts <- so$par[counts]
    return(counts)

}



get_mu <- function(obj) {
    depvar_levels <- unique(obj$area_smry$y.value)
    if (length(depvar_levels) > 2) {
        f <- get_mu.categorical
    } else {
        f <- get_mu.binary
    }
    f(obj)
}

get_mu.binary <- function(obj) {

    param_names <- names(obj$fit)
    alpha <- rstan::extract(obj$fit,
                            pars = grep("Intercept_mu", param_names, value = TRUE))[[1]]
    
    ## Get continuous predictors
    beta <- rstan::extract(obj$fit,
                           pars = grep("^b_", param_names, value = TRUE))
    beta <- do.call("cbind", beta)
    
    ## Get categorical variables
    ## First get their number
    rvars <- grep("^r_", param_names, value = TRUE)
    the_rvar <- sub("r_([0-9]+)_.*", "\\1", rvars)
    ncatvars <- max(as.numeric(the_rvar))
    cat_coefs <- list()
    for (i in 1:ncatvars) {
        coefs <- rstan::extract(obj$fit,
                                pars = grep(paste0("r_", i),
                                            param_names,
                                            value = TRUE))
        cat_coefs[[i]] <- do.call("cbind", coefs)
    }
    
### Construct this using the stuff in the data object
    nIter <- length(alpha)
    mu <- matrix(0,
                 nrow = nIter,
                 ncol = obj$data$ps_N)
### Construct by iteration
    for (i in 1:nIter) {
        mu[i,] <- alpha[i] +
            tcrossprod(beta[i, ], obj$data$ps_X)
        for (j in 1:ncatvars) {
            mu[i,] <- mu[i,] +
                cat_coefs[[j]][i, obj$data[[paste0("ps_J_", j)]]]
        }
    }

    return(mu)
}


get_mu.categorical <- function(obj) {
    message("Getting predictions on latent scale (categorical)")
    param_names <- names(obj$fit)

    ### Start with the intercept
    alpha <- rstan::extract(obj$fit,
                            pars = grep("Intercept_mu", param_names, value = TRUE))
    alpha <- do.call("cbind", alpha)

    nIter <- nrow(alpha)

### Start learning about the model we've got
### Here, nCats is number of modelled categories
### nCats = nParties - 1
    nCats <- ncol(alpha)
    zvars <- grep("^z_", param_names, value = TRUE)
    the_zvar <- as.numeric(sub("z_([0-9]+)_.*", "\\1", zvars))
    ncatvars <- max(the_zvar)

### Get the dep var levels
### This has length nCats - 1 (b/c we don't know the reference cat.)
    depvar_levels <- sub("z_[0-9]+_(.*?)\\[.*", "\\1", zvars)
    depvar_levels <- unique(depvar_levels)

    ## Get continuous predictors, storing them as a named list, where
    ## names are depvar_levels
    beta <- vector(mode = "list", length = length(depvar_levels))
    names(beta) <- depvar_levels
    for (i in depvar_levels) {
        matches <- grep(paste0("^b_mu", i),
                         param_names,
                         value = TRUE)
        beta[[i]] <- do.call("cbind",
                             rstan::extract(obj$fit, pars = matches))
    }
 
### Get categorical variables
### These will be stored in a list-of-lists
### first-level: depvar_levels
    ## second-level: number of variable
    eta <- vector(mode = "list", length = length(depvar_levels))
    names(eta) <- depvar_levels

    for (i in depvar_levels) {
        matches <- grep(paste0("^r_[0-9]+_", i),
                         param_names,
                        value = TRUE)
        n_cat_vars <- sub("r_([0-9]+)_.*", "\\1", matches)
        n_cat_vars <- max(as.numeric(n_cat_vars))

        tmp <- list()
        for (j in 1:n_cat_vars) {
            matches <- grep(paste0("^r_", j, "_", i),
                            param_names,
                            value = TRUE)
            tmp[[j]] <- do.call("cbind",
                                rstan::extract(obj$fit, pars = matches))
        }
        eta[[i]] <- tmp
    }


### Start assembling things
    mu <- array(0, dim = c(nIter, obj$data$ps_N, nCats))
    
    for (i in depvar_levels) {
        cat_posn <- which(depvar_levels == i)
### Promote alpha to matrix
        alpha_mat <- matrix(alpha[,cat_posn],
                            nrow = nIter,
                            ncol = obj$data$ps_N,
                            byrow = FALSE)
        mu[,,cat_posn] <- alpha_mat +
            tcrossprod(beta[[i]], obj$data$ps_X)
        for (j in 1:n_cat_vars) {
            cat_level <- obj$data[[paste0("ps_J_", j)]]
            mu[,,cat_posn] <- mu[,,cat_posn] +
                eta[[i]][[j]][, cat_level]
        }
    }

### Add on the zero column for the reference category
    ### 
    return(mu)
}


recreate_support <- function(obj) {
    depvar_levels <- unique(obj$area_smry$y.value)
    if (length(depvar_levels) > 2) {
        f <- recreate_support.categorical
    } else {
        f <- recreate_support.binary
    }
    f(obj)
}

recreate_support.binary <- function(obj) {

    mu <- get_mu(obj)
    ## convert to counts
    ## latent to predicted probability scale
    pr <- plogis(mu)
    counts <- pr * matrix(obj$data$ps_counts,
                          nrow = nrow(mu),
                          ncol = length(obj$data$ps_counts),
                          byrow = TRUE)
    
### Summarize these counts by area and group
    nAreas <- length(unique(obj$area_smry$area))
    nGroups <- length(unique(obj$grp_smry$var_name))
    
    area_props <- apply(counts, 1, function(x) {
        counts <- tapply(x,
                         obj$data$ps_area,
                         sum)
        totals <- tapply(obj$data$ps_counts,
                         obj$data$ps_area,
                         sum)
        counts / totals
    })

    group_props <- vector(mode = "list",
                          length = nGroups)
    for (i in 1:nGroups) {
        ps_var <- paste0("ps_J_", i)
        group_props[[i]] <- apply(counts, 1, function(x) {
            counts <- tapply(x,
                             obj$data[[ps_var]],
                             sum)
            totals <- tapply(obj$data$ps_counts,
                             obj$data[[ps_var]],
                             sum)
            counts / totals
        })
    }

    ## Tidy these properly
    ## Tidy area proportions
    ## At the moment, area_props is a matrix with dimensions nAreas by nSims
    ## area_smry should be a data frame with columns
    ## area, y.value, mean, sd, 2.5%, 50%, 97.5%

    probs <- grep("%", names(obj$area_smry), value = TRUE)
    probs_labels <- probs
    probs <- as.numeric(sub("%", "", probs)) / 100
    area_props <- apply(area_props, 1, function(x) {
        tmp <- data.frame(mean = mean(x),
                   sd = sd(x))
        for (p in 1:length(probs)) {
            tmp[[probs_labels[[p]]]] <- quantile(x, probs[p], na.rm = TRUE)
        }
        tmp
    })

    area_props <- do.call("rbind", area_props)
    area_props$area <- obj$areas
    ## Assume that the y-level is the second such in obj$data$aggy
    area_props$y.value <- colnames(obj$data$aggy)[2]

### Create a duplicate which has all this but reversed
    alter <- area_props
    for (i in c("mean", probs_labels)) {
        alter[, i] <- 1 - alter[, i]
    }
    alter$y.value <- colnames(obj$data$aggy)[1]

### Also need to reverse the order of the probs_labels
    names(alter)[match(probs_labels, names(alter))] <-
        rev(probs_labels)
    
    area_props <- merge(area_props,
                        alter,
                        all = TRUE)

    ## Rearrange variables
    area_props <- area_props[, c("area",
                                 "y.value",
                                 "mean", "sd",
                                 probs_labels)]

    for (i in 1:nGroups) {
        group_props[[i]] <- apply(group_props[[i]], 1, function(x) {
            tmp <- data.frame(mean = mean(x),
                              sd = sd(x))
            for (p in 1:length(probs)) {
                tmp[[probs_labels[[p]]]] <- quantile(x, probs[p], na.rm = TRUE)
            }
            tmp
        })
        group_props[[i]] <- do.call("rbind",
                                    group_props[[i]])
        group_props[[i]]$var_name <- obj$catlu[[i]]$var
        group_props[[i]]$var_level <- obj$catlu[[i]]$levels
        group_props[[i]]$y.value <- colnames(obj$data$aggy)[2]

        alter <- group_props[[i]]
        for (j in c("mean", probs_labels)) {
            alter[, j] <- 1 - alter[, j]
        }
        alter$y.value <- colnames(obj$data$aggy)[1]
        
        group_props[[i]] <- merge(group_props[[i]],
                                  alter,
                                  all = TRUE)

        ## Rearrange variables
        group_props[[i]] <- group_props[[i]][, c("var_name",
                                                 "var_level",
                                                 "y.value",
                                                 "mean", "sd",
                                                 probs_labels)]
    }

    group_props <- do.call("rbind", group_props)


    return(list(area_smry = area_props,
                group_smry = group_props))
}

recreate_support.categorical <- function(obj) {

    mu <- get_mu(obj)
    nIter <- dim(mu)[1]
    ## Convert to matrix to speed things up
    ## mumat now has dimensions (nthin * npsw) by nParties
    ## mumat works through iterations first, then rows
    ## so we have iter 1, row 1, iter 2, row 1, iter 3, row
    ## before we have anything from row 2
    mumat <- matrix(mu, ncol = dim(mu)[3])

### Add on reference category
    mumat <- cbind(0, mumat)
    emu <- exp(mumat)
    pr <- emu / rowSums(emu)
    
    ## convert to counts
    count_mat <- matrix(rep(obj$data$ps_counts,
                            each = nIter),
                        nrow = nrow(pr),
                        ncol = ncol(pr),
                        byrow = FALSE)
    counts <- pr * count_mat
    
    rm(pr); rm(mu); rm(mumat); rm(count_mat); rm(emu)

    ### Summarize these counts by area and group
    nAreas <- length(unique(obj$area_smry$area))
    nGroups <- length(unique(obj$grp_smry$var_name))
    
    area_counts <- aggregate(counts,
                             list(area = rep(obj$data$ps_area, each = nIter),
                                  iter = rep(1:nIter, times = length(obj$data$ps_area))),
                             sum)
    area_totals <- aggregate(obj$data$ps_counts,
                             list(area = obj$data$ps_area),
                             sum)
    
    for (i in 3:ncol(area_counts)) {
        area_counts[,i] <- area_counts[,i] / area_totals$x[area_counts$area]
    }

### Now tidy this up for formatting
### Get the names from the aggy variable
    colnames(area_counts)[3:ncol(area_counts)] <- colnames(obj$data$aggy)

###
    probs <- grep("%", names(obj$area_smry), value = TRUE)
    probs_labels <- probs
    probs <- as.numeric(sub("%", "", probs)) / 100

    means <- aggregate(area_counts[,3:ncol(area_counts)],
                       list(area = area_counts$area),
                       mean)
    means$quantity <- "mean"
    
    sds <- aggregate(area_counts[,3:ncol(area_counts)],
                       list(area = area_counts$area),
                     sd)
    sds$quantity <- "sd"
    
    quantiles <- list()
    for (p in 1:length(probs)) {
        quantiles[[p]] <- aggregate(area_counts[,3:ncol(area_counts)],
                                    list(area = area_counts$area),
                                    quantile,
                                    probs = probs[p],
                                    na.rm = TRUE)
        quantiles[[p]]$quantity <- probs_labels[p]
    }

    rm(area_counts)
    
    area_smry <- do.call("rbind", c(list(means, sds), quantiles))
    area_smry$area <- unique(obj$area_smry$area)[area_smry$area]

    ## Reshape this
    parties <- colnames(obj$data$aggy)
    holder <- list()
    for (p in parties) {
        tmp <- area_smry[,c("area", p, "quantity")]
        tmp <- reshape(tmp,
                       direction = "wide",
                       timevar = "quantity",
                       idvar = "area",
                       sep = "|")
        names(tmp) <- sub(".*?\\|", "", names(tmp))
        tmp$y.value <- p
        holder[[p]] <- tmp
    }
    area_smry <- do.call("rbind", holder)
    rownames(area_smry) <- NULL

    area_smry <- area_smry[,c("area", "y.value",
                              "mean", "sd",
                              probs_labels)]

### Now group summaries
    group_holder <- list()
    for (g in 1:nGroups) {
        grp_var <- paste0("ps_J_", g)
        group_counts <- aggregate(counts,
                                  list(group = rep(obj$data[[grp_var]], each = nIter),
                                       iter = rep(1:nIter, times = length(obj$data[[grp_var]]))),
                                  sum)
        
        group_totals <- aggregate(obj$data$ps_counts,
                                  list(group = obj$data[[grp_var]]),
                                  sum)
        
        for (i in 3:ncol(group_counts)) {
            group_counts[,i] <- group_counts[,i] / group_totals$x[group_counts$group]
        }
        
### Now tidy this up for formatting
### Get the names from the aggy variable
        colnames(group_counts)[3:ncol(group_counts)] <- colnames(obj$data$aggy)
        means <- aggregate(group_counts[,3:ncol(group_counts)],
                           list(group = group_counts$group),
                           mean)
        means$quantity <- "mean"
        
        sds <- aggregate(group_counts[,3:ncol(group_counts)],
                         list(group = group_counts$group),
                         sd)
        sds$quantity <- "sd"
        
        quantiles <- list()
        for (p in 1:length(probs)) {
            quantiles[[p]] <- aggregate(group_counts[,3:ncol(group_counts)],
                                        list(group = group_counts$group),
                                        quantile,
                                        probs = probs[p],
                                        na.rm = TRUE)
            quantiles[[p]]$quantity <- probs_labels[p]
        }
        
        rm(group_counts)
        
        group_smry <- do.call("rbind", c(list(means, sds), quantiles))
        group_smry$var_level <- obj$catlu[[g]]$levels[group_smry$group]

        ## Reshape this
        parties <- colnames(obj$data$aggy)
        holder <- list()
        for (p in parties) {
            tmp <- group_smry[,c("var_level", p, "quantity")]
            tmp <- reshape(tmp,
                           direction = "wide",
                           timevar = "quantity",
                           idvar = "var_level",
                           sep = "|")
            names(tmp) <- sub(".*?\\|", "", names(tmp))
            tmp$y.value <- p
            holder[[p]] <- tmp
        }
        group_smry <- do.call("rbind", holder)
        rownames(group_smry) <- NULL
        group_smry$var_name <- obj$catlu[[g]]$var
        group_smry <- group_smry[,c("var_name",
                                    "var_level", "y.value",
                                    "mean", "sd",
                                    probs_labels)]
        group_holder[[g]] <- group_smry
    }

    group_smry <- do.call("rbind", group_holder)
    return(list(area_smry = area_smry,
                group_smry = group_smry))
    
}
