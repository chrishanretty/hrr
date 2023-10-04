#' Adjust an MRP model to match known aggregate outcomes
#'
#' \code{adjust_mrp} overwrites the random area intercepts contained
#' in an existing model.
#'
#' @param obj an object created by calling `hrr` with argument `mrp_only` set to TRUE
#'
#' @export
adjust_mrp <- function(obj) {

    n_areas <- nrow(obj$data$aggy)
    n_chains <- length(obj$fit@stan_args)
    n_iter <- (obj$fit@sim$iter - obj$fit@sim$warmup) *
        n_chains
    
    n_parties <- ncol(obj$data$aggy)
    
### #################################
### ASSEMBLE THE VARIABLES FROM THE FIT
### #################################
    param_names <- names(obj$fit)

    alpha_names <- grep("Intercept_mu",
                        param_names,
                        fixed = TRUE,
                        value = TRUE)

    alpha <- rstan::extract(obj$fit,
                            pars = alpha_names,
                            permute = FALSE)

    alpha <- hrr:::collapse_chains(alpha)    

    nCats <- ncol(alpha)
    zvars <- grep("^z_", param_names, value = TRUE)
    the_zvar <- as.numeric(sub("z_([0-9]+)_.*", "\\1", zvars))
    ncatvars <- max(the_zvar)

    depvar_levels <- sub("z_[0-9]+_(.*?)\\[.*", "\\1", zvars)
    depvar_levels <- unique(depvar_levels)

    ## Get continuous predictors, storing them as a named list, where
    ## names are depvar_levels
    beta <- lapply(depvar_levels, function(i) {
        matches <- grep(paste0("^b_mu", i),
                        param_names,
                        value = TRUE)
        hrr:::collapse_chains(
                  rstan::extract(obj$fit,
                                 pars = matches,
                                 permute = FALSE))
        
    })
    
    names(beta) <- depvar_levels

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

        tmp <- vector(mode = "list", length = n_cat_vars)
        for (j in 1:n_cat_vars) {
            matches <- grep(paste0("^r_", j, "_", i),
                            param_names,
                            value = TRUE)
            tmp[[j]] <- hrr:::collapse_chains(
                rstan::extract(obj$fit,
                               pars = matches,
                               permute = FALSE))
        }
        eta[[i]] <- tmp
    }

### Handle outcomes
    outcomes <- obj$data$aggy
    outcomes <- replace(outcomes,
                        outcomes == 0,
                        .Machine$double.eps^0.25)
    
### #################################
### Create the linear predictor by iter
### #################################

### Create holders for the area and group summaries
    area_holder <- vector(mode = "list", length = n_iter)
    grp_holder <- vector(mode = "list", length = n_iter)
    adj_holder <- array(0,
                        dim = c(n_iter, n_areas,
                                length(depvar_levels)))
    for (idx in seq_len(n_iter)) {
        ##     ### alpha is an array with dimensions n_iter by nChoices - 1
        local_alpha <- alpha[idx, ]
### beta is a list of size nChoices with entries arrays n_iter by dim(X)
        local_beta <- lapply(beta, function(x)x[idx,])
### eta is a list of lists
        local_eta <- lapply(eta, function(x) lapply(x, function(z)z[idx,]))
        
        mu <- matrix(0,
                     nrow = obj$data$ps_N,
                     ncol = length(depvar_levels))
        
        for (i in depvar_levels) {
            cat_posn <- which(depvar_levels == i)
            mu[,cat_posn] <- local_alpha[cat_posn] +
                tcrossprod(local_beta[[i]], obj$data$ps_X)
            for (j in 1:n_cat_vars) {
                cat_level <- obj$data[[paste0("ps_J_", j)]]
                mu[,cat_posn] <- mu[,cat_posn] +
                    local_eta[[i]][[j]][cat_level]
            }
        }
        
        ### find the per-iter adjustment
        adj <- matrix(0, nrow = n_areas, ncol = length(depvar_levels))
        
        for (j in 1:n_areas) {
            outcome <- outcomes[j, ]
            this_area <- seq.int(from = obj$data$areastart[j],
                                 to = obj$data$areastop[j],
                                 by = 1)
            this_area_mu <- mu[this_area, ]
            res <- logit_swing(outcome = unlist(outcome),
                               linpreds = this_area_mu,
                               weights = obj$data$ps_counts[this_area])
            adj[j, ] <- res
            ## Store also in the global holder
            adj_holder[idx, j, ] <- res
        }
        
        mu <- mu + adj[obj$data$ps_area, ]
        mu <- cbind(0, mu)
        ## calculate the group support        
        emu <- exp(mu)
        pr <- emu / rowSums(emu)
        
        ## convert to counts
        count_mat <- matrix(obj$data$ps_counts,
                            nrow = nrow(pr),
                            ncol = ncol(pr),
                            byrow = FALSE)
        counts <- pr * count_mat
        
        rm(pr); rm(mu); rm(count_mat); rm(emu)

            ### Summarize these counts by area and group
        nAreas <- length(unique(obj$area_smry$area))
        nGroups <- length(unique(obj$grp_smry$var_name))
        
        area_counts <- aggregate(counts,
                                 list(area = obj$data$ps_area),
                                 sum)
        area_totals <- aggregate(obj$data$ps_counts,
                                 list(area = obj$data$ps_area),
                                 sum)
        
        for (i in 2:ncol(area_counts)) {
            area_counts[,i] <- area_counts[,i] / area_totals$x[area_counts$area]
        }
        
### Now tidy this up for formatting
        colnames(area_counts)[2:ncol(area_counts)] <- colnames(obj$data$aggy)

        area_holder[[idx]] <- area_counts
        
        all_group_counts <- list()
        for (g in 1:nGroups) {
            grp_var <- paste0("ps_J_", g)
            group_counts <- aggregate(counts,
                                      list(group = obj$data[[grp_var]]),
                                      sum)
            
            group_totals <- aggregate(obj$data$ps_counts,
                                      list(group = obj$data[[grp_var]]),
                                      sum)
            
            for (i in 2:ncol(group_counts)) {
                group_counts[,i] <- group_counts[,i] / group_totals$x[group_counts$group]
            }
            all_group_counts[[g]] <- group_counts
        }
        
        grp_holder[[idx]] <- all_group_counts


    }

### add index to are and group holders
    area_holder <- lapply(seq_len(length(area_holder)), function(i) {
        tmp <- area_holder[[i]]
        tmp$iter <- i
        tmp$area <- obj$areas[1:nrow(tmp)]
        return(tmp)
    })
    area_holder <- do.call("rbind", area_holder)
    
    grp_holder <- lapply(seq_len(length(grp_holder)), function(i) {
        ### tmp is a list of lists
        tmp <- grp_holder[[i]]
        for (j in seq_len(length(tmp))) {
            names(tmp[[j]]) <- c("group",
                                 colnames(obj$data$aggy))
            tmp[[j]]$var_name <- obj$catlu[[j]]$var
            tmp[[j]]$var_level <- obj$catlu[[j]]$levels[tmp[[j]]$group]
            tmp[[j]]$iter <- i
        }
        return(do.call("rbind", tmp))
    })

    grp_holder <- do.call("rbind", grp_holder)
    
## ### Fit the adjustments back in the object?
    area_re_counter <- length(obj$catlu) + 1
    for (i in seq_len(n_chains)) {
        for (j in seq_len(nAreas)) {
            for (k in 2:n_parties) {
                local_adj <- adj_holder[,j,k-1]
### Which iters come from chain i?
                chain <- rep(1:n_chains,
                             each = length(local_adj) / n_chains)
                local_adj <- local_adj[which(chain == i)]
### Which slots in the model object must we overwrite?
                tgt_var <- paste0("r_",
                              area_re_counter,
                              "_",
                              j,
                              "[",
                              k,
                              "]")
                obj$fit@sim$samples[[i]][[tgt_var]] <- obj$fit@sim$samples[[i]][[tgt_var]] + local_adj
            }
        }
    }
    
    retval <- list("Area" = area_holder,
                   "Group" = grp_holder,
                   "adjustments" = adj_holder,
                   "model" = obj)

    return(retval)
}


logit_swing <- function(outcome, linpreds, weights) {
    stopifnot(is.vector(outcome))
    if (length(outcome) == 1) {
        stop("Outcome has to be a vector")
    }
    if (sum(outcome) != 1) {
        pr <- outcome / sum(outcome)
    }
    
    targets <- outcome[-1] / outcome[1]
    if (ncol(linpreds) != (length(outcome) - 1)) {
        stop("linear predictor should have one fewer column than outcomes")
    }
    
    objfunc <- function(delta, tgt, mu, w8) {
        pr <- exp(mu + delta)
        actual <- sum(pr * w8) / sum(w8)
        abserr <- abs(actual - tgt)
        return(abserr)
    }
    results <- rep(NA, length(targets))
    for (t in targets) {
        idx <- which(targets == t)
        oobj <- optimize(objfunc,
                                 interval = c(-3, 3),
                                 tgt = t,
                                 mu = linpreds[,idx],
                         w8 = weights)
        results[idx] <- oobj$minimum
        ### what to do when convergence not achieved?
    }
    
    return(results)

}
