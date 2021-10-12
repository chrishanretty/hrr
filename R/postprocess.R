postprocess <- function(obj, f, data, ps, areavar, weightvar, catlu, probs = c(0.05, 0.95)) {
    fit <- obj$fit
    
    retval <- list()
    retval$fit <- fit

    
### ##################    
### (1) Area
### ##################
    
    if (inherits(fit, "stanfit")) {
        print("Summarizing...")
        area_smry <- rstan::summary(fit, probs = probs)
        area_smry <- area_smry$summary
        area_smry <- as.data.frame(area_smry)
        area_smry$param <- rownames(area_smry)
        rownames(area_smry) <- NULL
    } else {
        stop("Post-processing cmdstanr objects not yet supported")
    }
    area_smry <- area_smry[grepl("ps_area_counts", area_smry$param), ]
    area_smry$area_idx <- sub(".*\\[(.*),.*", "\\1", area_smry$param)
    area_smry$area_idx <- as.numeric(area_smry$area_idx)
    area_smry$yvalue_idx <- sub(".*\\[.*?,(.*)\\]", "\\1", area_smry$param)
    area_smry$yvalue_idx <- as.numeric(area_smry$yvalue_idx)
    
### Replace index entries with levels
    areas <- obj$areas
    area_smry$area <- areas[area_smry$area_idx]

    depvar_levels <- obj$depvar_levels
    area_smry$y.value <- depvar_levels[area_smry$yvalue_idx]

### Add on the total count (from ps) for this area
    cts <- aggregate.data.frame(ps[, weightvar], list(area_idx = ps[,areavar]), sum)
    area_smry <- merge(area_smry, cts,
                       by.x = "area_idx",
                       by.y = "area_idx",
                       all.x = TRUE,
                       all.y = FALSE)


    prob_cols <- paste0(prettyNum(probs * 100), "%")
    if (any(!is.element(prob_cols, colnames(area_smry)))) {
        stop("A formatting error happened trying to access summary values using values in probs. ")
    }
    
    for (v in c("mean", "sd", prob_cols)) {
        area_smry[, v] <- area_smry[, v] / area_smry$x
    }

    retval$area_smry <- area_smry[,c("area", "y.value",
                                     "mean", "sd",
                                     prob_cols,
                                     "n_eff", "Rhat")]



### ###################   
### (2) Characteristics
### ###################
    if (inherits(fit, "stanfit")) { 
        grp_smry <- rstan::summary(fit, probs = probs)
        grp_smry <- grp_smry$summary
        grp_smry <- as.data.frame(grp_smry)
        grp_smry$param <- rownames(grp_smry)
        rownames(grp_smry) <- NULL
        
    } else {
        stop("Post-processing cmdstanr objects not yet supported")
    }
    grp_smry <- grp_smry[grepl("ps_J_\\d+_counts", grp_smry$param), ]
    grp_smry$var_idx <- sub("ps_J_(\\d+).*", "\\1", grp_smry$param)
    grp_smry$var_idx <- as.numeric(grp_smry$var_idx)
    grp_smry$level_idx <- sub(".*counts\\[(.*),.*", "\\1", grp_smry$param)
    grp_smry$level_idx <- as.numeric(grp_smry$level_idx)
    grp_smry$yvalue_idx <- sub(".*\\[.*?,(.*)\\]", "\\1", grp_smry$param)
    grp_smry$yvalue_idx <- as.numeric(grp_smry$yvalue_idx)
    
### Replace index entries with levels
    holder <- list()
    for (i in 1:length(catlu)) {
        holder[[i]] <- data.frame(var_idx = i,
                                  var_name = catlu[[i]]$var,
                                  var_level = catlu[[i]]$levels,
                                  level_idx = 1:length(catlu[[i]]$levels))
    }
    catlu_df <- do.call("rbind", holder)

    grp_smry <- merge(grp_smry,
                       catlu_df,
                       by = c("var_idx", "level_idx"),
                       all.x = TRUE,
                       all.y = FALSE)
    
    levels <- obj$depvar_levels
    grp_smry$y.value <- depvar_levels[grp_smry$yvalue_idx]

### We have counts for variable/level combinations across the different outcomes
    ### Let's turn these into proportions by dividing each count by the sum
    cts <- aggregate.data.frame(grp_smry$mean,
                                list(var_idx = grp_smry$var_idx,
                                     level_idx = grp_smry$level_idx),
                                sum)
    
    grp_smry <- merge(grp_smry, cts,
                       by.x = c("var_idx", "level_idx"),
                       by.y = c("var_idx", "level_idx"),
                       all.x = TRUE,
                       all.y = FALSE)


    prob_cols <- paste0(prettyNum(probs * 100), "%")
    if (any(!is.element(prob_cols, colnames(grp_smry)))) {
        stop("A formatting error happened trying to access summary values using values in probs. ")
    }
    
    for (v in c("mean", "sd", prob_cols)) {
        grp_smry[, v] <- grp_smry[, v] / grp_smry$x
    }

    retval$grp_smry <- grp_smry[,c("var_name", "var_level",
                                     "y.value",
                                     "mean", "sd",
                                     prob_cols,
                                     "n_eff", "Rhat")]

### (3) Coefficients
### We need to get the coefficients
###
    if (inherits(fit, "stanfit")) { 
        coef_smry <- rstan::summary(fit, probs = probs)
        coef_smry <- coef_smry$summary
        coef_smry <- as.data.frame(coef_smry)
        coef_smry$param <- rownames(coef_smry)
        rownames(coef_smry) <- NULL
        
    } else {
        stop("Post-processing cmdstanr objects not yet supported")
    }

    ### (a) Replace b_mu[depvar_level][1], b_mu[depvar_level][2], etc.,
### with names of continuous variables
    cont_predictors <- terms(f, lhs = FALSE, rhs = c(FALSE, TRUE))
    cont_predictors <- all.vars(cont_predictors)

    betas <- coef_smry[grepl("^b_mu", param), ]

    for (i in 1:length(cont_predictors)) {
        betas$param <- sub(paste0("[", i, "]"),
                           paste0("_", cont_predictors[i]),
                           betas$param,
                           fixed = TRUE)
    }

### (b) Replace z_1_[depvar_level][catlevel] with names of categorical variables
    ### and levels that are stored in catlu
    zs <- coef_smry[grepl("^z_", param), ]
    
    for (i in 1:length(catlu)) {
### First replace the entry in square brackets with the level
        for (j in 1:length(catlu[[i]]$levels)) { 
            zs$param <- sub(paste0("z_", i, "_(.*)\\[", j, "\\]"),
                            paste0("z_", i, "_\\1\\[", catlu[[i]]$levels[j], "\\]"),
                            zs$param)
        }
        
### Replace z_i_ with z_[varname]        
        zs$param <- sub(paste0("z_", i),
                        paste0("z_", catlu[[i]]$var),
                        zs$param,
                        fixed = TRUE)
        
    }

    coefs <- merge(betas, zs, all = TRUE)
    coefs <- coefs[order(coefs$param),]
    retval$coefs <- coefs
    
    return(retval)
    
}
