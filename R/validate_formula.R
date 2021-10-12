### Should return Formula object

validate_formula <- function(f, areavar, weightvar) {
    require(Formula)
    f <- as.Formula(f)
    
### (1) Do it have an outcome variable?
    lhs <- attr(f, "lhs")
    if (is.null(lhs)) {
        stop("Formula must have an outcome variable")
    }
    
    ## (2) Does it have two-parts on the RHS?
    rhs <- attr(f, "rhs")
    if (is.null(rhs)) {
        stop("Formula must have explanatory variables")
    }
    if (length(rhs) != 2) {
        stop("Formula must have a two-part right-hand side")
    }
    
    ## (3) Does it have any (g)lmer style notation?
    rhs_terms <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    rhs_labels <- attr(rhs_terms, "term.labels")
    if (any(grepl("|", rhs_labels, fixed = TRUE))) {
        stop("(G)lmer style notation `(1|grp)` is not needed for random intercepts in the first part of the formula")
    }

    ## (4) Does it contain the constituency identifier or the weights variable
    rhs_terms <- terms(f, lhs = FALSE)
    rhs_labels <- attr(rhs_terms, "term.labels")
    if (any(rhs_labels == areavar)) {
        stop("Model formula should not include small area identifier")
    }
    if (any(rhs_labels == weightvar)) {
        stop("Model formula should not include weighting variable")
    }
    
    return(f)
}
