compare_dfs <- function(formula,
                        df1,
                        df2) {

    ## Does formula have class formula?
    if(class(formula) != "formula") {
        formula <- as.formula(formula)
    }

    ## What about the data frames?
    ## Let's change from tibbles to base data frames
    df1 <- as.data.frame(df1)
    df2 <- as.data.frame(df2)
    
    ## Check all the variables are present
    formula_vars <- all.vars(formula)

    variables_in_df(formula_vars, df1, "df1")
    variables_in_df(formula_vars, df2, "df2")

    ## Check their class is the same
    ## Taking the first class
    ## (This deals with ordered factors, which have class "ordered", "factor"
    var_classes1 <- sapply(formula_vars, function(v) class(df1[,v])[1])
    var_classes2 <- sapply(formula_vars, function(v) class(df2[,v])[1])

    ## Check the variables have the same classes
    if (any(var_classes1 != var_classes2)) {
        clashing_classes <- which(var_classes1 != var_classes2)
        clashing_classes <- paste0(formula_vars[clashing_classes], collapse = ", ")
        
        stop(paste0("Variables ",
                    clashing_classes,
                    " have different classes"))
        
    }
    
    ## For factor variables, check the levels match
    factor_vars <- which(var_classes %in% c("ordered", "factor"))
    if (length(factor_vars) > 0) { 
        factor_vars <- formula_vars[factor_vars]
        
        l1 <- sapply(df1[,factor_vars], levels)
        l2 <- sapply(df2[,factor_vars], levels)
        test <- all.equal(l1, l2)

        if(!test) {
            stop(paste0("Variables have different levels: ",
                        paste0(test, collapse = "\n")))
        }
    }

    ## What about interactions?
    
    return(TRUE)
}

variables_in_df <- function(variables,
                            df,
                            label) {
    if (any(!is.element(variables, names(df)))) {
        missing_vars <- setdiff(variables, names(df))
        stop(paste0(label, " is missing variables: ",
                    paste0(missing_vars, collapse = ", ")))
    }
}
