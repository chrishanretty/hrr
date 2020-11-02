#' Estimate hierarchical related regression model
#'
#' @param formula a formula object, or something that can be coerced to a formula object
#' @param data data frame containing the individual level data
#' @param ps data frame containing the post-stratification data; must contain variable `count`
#' @param result data frame containing the results; must contain column names equal to levels of the dependent variable in `formula`
#' @param areavar a character containing the name of the variable giving the area
#' @return returns TRUE or reports an error
#'
hrr <- function(formula, data, ps, result) {
    return(TRUE)
}

hrr_code_func <- function(code, dat, depvar) {
### Purpose: pull together all the stanvars
### Input: code and data
    ### Output: stanvars
    require(brms)
    pll_to_pred_prob(code) +
        add_ps_data_code(code) +
        add_ps_tdata_code(levels(factor(dat[,depvar]))) +
        add_tpars_code(code) +
        add_model_code(code)
}

add_ps_data_code <- function(code) {
### Purpose: produce a stanvar object with data declarations for the post-strat. data
### Input: stan code
### Output: stanvar object
    require(stringr)
    require(brms)
### We need to get the existing data block
    data_block <- stringr::str_extract(code,
                              pattern = regex("data \\{.*\\}\ntransformed data \\{",
                                              multiline = TRUE,
                                              dotall = TRUE))

### Keep those lines which begin with X_, Z_, or J_
    data_block <- stringr::str_extract_all(data_block,
                              pattern = regex(".*(X_|Z_|J_).*",
                                              multiline = FALSE))[[1]]

### Replace these entries with ps_X_, etc.,
    data_block <- stringr::str_replace(data_block,
                              pattern = "(X_|Z_|J_)",
                              replacement = "ps_\\1")

### Replace references to N with ps_N
    data_block <- stringr::str_replace(data_block,
                              pattern = fixed("[N]"),
                              replacement = "[ps_N]")
    
    data_block <- stringr::str_replace(data_block,
                              pattern = fixed("[N, "),
                              replacement = "[ps_N, ")

### Add on the other things we need in this block;
### ps_N, nAreas, areastart, areastop
    top_block = c("
   int ps_N;
   int nAreas;
   int areastart[nAreas];
   int areastop[nAreas];
   int<lower=1> ps_counts[ps_N];
   int aggy[nAreas, ncat];
")

### Combine and return
    data_block <- c(top_block, data_block)
    data_block <- paste0(data_block, collapse = "\n")
    foo <- 1
    brms::stanvar(x = foo,
            scode = data_block,
            block = "data",
            position = "start")

}



add_ps_tdata_code <- function(cats) {
### Purpose: add data transformations to the transformed data block
### Input: list of output levels
### Output: stanvar
    require(glue)
    require(brms)
    cats <- cats[-1]
    if (any(grepl(" ", cats))) {
        stop("Don't yet know how to handle outcome levels with whitespace")
    }
    
    pscode_tdata_top <- glue::glue('
  matrix[ps_N, K_mu<<cats>> - 1] ps_Xc_mu<<cats>>;  // 
', .open = "<<", .close = ">>")
    
pscode_tdata_bottom <- glue::glue('
  for (i in 2:K_mu<<cats>>) {
    ps_Xc_mu<<cats>>[, i - 1] = ps_X_mu<<cats>>[, i] - means_X_mu<<cats>>[i - 1];
  } 
', .open = "<<", .close = ">>")

    brms::stanvar(pscode_tdata_top,
            scode = pscode_tdata_top,
            block = "tdata",
            position = "start") +
    brms::stanvar(pscode_tdata_top,
            scode = pscode_tdata_bottom,
            block = "tdata",
            position = "end")

}

pll_to_pred_prob<- function(code) {
### Purpose: create a new stan function
### Input: preliminary code
### output: stanvar
    require(stringr)
    require(brms)
    ## Get the function code
    code <- stringr::str_extract(code,
                pattern = regex("real partial_log_lik.*return ptarget;\\s+\\}",
                                dotall = TRUE,
                                multiline = TRUE))
### Amend the function declaration: name and return type
    code <- stringr::str_replace(code,
                        fixed("real partial_log_lik(int[] seq, "),
                        "vector pred_prob(")

### Amend the function declaration: additional argument
    code <- stringr::str_replace(code,
                        "vector pred_prob(.*?)\\) \\{",
                        "vector pred_prob\\1, int[] psweights) {")
    
### Replace the current declaration of the return variable
    code <- stringr::str_replace(code,
                        "real ptarget = 0;",
                        "vector[ncat] pp = rep_vector(0, ncat);")

### Amend the calculation
    code <- stringr::str_replace(code,
                        fixed("ptarget += categorical_logit_lpmf(Y[nn] | mu[n]);"),
                        "pp += softmax(mu[n]) * psweights[nn];")

### Amend the return value
    code <- stringr::str_replace(code,
                        fixed("return ptarget;"),
                        "return (pp / sum(pp));")

    brms::stanvar(scode = code,
            block = "functions")
} 

add_tpars_code <- function(code) {
### Purpose: add transformed parameters declaration
### Input: stan code
### Output: stanvar object
    
    top <- brms::stanvar(scode = "matrix[nAreas, ncat] aggmu;",
                   block = "tparameters",
                   position = "start")

    old_func <- stringr::str_extract(code,
                            pattern = "reduce_sum\\(partial_log_lik.*")

### Amend the beginning of the function call
    new_func <- stringr::str_replace(old_func,
                            fixed("reduce_sum(partial_log_lik, seq, grainsize, "),
                            "pred_prob(areastart[i], areastop[i], ")

### Amend the variable calls in the middle

    new_func <- stringr::str_replace_all(new_func,
                                fixed("Xc_"),
                                "ps_Xc_")
    
    new_func <- stringr::str_replace_all(new_func,
                                fixed("J_"),
                                "ps_J_")
    
    new_func <- stringr::str_replace_all(new_func,
                                fixed("Z_"),
                                "ps_Z_")

### Amend the end by adding psweights arg);
    new_func <- stringr::str_replace(new_func,
                        fixed(");"),
                        ", ps_counts)';")
    
    bottom <- paste0("
    for (i in 1:nAreas) {
        aggmu[i] = ",
    new_func,
    "
        }")

    bottom <- brms::stanvar(scode = bottom,
                      block = "tparameters",
                      position = "end")

    top + bottom
    
}

add_model_code <- function(code) {
### Purpose: add a sampling statement
### Input: preliminary code
### Output: stanvar
    scode <- "
  if (!prior_only) {
   for (i in 1:nAreas) { 
      aggy[i] ~ multinomial(to_vector(aggmu[i]));
   }
  }
"
    brms::stanvar(scode = scode,
            block = "model")
}

hrr_data_func <- function(formula, dat, ps, res, cats) { 
    prelim_data <- brms::make_standata(formula,
                          dat,
                          family = "categorical",
                          threads = threading(4))
    
    ps_data <- brms::make_standata(formula,
                             ps %>%
                             mutate(vi = sample(dat$vi,
                                                size = n(),
                                                replace = TRUE)),
                             family = "categorical",
                             threads = threading(4))

### Retain, from ps_data, anything which begins with an X_, a Z_, or or a J_
    new_ps_data <- ps_data[grep("X_|Z_|J_", names(ps_data))]
### Rename these variables
    names(new_ps_data) <- paste0("ps_", names(new_ps_data))
    new_ps_data$ps_N <- ps_data$N
    new_ps_data$nAreas <- length(unique(ps$ONSCode))
    areapos <- ps %>%
        dplyr::mutate(rowno = 1:dplyr::n()) %>%
        dplyr::group_by(ONSCode) %>%
        dplyr::summarize(start = min(rowno),
                  stop = max(rowno))

    new_ps_data$areastart <- areapos$start
    new_ps_data$areastop <- areapos$stop
    new_ps_data$ps_counts <- ps$count
### Combine with the other data
    data <- c(prelim_data, new_ps_data)
    data$prior_only <- 1

### Get aggregate results out
### Make sure they follow the same order as the ps data
    data$aggy <- as.matrix(res[,cats])
    return(data)
}

