library(devtools)
build()
load_all()

### library(hrr)
data("toydata")
data("toyps")

toydata <- as.data.frame(toydata)
toyps <- as.data.frame(toyps)

### Tidy
fctrs <- c("area", "cat1", "cat2", "cat3")
for (f in fctrs) toydata[,f] <- factor(toydata[,f])
for (f in fctrs) toyps[,f] <- factor(toyps[,f])

### Create the results
res <- unique(toydata[,c("area",
                         "red", "green", "blue")])


### With adapt_delta = 0.95 and max_treedepth = 15,
### This takes 6 hours with no dirichlet sampling
### Drop down to threedepth = 14 and you get 41% hitting the max.
### but ESS is mostly okay. Saves you 20 minutes or so.
### With dirichlet sampling, then we're laughing, just 20 minutes
## test <- hrr(vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
##         cont1 + cont2,
##     data = toydata,
##     ps = toyps,
##     result = res,
##     areavar = "area",
##     chains = 3,
##     cores = 8,
##     adjust = FALSE,
##     threads_per_chain = 8,
##     testing = FALSE,
##     dirichlet = TRUE,
##     control = list(max_treedepth = 14,
##                    adapt_delta = 0.95))

test <- hrr(vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
        cont1 + cont2,
    data = toydata,
    ps = toyps,
    result = res,
    areavar = "area",
    chains = 3,
    cores = 8,
    adjust = FALSE,
    testing = FALSE,
    algorithm = "meanfield",
    dirichlet = TRUE)

### Do the area levels of support match the recovered levels of support?
yhat <- area_support(test)
y <- res %>%
    tidyr::pivot_longer(cols = -area, names_to = "outcome") %>%
    dplyr::group_by(area) %>%
    dplyr::mutate(known_share = value / sum(value))

plot.df <- merge(y, yhat,
                 by = c("area", "outcome"),
                 all = TRUE) %>%
    dplyr::mutate(col = ifelse((q5 > known_share | q95 < known_share),
                        "red",
                        "black"))

ggplot(plot.df, aes(x = known_share, y = mean_share,
                    ymin = q5, ymax = q95,
                    colour = col)) +
    geom_abline(slope = 1, intercept = 0) +
    scale_colour_identity(guide = FALSE) + 
    geom_pointrange()

### What about group support?
yhat <- group_support(test)

### saveRDS(test, file = "../inst/extdata/toy_run.rds")

## opt_test <- hrr(vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
##         cont1 + cont2,
##     data = toydata,
##     ps = toyps,
##     result = res,
##     areavar = "area",
##     chains = 3,
##     cores = 8,
##     adjust = FALSE,
##     threads_per_chain = 8,
##     testing = TRUE,
##     dirichlet = TRUE,
##     control = list(max_treedepth = 14,
##                    adapt_delta = 0.95))

## library(cmdstanr)
## code <- opt_test$code
## sf <- paste0(tempfile(), ".stan")
## writeLines(code, con = sf)
## mod <- cmdstan_model(sf)
## data <- opt_test$data
## class(data) <- "list"
## fit_optim <- mod$optimize(data = data,
##                           seed = 123,
##                           iter = 1000)

## fit_vb <- mod$variational(data = data,
##                        seed = 123,
##                        iter = 1000)

## s <- fit_optim$summary()

## ### Now do london run
## data("londonmayor")
## data("londonps")

## ### Check things are factorized properly
## for (v in c("ageGroup", "education", "ethnicity", "gender")) {
##     londonmayor[,v] <- factor(londonmayor[,v])
##     londonps[,v] <- factor(londonps[,v])   
## }

## ## Get results out
## res <- unique(londonmayor[,c("ONSCode",
##                              "Con_count_16",
##                              "DNV_count_16",
##                              "Green_count_16",
##                              "Lab_count_16",
##                              "LDem_count_16",
##                              "Other_count_16",
##                              "UKIP_count_16")])
## names(res)[-1] <- levels(factor(londonmayor$vi))


## londonps$days_after_elex <- 0
## londonps$days_after_elex[1] <- 1
## londonps <- merge(londonps,
##                   unique(londonmayor[,c("ONSCode", "LabPct", "LDemPct", "GreenPct")]),
##                   by = "ONSCode",
##                   all.x = TRUE,
##                   all.y = FALSE)

## londonps$count <- as.integer(londonps$count)

## ### With dirichlet = TRUE, we're talking just under 4 hours
## test2 <- hrr(vi ~ (1|ONSCode) + (1|ageGroup) + (1|education) +
##                  (1|ethnicity) +  gender + log1p(days_after_elex) + 
##                  LabPct + LDemPct + GreenPct,
##              data = londonmayor,
##              ps = londonps,
##              result = res,
##              areavar = "ONSCode",
##              chains = 3,
##              cores = 8,
##              adjust = FALSE,
##              threads_per_chain = 8,
##              testing = FALSE,
##              dirichlet = TRUE,
##              control = list(max_treedepth = 12,
##                             adapt_delta = 0.95))

## saveRDS(test2, file = "../inst/extdata/toy_run.rds")
