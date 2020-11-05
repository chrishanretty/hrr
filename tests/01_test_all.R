library(devtools)
build()
load_all()

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
### This takes 6 hours
### Drop down to threedepth = 14 and you get 41% hitting the max.
### but ESS is mostly okay. Saves you 20 minutes or so.
test <- hrr(vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
        cont1 + cont2,
    data = toydata,
    ps = toyps,
    result = res,
    areavar = "area",
    chains = 3,
    cores = 8,
    adjust = FALSE,
    threads_per_chain = 8,
    testing = FALSE,
    control = list(max_treedepth = 14,
                   adapt_delta = 0.95))

saveRDS(test, file = "../inst/extdata/toy_run.rds")

### Now do london run
data("londonmayor")
data("londonps")

### Check things are factorized properly
for (v in c("ageGroup", "education", "ethnicity", "gender")) {
    londonmayor[,v] <- factor(londonmayor[,v])
    londonps[,v] <- factor(londonps[,v])   
}

## Get results out
res <- unique(londonmayor[,c("ONSCode",
                             "Con_count_16",
                             "DNV_count_16",
                             "Green_count_16",
                             "Lab_count_16",
                             "LDem_count_16",
                             "Other_count_16",
                             "UKIP_count_16")])
names(res)[-1] <- levels(factor(londonmayor$vi))


londonps$days_after_elex <- 0
londonps$days_after_elex[1] <- 1
londonps <- merge(londonps,
                  unique(londonmayor[,c("ONSCode", "LabPct", "LDemPct", "GreenPct")]),
                  by = "ONSCode",
                  all.x = TRUE,
                  all.y = FALSE)

londonps$count <- as.integer(londonps$count)

test2 <- hrr(vi ~ (1|ONSCode) + (1|ageGroup) + (1|education) +
                 (1|ethnicity) +  gender + log1p(days_after_elex) + 
                 LabPct + LDemPct + GreenPct,
             data = londonmayor,
             ps = londonps,
             result = res,
             areavar = "ONSCode",
             chains = 3,
             cores = 8,
             adjust = FALSE,
             threads_per_chain = 8,
             testing = FALSE,
             control = list(max_treedepth = 15,
                            adapt_delta = 0.95))

saveRDS(test2, file = "../inst/extdata/toy_run.rds")
