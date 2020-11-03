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
    
test <- hrr(vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
        cont1 + cont2,
    data = toydata,
    ps = toyps,
    result = res,
    areavar = "area",
    chains = 3,
    adjust = TRUE,
    parallel_chains = 3,
    threads_per_chain = 8,
    testing = TRUE)

sf <- paste0(tempfile(), ".stan")
writeLines(test[[1]]$code(), con = sf)

mod <- cmdstanr::cmdstan_model(sf,
                               cpp_options = list(stan_threads = TRUE))

out <- mod$sample(data = test[[2]],
                  parallel_chains = 3,
                  chains = 3,
                  threads_per_chain = 8)

### This takes about ten minutes.
test <- hrr(vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
        cont1 + cont2,
    data = toydata,
    ps = toyps,
    result = res,
    areavar = "area",
    adjust = TRUE,
    chains = 3,
    parallel_chains = 3,
    threads_per_chain = 8,
    iter_warmup = 250,
    iter_sampling = 250,
    init = 0,
    testing = FALSE)

s <- test$summary()
aggs <- test$summary("aggmu")

### Do these match what we know the result to be?
aggs <- aggs %>%
    tidyr::separate(variable, into = c("var1", "var2"), sep = ",") %>%
    dplyr::mutate(var1 = gsub("[^0-9]", "", var1),
           var1 = as.numeric(var1),
           var2 = gsub("[^0-9]", "", var2),
           var2 = as.numeric(var2)) %>%
    dplyr::mutate(ONSCode = levels(londonps$ONSCode)[var1],
                  party = levels(londonmayor$vi)[var2]) %>%
    dplyr::select(ONSCode, party, mean, q5, q95)

### Merge this with the results
res[,-1] <- res[,-1] / rowSums(res[,-1])
res <- res %>%
    tidyr::pivot_longer(cols = -ONSCode,
                        names_to = "party")

comb <- merge(aggs, res,
              by = c("ONSCode", "party"),
              all = TRUE)
