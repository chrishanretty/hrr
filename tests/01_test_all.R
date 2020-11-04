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

