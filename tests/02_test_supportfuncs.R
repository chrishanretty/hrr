### Purpose of this code
### To check whether we can successfully create stancode and standata
### using a variety of specifications (adjust, dirichlet, etc.,)
### on toydata and londonmayoral data

### Do the support functions work?
library(hrr)
library(dplyr)

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

f <- vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
    cont1 + cont2

fit1 <- hrr(f,
            data = toydata,
            ps = toyps,
            result = res,
            areavar = "area",
            chains = 3,
            cores = 8,
            adjust = FALSE,
            dirichlet = TRUE,
            testing = FALSE,
            algorithm = "meanfield")

g1 <- group_support(fit1)

aggs <- g1 %>%
    dplyr::group_by(grp_factor, grp_level) %>%
    dplyr::summarize(tots = sum(mean_share))

stopifnot(all(round(aggs$tots, 2) == 1))

fit2 <- hrr(f,
            data = toydata,
            ps = toyps,
            result = res,
            areavar = "area",
            chains = 3,
            cores = 8,
            adjust = FALSE,
            dirichlet = TRUE,
            testing = FALSE,
            iter = 800,
            warmup = 400,
            algorithm = "sampling")

g2 <- group_support(fit2)

aggs <- g2 %>%
    dplyr::group_by(grp_factor, grp_level) %>%
    dplyr::summarize(tots = sum(mean_share))

stopifnot(all(round(aggs$tots, 2) == 1))
