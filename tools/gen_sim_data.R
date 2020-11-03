### Create some population data
### And sample from it
library(tidyverse)
set.seed(431)
npop <- 1e6
nobs <- 800
nAreas <- 16
ind <- data.frame(cat1 = sample(LETTERS[1:5], size = npop, replace = TRUE),
           cat2 = sample(LETTERS[6:10], size = npop,
                         prob = c(0.5, 0.3, 0.1, 0.1, 0.1), replace = TRUE),
           cat3 = sample(LETTERS[11:13], size = npop,
                         prob = c(0.8, 0.1, 0.1), replace = TRUE),
           area = sample(letters[1:nAreas], size = npop, replace = TRUE))

aux <- data.frame(area = letters[1:nAreas],
                  cont1 = rnorm(nAreas),
                  cont2 = rnorm(nAreas))

dat <- merge(ind, aux, by = "area", all = TRUE)
let2num <- function(x)match(toupper(x), LETTERS)
### Create mu
mu_cat1 <- rnorm(n = 13, mean = 0, sd = 1)
mu_cat2 <- rnorm(n = 13, mean = 0, sd = 2.5)
mu_cat3 <- rnorm(n = 13, mean = 0, sd = 5)
beta <- rnorm(n = 2)
eta <- rnorm(n = nAreas)

dat$mu <- mu_cat1[let2num(dat$cat1)] +
    mu_cat2[let2num(dat$cat2)] +
    mu_cat3[let2num(dat$cat3)] +
    beta[1] * dat$cont1 +
    beta[2] * dat$cont2 + 
    eta[let2num(dat$area)]

dat$vi <- cut(dat$mu,
              breaks = quantile(dat$mu, probs = c(0, 1/3, 2/3, 1)),
              labels = c("red", "green", "blue"))

dat$vi[is.na(dat$vi)] <- "red"

### Take from this the aggregate results
res <- dat %>%
    group_by(area, vi) %>%
    summarize(count = n()) %>%
    pivot_wider(id_cols = area, values_from = count, names_from = vi)

### Take from this the post-strat frame
ps <- dat %>%
    group_by(area, cat1, cat2, cat3) %>%
    summarize(count = n(),
              cont1 = unique(cont1),
              cont2 = unique(cont2))

### Take from this the individual level data
ind <- dat[sample(1:nrow(dat), size = nobs, replace = FALSE),]

save(res, ps, ind, file = "toy_data.rda")
