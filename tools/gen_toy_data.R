
### Code to generate toy data (n = 800) and a toy post-stratification
### frame (n = 1e6) with binary and categorical outcomes.

library(devtools)
library(usethis)

set.seed(7113)
logit <- qlogis
inv.logit <- plogis
nPop <- 1e6
nObs <- 800
nAreas <- 26

area <- sample(LETTERS, size = nPop, replace = TRUE)

k_1 <- sample(letters[1:5],
              size = nPop,
              replace = TRUE,
              prob = rep(1/5, 5))

k_2 <- sample(letters[6:10],
              size = nPop,
              replace = TRUE,
              prob = c(0.5, 0.3, 0.1, 0.1, 0.1))

k_3 <- sample(letters[11:13],
              size = nPop,
              replace = TRUE,
              prob = c(0.8, 0.1, 0.1))

x_1 <- rnorm(nAreas)
x_2 <- rnorm(nAreas)

### Permute areas to match area in population data frame
x_1 <- x_1[match(area, LETTERS)]
x_2 <- x_2[match(area, LETTERS)]

### Group this in a data frame
psf <- data.frame(area = area,
                  k_1 = k_1,
                  k_2 = k_2,
                  k_3 = k_3,
                  x_1 = x_1,
                  x_2 = x_2)

for (v in c("k_1", "k_2", "k_3", "area")) {
    psf[, v] <- factor(psf[,v])
}

### Generate fake binary data
bin_alpha <- -0.5
bin_beta_1 <- rnorm(1)
bin_beta_2 <- rnorm(1)
bin_eta <- rnorm(nAreas)
bin_gamma_1 <- rnorm(length(unique(k_1)))
bin_gamma_2 <- rnorm(length(unique(k_2)), sd = 0.5)
bin_gamma_3 <- rnorm(length(unique(k_3)), sd = 2)

### Binary latent variable
psf$bin_mu <- with(psf,
                   bin_alpha +
                   bin_beta_1 * x_1 +
                   bin_beta_2 * x_2 +
                   bin_gamma_1[as.numeric(k_1)] +
                   bin_gamma_2[as.numeric(k_2)] +
                   bin_gamma_3[as.numeric(k_3)] +
                   bin_eta[as.numeric(area)])

psf$bin_y <- rbinom(nPop, size = 1, prob = inv.logit(psf$bin_mu))

### Generate fake multinomial data
cat_green_alpha <- 2/3
cat_blue_alpha <- -1/3
cat_green_beta_1 <- rnorm(1)
cat_blue_beta_1 <- rnorm(1)
cat_green_beta_2 <- rnorm(1)
cat_blue_beta_2 <- rnorm(1)
cat_green_eta <- rnorm(nAreas)
cat_blue_eta <- rnorm(nAreas)
cat_green_gamma_1 <- rnorm(length(unique(k_1)))
cat_green_gamma_2 <- rnorm(length(unique(k_2)))
cat_green_gamma_3 <- rnorm(length(unique(k_3)))
cat_blue_gamma_1 <- rnorm(length(unique(k_1)))
cat_blue_gamma_2 <- rnorm(length(unique(k_2)))
cat_blue_gamma_3 <- rnorm(length(unique(k_3)))

psf$cat_mu_red <- 0
psf$cat_mu_green <- with(psf,
                   cat_green_alpha +
                   cat_green_beta_1 * x_1 +
                   cat_green_beta_2 * x_2 +
                   cat_green_gamma_1[as.numeric(k_1)] +
                   cat_green_gamma_2[as.numeric(k_2)] +
                   cat_green_gamma_3[as.numeric(k_3)] +
                   cat_green_eta[as.numeric(area)])


psf$cat_mu_blue <- with(psf,
                   cat_blue_alpha +
                   cat_blue_beta_1 * x_1 +
                   cat_blue_beta_2 * x_2 +
                   cat_blue_gamma_1[as.numeric(k_1)] +
                   cat_blue_gamma_2[as.numeric(k_2)] +
                   cat_blue_gamma_3[as.numeric(k_3)] +
                   cat_blue_eta[as.numeric(area)])

mumat <- with(psf, cbind(exp(cat_mu_red), exp(cat_mu_green), exp(cat_mu_blue)))
draws <- t(apply(mumat, 1, rmultinom, n = 1, size = 1))
psf$cat_y <- apply(draws, 1, function(x)c("red", "green", "blue")[which(x == 1)])

### Sample from this to create the toy data
dat <- psf[sample(1:nPop, size = nObs, replace = FALSE),]

### Add on the actual results to the toydata frame
res <- tapply(psf$cat_y, list(area = psf$area), table, simplify = TRUE)
res <- do.call("rbind", res)
res <- data.frame(res)
res$area <- rownames(res)

dat <- merge(dat, res, by = "area", all = TRUE)

### Convert the data into a post-stratification frame with counts
psf$value <- 1
psf2 <- with(psf, aggregate(list("count" = value),
                            list("area" = area,
                                 "k_1" = k_1,
                                 "k_2" = k_2,
                                 "k_3" = k_3,
                                 "cat_y" = cat_y),
                            sum))
             
### Save these objects for future use
### The psf, the data, and the covariate vectors?
toydat <- dat
toypsf <- psf2

use_data(toydat, overwrite = TRUE)
use_data(toypsf, overwrite = TRUE)
