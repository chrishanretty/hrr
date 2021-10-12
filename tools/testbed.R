library(devtools)
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
res2 <- res
res2$purple <- res2$red + res2$green
res2 <- unique(res2[,c("area",
                       "blue", "purple")])

### Create a weight var
toydata$w8 <- 1 + rnorm(nrow(toydata), sd = 0.1)

### What does the code look like for a multithreaded binomial model
toydata$vi2 <- ifelse(toydata$vi %in% c("red", "green"),
                     "purple",
                     "blue")

f <- vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
    cont1 + cont2
f2 <- vi2 ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
    cont1 + cont2

foo <- hrr(f, data = toydata, ps = toyps, result = res,
           areavar = "area",
           dirichlet = TRUE,
           algorithm = "meanfield",
           testing = FALSE)


bar <- hrr(f2, data = toydata, ps = toyps, result = res2,
           areavar = "area",
           dirichlet = FALSE,
           adjust = TRUE,
           algorithm = "sampling",
           testing = FALSE)

prelim_code <- make_stancode(f2, data = toydata, family = bernoulli,
                             backend = "cmdstanr",
                             threads= threading(2))


cat(hrr:::add_ps_data_code(prelim_code)[[1]]$scode, sep = "\n")

cat(hrr:::add_genquant_code(prelim_code, f2, toydata, adjust = FALSE, dirichlet = TRUE)[[1]]$scode)

### Can I do optimizing
### Write the code out


tf <- tempfile()
tf <- paste0(tf, ".stan")
writeLines(foo$code, con = tf)

mod <- cmdstan_model(tf)
standata <- list()

for (t in names(foo$data)) {
  standata[[t]] <- foo$data[[t]]
}

nReps <- 100
reps <- lapply(1:nReps, function(i)mod$optimize(data = standata, iter = 1e5))

### Are there draws?
repo <- lapply(reps, function(x) try(dim(x$draws())))
repc <- lapply(repo, class)

### 
pos <- which(unlist(repc) != "try-error")
successful_draws <- lapply(reps[pos], function(x)x$draws())
successful_draws <- do.call("rbind", successful_draws)
### Just params
successful_draws <- successful_draws[, grep("counts",
                                            colnames(successful_draws),
                                            invert = TRUE)]
successful_draws <- successful_draws[, grep("aggmu",
                                            colnames(successful_draws),
                                            invert = TRUE)]
### What's the variation like?
apply(successful_draws, 2, function(x) sd(x))

newton_fit <- mod$optimize(data = standata, seed = 456,
                     iter = 1e5,
                     algorithm = "newton")

lbfgs_fit <- mod$optimize(data = standata, seed = 456,
                     iter = 1e5,
                     algorithm = "lbfgs")
