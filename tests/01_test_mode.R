### Purpose of this code
### To check whether we can successfully create stancode and standata
### using a variety of specifications (adjust, dirichlet, etc.,)
### on toydata

### Toy data
library(hrr)

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

### Create a weight var
toydata$w8 <- 1 + rnorm(nrow(toydata), sd = 0.1)

f <- vi ~ (1|area) + (1|cat1) + (1|cat2) + (1|cat3) + 
    cont1 + cont2

test_func <- function(testing, adjust, dirichlet, weights) {
    if (weights) { 
        hrr(f,
            data = toydata,
            ps = toyps,
            result = res,
            areavar = "area",
            chains = 3,
            cores = 8,
            adjust = adjust,
            dirichlet = dirichlet,
            weights = "w8",
            testing = testing,
            algorithm = "meanfield")
    } else {
        hrr(f,
            data = toydata,
            ps = toyps,
            result = res,
            areavar = "area",
            chains = 3,
            cores = 8,
            adjust = adjust,
            dirichlet = dirichlet,
            testing = testing,
            algorithm = "meanfield")
    }
    
}

test_inputs <- expand.grid(testing = c(FALSE, TRUE),
                           adjust = c(FALSE, TRUE),
                           dirichlet = c(FALSE, TRUE),
                           weights = c(FALSE, TRUE))


a <- test_func(TRUE, TRUE, FALSE, TRUE)

results <- apply(test_inputs, 1, function(x){
    print(x)
    test_func(testing = x["testing"],
              adjust = x["adjust"],
              dirichlet = x["dirichlet"],
              weights = x["weights"])
})

