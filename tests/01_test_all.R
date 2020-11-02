library(devtools)
build()
load_all()

data("londonmayor")
data("londonps")


### Tidy
fctrs <- c("ONSCode", "ageGroup", "education", "ethnicity", "gender", "vi")
for (f in fctrs) londonmayor[,f] <- factor(londonmayor[,f])

fctrs <- c("ONSCode", "ageGroup", "education", "ethnicity", "gender")
for (f in fctrs) londonps[,f] <- factor(londonps[,f])

### Add on some variables to the post-strat frame
londonps <- merge(londonps,
                  unique(londonmayor[,c("ONSCode", "LabPct_sc", "GreenPct_sc")]),
                  by = "ONSCode",
                  all.x = TRUE,
                  all.y = FALSE)

### Create the results
res <- unique(londonmayor[,c("ONSCode",
                             "Con_count_16",
                             "DNV_count_16",
                             "Green_count_16",
                             "Lab_count_16",
                             "LDem_count_16",
                             "Other_count_16",
                             "UKIP_count_16")])
names(res) <- c("ONSCode", "Con", "DNV", "Green",
                "Lab", "LDem", "Other", "UKIP")

### Set up the priors
library(brms)
### Start with one
bp <- set_prior("student_t(3, 0, 2.5", class = "Intercept", dpar = "muDNV")

### Continuous variables, plus gender
dpars <- paste0("mu", levels(factor(dat$vi))[-1])
vars <- c("gender", "LabPct_sc", "GreenPct_sc")
for (d in dpars) {
    for (v in vars) {
        bp <- c(bp,
                set_prior("normal(0, 1)",
                      class = "b",
                      coef = v,
                      dpar = d))
    }
}

### Categorical variables
vars <- c("ageGroup", "education", "ethnicity")
for (d in dpars) {
    for (v in vars) {
        bp <- c(bp,
                set_prior("normal(0, 5)",
                      class = "sd",
                      coef = v,
                      dpar = d))
    }
}

### Remove duplicates
bp <- bp[duplicated(bp),]

    
test <- hrr(vi ~ (1|ONSCode) + (1|ageGroup) + (1|education) + (1|ethnicity) + gender +
        LabPct_sc + GreenPct_sc,
    data = londonmayor,
    ps = londonps,
    result = res,
    areavar = "ONSCode",
    chains = 3,
    parallel_chains = 3,
    threads_per_chain = 4,
    testing = FALSE)
        
