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
                             "Con_16_count",
                             "DNV_16_count",
                             "Green_16_count",
                             "Lab_16_count",
                             "LDem_16_count",
                             "Other_16_count",
                             "UKIP_16_count")])
names(res) <- c("ONSCode", "Con", "DNV", "Green",
                "Lab", "LDem", "Other", "UKIP")

compare_dfs(vi ~ ONSCode + ageGroup + education + ethnicity + gender +
           LabPct_sc + GreenPct_sc,
            londonmayor, londonps)

test <- hrr(vi ~ (1|ONSCode) + (1|ageGroup) + (1|education) + (1|ethnicity) + gender +
        LabPct_sc + GreenPct_sc,
    data = londonmayor,
    ps = londonps,
    result = res,
    areavar = "ONSCode",
    chains = 3,
    parallel_chains = 3,
    threads_per_chain = 4,
    testing = TRUE)
        
