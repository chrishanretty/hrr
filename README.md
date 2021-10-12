# The hrr package

`hrr` is an R package that estimates hierarchical related regression models. 

Hierarchical related regression models are described in:

> "Jackson, Christopher, And Nicky Best, and Sylvia Richardson. "Hierarchical related regression for combining aggregate and individual data in studies of socio‐economic disease risk factors." Journal of the Royal Statistical Society: Series A (Statistics in Society) 171.1 (2008): 159-178."

`hrr` is built on top of (R)Stan. The current package is experimental. Not all options are fully supported, and code may break regularly.

`hrr` is used in the following way:

```
library(hrr)
data("londonmayor")
data("londonps")

f <- vi ~ ageGroup + education + ethnicity + gender | LabPct_sc + LDemPct_sc + GreenPct_sc + OtherPct_sc

aux <- unique(londonmayor[,c("ONSCode", "LabPct_sc", "LDemPct_sc", "GreenPct_sc",
 "OtherPct_sc")])
res <- unique(londonmayor[,c("ONSCode", "Con_counts_2016", "Lab_counts_2016", "Green_counts_2016",
 "LDem_counts_2016", "UKIP_counts_2016", "Other_counts_2016", "DNV_counts_2016")])
names(res) <- c("ONSCode", "Con", "Lab", "Green",
"LDem", "UKIP", "Other", "DNV")


mod <- hrr(f, data = londonmayor, ps = londonps, aux = aux,
   res = res, areavar = "ONSCode", weightvar = "count",
   testing = FALSE, adjust = FALSE, overdispersed = TRUE,
   iter = 320, chains = 4, cores = 4)


```

