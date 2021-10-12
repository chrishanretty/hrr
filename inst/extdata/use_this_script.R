library(devtools)
library(usethis)

### Want to create two objects, londonmayor and londonps

ind <- readRDS("~/Dropbox/leverhulme/methods_paper/working/tidy_ind_data.rds")
ind$londonTurnoutW7 <- ind$londonFirstW7 <- NULL

res <- readRDS("~/Dropbox/leverhulme/methods_paper/working/tidy_results_all.rds")
names(res)[-1] <- paste0(names(res)[-1], "_counts_2016")

aux <- readRDS("~/Dropbox/leverhulme/methods_paper/working/tidy_aux_data.rds") %>%
    dplyr::select(ONSCode,
                  Ward,
                  Total,
                  ends_with("Pct"),
                  ends_with("Pct_mean"),
                  ends_with("Pct_sd"),
                  ends_with("Pct_sc"))


### Merge this together
londonmayor <- merge(ind, res, by = "ONSCode", all = TRUE)
londonmayor <- merge(londonmayor, aux, by = "ONSCode", all = TRUE)

londonmayor <- londonmayor %>%
    dplyr::filter(ONSCode != "E09000001")

londonps <- readRDS("~/Dropbox/leverhulme/methods_paper/working/raked_joint.rds")

usethis::use_data(londonmayor, overwrite = TRUE)
usethis::use_data(londonps, overwrite = TRUE)
