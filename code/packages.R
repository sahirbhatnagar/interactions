## ---- required-packages ----

getPckg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")

pckg = try(require(knitr))
if(!pckg) {
    cat("Installing 'knitr' from CRAN\n")
    getPckg("knitr")
    require(knitr)
}

pckg = try(require(data.table))
if(!pckg) {
    cat("Installing 'data.table' from CRAN\n")
    getPckg("data.table")
    require(data.table)
}

pckg = try(require(magrittr))
if(!pckg) {
    cat("Installing 'magrittr' from CRAN\n")
    getPckg("magrittr")
    require(magrittr)
}

pckg = try(require(glmnet))
if(!pckg) {
    cat("Installing 'glmnet' from CRAN\n")
    getPckg("glmnet")
    require(glmnet)
}

pckg = try(require(stringr))
if(!pckg) {
    cat("Installing 'stringr' from CRAN\n")
    getPckg("stringr")
    require(stringr)
}

pckg = try(require(plyr))
if(!pckg) {
    cat("Installing 'plyr' from CRAN\n")
    getPckg("plyr")
    require(plyr)
}

pckg = try(require(dplyr))
if(!pckg) {
    cat("Installing 'dplyr' from CRAN\n")
    getPckg("dplyr")
    require(dplyr)
}

pckg = try(require(gcdnet))
if(!pckg) {
    cat("Installing 'gcdnet' from CRAN\n")
    getPckg("gcdnet")
    require(gcdnet)
}



