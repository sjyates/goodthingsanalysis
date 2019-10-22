## load necessary libraries
library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
library("ggparallel")
library("igraph")
library("tidyr")
library("knitr")
library("xtable")
library("forcats")
library("lattice")
library("Hmisc")
library("corrplot")
library("haven")

new.probs.start <-
  poLCA.reorder(UserTypesLCAresults$best$probs.start, c(6, 5, 1, 2, 3, 4))

##run polca with adjusted ordering

UserTypesLCAresults$bestreordered <-
  poLCA(
    UsesBind,
    X2018.uses.data.naomit,
    nclass = 6,
    probs.start = new.probs.start,
    maxiter = 5000,
    tol = 1e-5,
    nrep = 20,
    na.rm = TRUE,
    verbose = TRUE,
    calc.se = TRUE
  )

graphsforlca.function(UserTypesLCAresults$bestreordered$probs,
                      "Internet user groups")

finalgraphforlca <- graphsforlca.function(UserTypesLCAresults$bestreordered$probs,
                                          "Internet user groups")
