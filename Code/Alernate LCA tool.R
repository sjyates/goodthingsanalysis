library("ca")
library("car")
library("corrplot")
library("DescTools")
library("depmixS4")
library("dplyr")
library("factoextra")
library("FactoMineR")
library("forcats")
library("foreign")
library("generalhoslem")
library("ggparallel")
library("ggplot2")
library("ggrepel")
library("gmodels")
library("haven")
library("Hmisc")
library("igraph")
library("knitr")
library("lattice")
library("logistf")
library("plyr")
library("poLCA")
library("RcmdrMisc")
library("reshape2")
library("sjlabelled")
library("tidyr")
library("xtable")

lca <- X2018_for_analysis

lca <- lca[, c("IN13A",
      "IN13B",
      "IN13C",
      "IN13D",
      "IN13E",
      "IN13F",
      "IN13G",
      "IN13H",
      "IN13I",
      "IN13J",
      "IN13K",
      "IN13L")]+1
mod1 <- mix(
  list(
    IN13A ~ 1,
    IN13B ~ 1,
    IN13C ~ 1,
    IN13D ~ 1,
    IN13E ~ 1,
    IN13F ~ 1,
    IN13G ~ 1,
    IN13H ~ 1,
    IN13I ~ 1,
    IN13J ~ 1,
    IN13K ~ 1,
    IN13L ~ 1
  ),
  data = lca,
  # the dataset to use
  nstates = 2,
  # the number of latent classes
  family = list(
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(),
    binomial(
    )
),
respstart = runif(24)
)

fmod1 <- fit(mod1, verbose = FALSE)
