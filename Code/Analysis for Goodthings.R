## load necessary libraries
library("corrplot")
library("DescTools")
library("dplyr")
library("factoextra")
library("FactoMineR")
library("forcats")
library("generalhoslem")
library("ggparallel")
library("ggplot2")
library("gmodels")
library("haven")
library("Hmisc")
library("igraph")
library("knitr")
library("lattice")
library("logistf")
library("plyr")
library("poLCA")
library("powerAnalysis")
library("psych")
library("RcmdrMisc")
library("reshape2")
library("sjlabelled")
library("tidyr")
library("xtable")


##Function to generate Cramers V as needed
cv.test = function(x, y) {
  CV = sqrt(chisq.test(x, y, correct = FALSE)$statistic /
              (length(x) * (min(
                length(unique(x)), length(unique(y))
              ) - 1)))
  print.noquote("Cramér V / Phi:")
  return(as.numeric(CV))
}

##Set up function for entropy calculation
entropy <- function (var.p) {
  sum(-var.p * log(var.p))
}

newdata <-
  cbind(X2018.uses.data.naomit$IOBS,
        UserTypesLCAresults$bestreordered$predclass)
colnames(newdata) <- c("IOBS", "Latent_Class")
X2018.uses.data.naomit.lcaresults <-
  merge(
    X2018_spss_processed_variable_names_GOODTHINGS,
    newdata,
    by = "IOBS",
    all.x = TRUE
  )
X2018.uses.data.naomit.lcaresults$Latent_Class[is.na(X2018.uses.data.naomit.lcaresults$Latent_Class)] <-
  7
write_sav(X2018.uses.data.naomit.lcaresults,
          "2018_data_with_user_lca.sav")

UserTypesLCAresults$Classchisq <-
  chisquaretest.predictions.function(
    X2018.uses.data.naomit.lcaresults$AWTV3,
    X2018.uses.data.naomit.lcaresults$Latent_Class,
    get_labels(X2018.uses.data.naomit.lcaresults$AWTV3),
    c(
      "Extensive political",
      "Extensive",
      "General (no social media)",
      "Social and media",
      "Social media limited",
      "Limited",
      "Non-users"
    )
  )

UserTypesLCAresults$Agechisq <-
  chisquaretest.predictions.function(
    X2018.uses.data.naomit.lcaresults$AWTV2,
    X2018.uses.data.naomit.lcaresults$Latent_Class,
    get_labels(X2018.uses.data.naomit.lcaresults$AWTV2),
    c(
      "Extensive political",
      "Extensive",
      "General (no social media)",
      "Social and media",
      "Social media limited",
      "Limited",
      "Non-users"
    )
  )

UserTypesLCAresults$SmartOnlychisq <-
  chisquaretest.predictions.function(
    X2018.uses.data.naomit.lcaresults$SUMDEVC,
    X2018.uses.data.naomit.lcaresults$Latent_Class,
    c("Uses laptop or desktop", "Phone or tablet only"),
    c(
      "Extensive political",
      "Extensive",
      "General (no social media)",
      "Social and media",
      "Social media limited",
      "Limited",
      "Non-users"
    )
  )

X2018.uses.data.naomit.lcaresults$C7 <-
  replace(
    X2018.uses.data.naomit.lcaresults$C7,
    X2018.uses.data.naomit.lcaresults$C7 == 5,
    NA
  )
X2018.uses.data.naomit.lcaresults$C7 <-
  replace(
    X2018.uses.data.naomit.lcaresults$C7,
    X2018.uses.data.naomit.lcaresults$C7 == 6,
    NA
  )


chisquaretest.predictions.function(
  X2018.uses.data.naomit.lcaresults.usersonly$SUMDEVC,
  X2018.uses.data.naomit.lcaresults.usersonly$Latent_Class,
  c("Uses laptop(+)", "Smart device only"),
  c(
    "Extensive political",
    "Extensive",
    "General (no social media)",
    "Social and media",
    "Social media limited",
    "Limited"
  )
)

X2018.uses.data.naomit.lcaresults.usersonly$TotalDevices <-
  X2018.uses.data.naomit.lcaresults.usersonly$IN2A +
  X2018.uses.data.naomit.lcaresults.usersonly$IN2B +
  X2018.uses.data.naomit.lcaresults.usersonly$IN2C +
  X2018.uses.data.naomit.lcaresults.usersonly$IN2D +
  X2018.uses.data.naomit.lcaresults.usersonly$IN2E +
  X2018.uses.data.naomit.lcaresults.usersonly$IN2F +
  X2018.uses.data.naomit.lcaresults.usersonly$IN2G

X2018.uses.data.naomit.lcaresults.usersonly$JustSmartphone <- 0
X2018.uses.data.naomit.lcaresults.usersonly$JustSmartphone[X2018.uses.data.naomit.lcaresults.usersonly$IN2A == 1 &
                                                             X2018.uses.data.naomit.lcaresults.usersonly$TotalDevices == 1] <-
  1

chisquaretest.predictions.function(
  X2018.uses.data.naomit.lcaresults.usersonly$JustSmartphone,
  X2018.uses.data.naomit.lcaresults.usersonly$Latent_Class,
  c("Multi-device", "Smart phone only"),
  c(
    "Extensive political",
    "Extensive",
    "General (no social media)",
    "Social and media",
    "Social media limited",
    "Limited"
  )
)
