## Sink output to ovearll results file
sink("GoodThingsAnalysis", append = FALSE, split = TRUE)
pdf("OverallresutsGoogThings.pdf")

##Load necessary libraries
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
library("psych")
library("powerAnalysis")
library("RcmdrMisc")
library("logistf")
library("gmodels")
library("generalhoslem")
library("DescTools")
library("sjlabelled")
library("FactoMineR")
library("factoextra")
library("rms")
library("nnet")
library("effects")
library("jtools")
library("huxtable")
library("officer")
library("flextable")
library("descr")
library("xtable")

##Set up functions to be used in the analysis

##Set up function for entropy calculation
entropy <- function (var.p) {
  sum(-var.p * log(var.p))
}

overalllatex <- list()

## Function to undertake chisquare analayis and plot graphs
chisquaretest.predictions.function <-
  function(indfactor.data,
           predclass.data,
           indfactor.names,
           predclass.names) {
    chisquare.results <- chisq.test(indfactor.data, predclass.data)
    residuals.data <- chisquare.results$residuals
    colnames(residuals.data) <- predclass.names
    rownames(residuals.data) <- indfactor.names
    corrplot(residuals.data, is.cor = FALSE)
    contrib.data <-
      100 * residuals.data ^ 2 / chisquare.results$statistic
    round(contrib.data, 3)
    corrplot(contrib.data, is.cor = FALSE)
    return(chisquare.results)
  }

##Start of specific analyses of Ofcom data
##Create non-user and limited user variables

X2018.uses.data.naomit.lcaresults$Nonuser <- 0
X2018.uses.data.naomit.lcaresults$Nonuser[X2018.uses.data.naomit.lcaresults$Latent_Class == 7] <-
  1

X2018.uses.data.naomit.lcaresults$LimitedUser <- 0
X2018.uses.data.naomit.lcaresults$LimitedUser[X2018.uses.data.naomit.lcaresults$Latent_Class == 6] <-
  1

X2018.uses.data.naomit.lcaresults$LimitedUserSocialMedia <- 0
X2018.uses.data.naomit.lcaresults$LimitedUserSocialMedia[X2018.uses.data.naomit.lcaresults$Latent_Class == 5] <-
  1

X2018.uses.data.naomit.lcaresults$CombinedLimitedUser <- 0
X2018.uses.data.naomit.lcaresults$CombinedLimitedUser[X2018.uses.data.naomit.lcaresults$Latent_Class == 6] <-
  1
X2018.uses.data.naomit.lcaresults$CombinedLimitedUser[X2018.uses.data.naomit.lcaresults$Latent_Class == 5] <-
  1

X2018.uses.data.naomit.lcaresults$SocialMediaOnly <- 0
X2018.uses.data.naomit.lcaresults$SocialMediaOnly[X2018.uses.data.naomit.lcaresults$Latent_Class == 4] <-
  1

X2018.uses.data.naomit.lcaresults$NotVeryConfident <- 0
X2018.uses.data.naomit.lcaresults$NotVeryConfident[X2018.uses.data.naomit.lcaresults$IN11AD == 1] <-
  1
X2018.uses.data.naomit.lcaresults$NotVeryConfident[X2018.uses.data.naomit.lcaresults$IN11AE == 1] <-
  1
X2018.uses.data.naomit.lcaresults$NotVeryConfident[X2018.uses.data.naomit.lcaresults$IN11AF == 1] <-
  1

X2018.uses.data.naomit.lcaresults$HealthImpact <- 0
X2018.uses.data.naomit.lcaresults$HealthImpact[X2018.uses.data.naomit.lcaresults$C12 == 1] <-
  1

X2018.uses.data.naomit.lcaresults$Education <-
  X2018.uses.data.naomit.lcaresults$C7
X2018.uses.data.naomit.lcaresults$Education[X2018.uses.data.naomit.lcaresults$Education >=
                                              5] <- 1

X2018.uses.data.naomit.lcaresults$NumberChildren <-
  X2018.uses.data.naomit.lcaresults$C3

X2018.uses.data.naomit.lcaresults$Income <-
  X2018.uses.data.naomit.lcaresults$C9
X2018.uses.data.naomit.lcaresults$Income[X2018.uses.data.naomit.lcaresults$Income == 7] <-
  0

levels(X2018.uses.data.naomit.lcaresults$Latent_Class) <-     c(
  "Extensive political",
  "Extensive",
  "General (no social media)",
  "Social and media",
  "Social media limited",
  "Limited",
  "Non-users"
)


print("Binary logisitc regression on non-users")
print("---------------------------------------")
print("                                       ")
fit.glmlogit.nonuser <- glm(
  Nonuser ~
    Education +
    DEP +
    HealthImpact +
    HTYPE +
    QLOC +
    AWTV3 +
    AWTV2,
  family = binomial(link = "logit"),
  data = X2018.uses.data.naomit.lcaresults
)

summary(fit.glmlogit.nonuser)
PseudoR2(fit.glmlogit.nonuser, which = "McFadden")
PseudoR2(fit.glmlogit.nonuser, which = "Nagelkerke")
PseudoR2(fit.glmlogit.nonuser, which = "CoxSnell")
logitgof(
  fit.glmlogit.nonuser$y,
  fitted(fit.glmlogit.nonuser),
  g = 10,
  ord = FAlSE
)

fit.glmlogit.nonuser0 <- glm(Nonuser ~ 1,
                             family = binomial(link = "logit"), data = X2018.uses.data.naomit.lcaresults)

confint <- exp(cbind(coef(fit.glmlogit.nonuser),
                     confint(fit.glmlogit.nonuser)))
confint
overalllatex$model1confint <-
  xtable(confint, caption = "Non-user Binary GLM confidence intervals")

anova(fit.glmlogit.nonuser0,
      fit.glmlogit.nonuser,
      test = "Chisq")

X2018.uses.data.naomit.lcaresults$modelnonuser <- 0
X2018.uses.data.naomit.lcaresults$modelnonuser[fit.glmlogit.nonuser$linear.predict > 0] <-
  1
CrossTable(
  X2018.uses.data.naomit.lcaresults$Nonuser,
  X2018.uses.data.naomit.lcaresults$modelnonuser
)
accuracy <-
  table(
    X2018.uses.data.naomit.lcaresults$Nonuser,
    X2018.uses.data.naomit.lcaresults$modelnonuser
  )
sum(diag(accuracy)) / sum(accuracy)


model1 <- summ(fit.glmlogit.nonuser, confint = TRUE, digits = 3)
model1
overalllatex$model1 <-
  xtable(model1$coeftable, caption = "Nonusers Binary GLM")
export_summs(fit.glmlogit.nonuser,
             to.file = "pdf",
             file.name = "Nonuser.pdf")

#Subset of users only
print("Binary logisitc regression on limited users - only comparing to users")
print("---------------------------------------------------------------------")
print("                                       ")

X2018.uses.data.naomit.lcaresults.usersonly <-
  subset(X2018.uses.data.naomit.lcaresults, Latent_Class < 7)

fit.glmlogit.limited <- glm(
  CombinedLimitedUser ~
    Education +
    DEP +
    HealthImpact +
    HTYPE +
    QLOC +
    AWTV3 +
    AWTV2,
  family = binomial(link = "logit"),
  data = X2018.uses.data.naomit.lcaresults.usersonly
)

summary(fit.glmlogit.limited)
PseudoR2(fit.glmlogit.limited, which = "McFadden")
PseudoR2(fit.glmlogit.limited, which = "Nagelkerke")
PseudoR2(fit.glmlogit.limited, which = "CoxSnell")
logitgof(
  fit.glmlogit.limited$y,
  fitted(fit.glmlogit.limited),
  g = 10,
  ord = FAlSE
)

fit.glmlogit.limited0 <- glm(CombinedLimitedUser ~ 1,
                             family = binomial(link = "logit"),
                             data = X2018.uses.data.naomit.lcaresults.usersonly)

exp(cbind(coef(fit.glmlogit.limited),
          confint(fit.glmlogit.limited)))
anova(fit.glmlogit.limited0,
      fit.glmlogit.limited,
      test = "Chisq")

X2018.uses.data.naomit.lcaresults.usersonly$modellimited <- 0
X2018.uses.data.naomit.lcaresults.usersonly$modellimited[fit.glmlogit.limited$linear.predict > 0] <-
  1
CrossTable(
  X2018.uses.data.naomit.lcaresults.usersonly$CombinedLimitedUser,
  X2018.uses.data.naomit.lcaresults.usersonly$modellimited
)
accuracy <-
  table(
    X2018.uses.data.naomit.lcaresults.usersonly$CombinedLimitedUser,
    X2018.uses.data.naomit.lcaresults.usersonly$modellimited
  )
sum(diag(accuracy)) / sum(accuracy)

model2 <- summ(fit.glmlogit.limited, confint = TRUE, digits = 3)
overalllatex$model2 <-
  xtable(model2$coeftable, caption = "Combined Limited Users Binary GLM")
model2
export_summs(fit.glmlogit.limited,
             to.file = "pdf",
             file.name = "CombinedLimitedUser.pdf")

print("Binary logisitc regression on limited social media users - only comparing to users")
print(
  "----------------------------------------------------------------------------------"
)
print("                                       ")

fit.glmlogit.limitedsoc <- glm(
  LimitedUserSocialMedia ~
    Education +
    DEP +
    HealthImpact +
    HTYPE +
    QLOC +
    AWTV3 +
    AWTV2,
  family = binomial(link = "logit"),
  data = X2018.uses.data.naomit.lcaresults.usersonly
)

summary(fit.glmlogit.limitedsoc)
PseudoR2(fit.glmlogit.limitedsoc, which = "McFadden")
PseudoR2(fit.glmlogit.limitedsoc, which = "Nagelkerke")
PseudoR2(fit.glmlogit.limitedsoc, which = "CoxSnell")
logitgof(
  fit.glmlogit.limitedsoc$y,
  fitted(fit.glmlogit.limitedsoc),
  g = 10,
  ord = FAlSE
)

fit.glmlogit.limitedsoc0 <- glm(LimitedUserSocialMedia ~ 1,
                                family = binomial(link = "logit"),
                                data = X2018.uses.data.naomit.lcaresults.usersonly)

exp(cbind(
  coef(fit.glmlogit.limitedsoc),
  confint(fit.glmlogit.limitedsoc)
))
anova(fit.glmlogit.limitedsoc0,
      fit.glmlogit.limitedsoc,
      test = "Chisq")

X2018.uses.data.naomit.lcaresults.usersonly$modellimitedsoc <- 0
X2018.uses.data.naomit.lcaresults.usersonly$modellimitedsoc[fit.glmlogit.limitedsoc$linear.predict > 0] <-
  1
CrossTable(
  X2018.uses.data.naomit.lcaresults.usersonly$LimitedUserSocialMedia,
  X2018.uses.data.naomit.lcaresults.usersonly$modellimitedsoc
)
accuracy <-
  table(
    X2018.uses.data.naomit.lcaresults.usersonly$LimitedUserSocialMedia,
    X2018.uses.data.naomit.lcaresults.usersonly$modellimitedsoc
  )
sum(diag(accuracy)) / sum(accuracy)

model3 <- summ(fit.glmlogit.limitedsoc,
               confint = TRUE,
               digits = 3)
model3
overalllatex$model3 <-
  xtable(model3$coeftable, caption = "Limited Users (Social Media) Binary GLM")
export_summs(fit.glmlogit.limitedsoc,
             to.file = "pdf",
             file.name = "LimitedUserSocialMedia.pdf")

print("Binary logisitc regression on social media only users - only comparing to users")
print(
  "-------------------------------------------------------------------------------"
)
print("                                       ")

fit.glmlogit.socmedonly <- glm(
  SocialMediaOnly ~
    Education +
    DEP +
    HealthImpact +
    HTYPE +
    QLOC +
    AWTV3 +
    AWTV2,
  family = binomial(link = "logit"),
  data = X2018.uses.data.naomit.lcaresults
)

summary(fit.glmlogit.socmedonly)
PseudoR2(fit.glmlogit.socmedonly, which = "McFadden")
PseudoR2(fit.glmlogit.socmedonly, which = "Nagelkerke")
PseudoR2(fit.glmlogit.socmedonly, which = "CoxSnell")
logitgof(
  fit.glmlogit.socmedonly$y,
  fitted(fit.glmlogit.socmedonly),
  g = 10,
  ord = FAlSE
)

fit.glmlogit.socmedonly0 <- glm(SocialMediaOnly ~ 1,
                                family = binomial(link = "logit"),
                                data = X2018.uses.data.naomit.lcaresults)

exp(cbind(
  coef(fit.glmlogit.socmedonly),
  confint(fit.glmlogit.socmedonly)
))
anova(fit.glmlogit.socmedonly0,
      fit.glmlogit.socmedonly,
      test = "Chisq")

X2018.uses.data.naomit.lcaresults$modelsocmedonly <- 0
X2018.uses.data.naomit.lcaresults$modelsocmedonly[fit.glmlogit.socmedonly$linear.predict > 0] <-
  1
CrossTable(
  X2018.uses.data.naomit.lcaresults.usersonly$SocialMediaOnly,
  X2018.uses.data.naomit.lcaresults.usersonly$modelsocmedonly
)
accuracy <-
  table(
    X2018.uses.data.naomit.lcaresults.usersonly$SocialMediaOnly,
    X2018.uses.data.naomit.lcaresults.usersonly$modelsocmedonly
  )
printprint(sum(diag(accuracy)) / sum(accuracy))

model4 <- summ(fit.glmlogit.socmedonly,
               confint = TRUE,
               digits = 3)
model4
overalllatex$model4 <-
  xtable(model4$coeftable, caption = "Limited Users (Social Media) Binary GLM")
export_summs(fit.glmlogit.socmedonly,
             to.file = "pdf",
             file.name = "SocialMediaOnly.pdf")


print("Multinominal logisitc regression Latent Classes")
print("-----------------------------------------------")
print("                                       ")

fit.multinom <- multinom(Latent_Class ~
                           Education +
                           DEP +
                           HealthImpact +
                           HTYPE +
                           QLOC +
                           AWTV3 +
                           AWTV2,
                         data = X2018.uses.data.naomit.lcaresults)

print("Multinominal B values")
print("---------------------")
print("                                       ")

model5 <- summary(fit.multinom)
model5
overalllatex$model5B <-
  xtable(model5$coefficients, caption = "Multinomonal regression coefficients")
overalllatex$model5SE <-
  xtable(model5$standard.errors, caption = "Multinomonal regression standard errors")

print("Multinominal z scores")
print("---------------------")
print("                                       ")

zscores <-
  summary(fit.multinom)$coefficients / summary(fit.multinom)$standard.errors
zscores
overalllatex$model5Z <-
  xtable(zscores, caption = "Multinomonal regression zscores")

print("Multinominal p values")
print("---------------------")
print("                                       ")

pvalues <- (1 - pnorm(abs(z), 0, 1)) * 2
pvalues
overalllatex$model5Z <-
  xtable(pvalues, caption = "Multinomonal pvalues")

print("Multinominal odds ratios")
print("------------------------")
print("                                       ")

oddsratios <- exp(coef(fit.multinom))
oddsratios
overalllatex$model5odds <-
  xtable(pvalues, caption = "Multinomonal odds ratios")

#Exploring Chi Squares post Multinom

gvarlist = c("Education",
             "DEP",
             "HealthImpact",
             "HTYPE",
             "QLOC",
             "AWTV3",
             "S2XX",
             "AWTV2")

pvarlist <- list()
pvarlist$Education <-
  c("16 or under", "17-18", "19-20", "21 or over")
pvarlist$S2XX <- c("25-34",
                   "16-24",
                   "35-44",
                   "45-54",
                   "55-64",
                   "65-74",
                   "75+")
pvarlist$AWTV2 <- c("16-34",
                    "35-54",
                    "55+")
pvarlist$HealthImpact <- c("No Impact", "Impact")
pvarlist$HTYPE <- c(
  "1 adult",
  "2 adults",
  "3+ adults",
  "1 adult & children",
  "2 adults & children",
  "3+ adults & children"
)
pvarlist$Income <-
  c("No data",
    "< £10400",
    "to £15599",
    "to £25999",
    "to £36399",
    "to £51999",
    "> £51999")
pvarlist$DEP <- c("Low", "Middle", "High")
pvarlist$QLOC <- c("Urban", "Rural")
pvarlist$AWTV3 <- c("A&B", "C1", "C2", "D&E")

levels(X2018.uses.data.naomit.lcaresults$HealthImpact) <-
  c("No Impact", "Impact")

overallchiresults <-
  setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Analysis", "chi", "df", "p", "cv"))

limitedchiresults <-
  setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Analysis", "chi", "df", "p", "cv"))

for (gitem in gvarlist) {
  predclass.data <-
    haven::as_factor(X2018.uses.data.naomit.lcaresults[[gitem]])
  indfactor.data <-
    haven::as_factor(X2018.uses.data.naomit.lcaresults$Latent_Class)
  chicols <- pvarlist[[gitem]]
  chirows <-     c(
    "Extensive political",
    "Extensive",
    "General (no social media)",
    "Social and media",
    "Social media limited",
    "Limited",
    "Non-users"
  )
  levels(indfactor.data) <- chirows
  overallchitables[[gitem]] <-
    CrossTable(
      indfactor.data,
      predclass.data,
      prop.r = TRUE,
      prop.c = TRUE,
      prop.t = TRUE,
      prop.chisq = FALSE,
      chisq = TRUE,
      dnn = c("Latent Class", gitem),
      format = "SPSS",
      row.labels = TRUE
    )
  overalllatex[[gitem]] <-
    xtable(overallchitables[[gitem]], caption = paste(gitem, " vs Latent Class"))
  chiresults <-
    chisquaretest.predictions.function(predclass.data,
                                       indfactor.data,
                                       chicols,
                                       chirows)
  overallchiresults[nrow(overallchiresults) + 1,] = list(
    paste(gitem, " vs ", "Latent Class"),
    chiresults$statistic,
    chiresults$parameter,
    chiresults$p.value,
    cv.test(indfactor.data, predclass.data)
  )
}

for (gitem in gvarlist) {
  predclass.data <-
    haven::as_factor(X2018.uses.data.naomit.lcaresults[[gitem]])
  indfactor.data <-
    haven::as_factor(X2018.uses.data.naomit.lcaresults$CombinedLimitedUser)
  chicols <- pvarlist[[gitem]]
  chirows <-     c("Other Users",
                   "Limited Users")
  levels(indfactor.data) <- chirows
  overallchitables$limited[[gitem]] <-
    CrossTable(
      indfactor.data,
      predclass.data,
      prop.r = TRUE,
      prop.c = TRUE,
      prop.t = TRUE,
      prop.chisq = FALSE,
      chisq = TRUE,
      dnn = c("Latent Class", gitem),
      format = "SPSS",
      row.labels = TRUE
    )
  print(overallchitables$limited[[gitem]])
  overalllatex$limited[[gitem]] <-
    xtable(overallchitables$limited[[gitem]],
           caption = paste(gitem, " vs Limited Users"))
  chiresults <-
    chisquaretest.predictions.function(predclass.data,
                                       indfactor.data,
                                       chicols,
                                       chirows)
  limitedchiresults[nrow(overallchiresults) + 1,] = list(
    paste(gitem, " vs ", "Latent Class"),
    chiresults$statistic,
    chiresults$parameter,
    chiresults$p.value,
    cv.test(indfactor.data, predclass.data)
  )
}


print("Overall Chi-results for multiple comparisions table")
print("---------------------------------------------------")
print("                                       ")

overallchiresults
overalllatex$chitables <-
  xtable(overallchiresults, caption = "Overall Chi results")
savefile <- file.path(getwd(), "chiresults.csv")
write.csv(overallchiresults, savefile)

limitedchiresults
overalllatex$limitedchitables <-
  xtable(limitedchiresults, caption = "Overall Chi results")

overalllatex

dev.off()
sink()