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

##Set up functions to be used in the analysis

##Set up function for entropy calculation
entropy <- function (var.p) {
  sum(-var.p * log(var.p))
}

##Function to find best result from 20 iterations
##to protect agsinst finding local minima
##Returns best LCA results and BIC values table
findbestlca.function <-
  function(analysis.vars,
           analysis.data,
           var.minclasses,
           var.maxclasses) {
    ##Set up temporary return, results and models frames
    bestnoclasses.var = 0
    lcaresults.list <- list()
    lcmodels.dataframe <- data.frame(
      NoClasses = 0,
      ll = 0,
      df = 0,
      BIC = 0,
      AIC = 0,
      ll_ratio = 0,
      Chi = 0,
      Entropy = 0
    )
    ##Loop to undertake LCA on min to max classes
    for (i in seq(var.minclasses, var.maxclasses, 1)) {
      str.message <- paste("Running model for ", i, "classes")
      writeLines(str.message)
      lc <-
        poLCA(
          analysis.vars,
          analysis.data,
          nclass = i,
          maxiter = 5000,
          tol = 1e-5,
          nrep = 20,
          na.rm = TRUE,
          verbose = TRUE,
          calc.se = TRUE
        )
      ##Check for best model based on BIC criterion
      #If this is first model set min_bic to model BIC
      if (min_bic == 0) {
        min_bic <- lc$bic
      }
      ##If only one model sought then dont worry about min_bic, else check
      if (var.minclasses == var.maxclasses) {
        lca.best.model <- lc
        bestnoclasses.var = i
      } else if (lc$bic < min_bic) {
        min_bic <- lc$bic
        lca.best.model <- lc
        bestnoclasses.var = i
      }
      ##Store criteria for comparative table
      lcmodels.dataframe[i - 1, 1] <- i
      lcmodels.dataframe[i - 1, 2] <- lc$llik
      lcmodels.dataframe[i - 1, 3] <- lc$resid.df
      lcmodels.dataframe[i - 1, 4] <- lc$bic
      lcmodels.dataframe[i - 1, 5] <- lc$aic
      lcmodels.dataframe[i - 1, 6] <- lc$Gsq
      lcmodels.dataframe[i - 1, 7] <- lc$Chisq
      if (i > var.minclasses) {
        error_post <- mean(apply(lc$posterior, 1, entropy), na.rm = TRUE)
        lcmodels.dataframe[i - 1 , 8] <-
          round(((error_prior - error_post) / error_prior), 3)
        error_prior <- entropy(lc$P) # class proportions model 2
      } else {
        lcmodels.dataframe[i, 8] <- c("-")
        error_prior <- entropy(lc$P) # class proportions model 1
      }
    }
    ##Force consistent order for classes
    probs.start.new <-
      poLCA.reorder(lca.best.model$probs.start,
                    order(lca.best.model$P, decreasing = TRUE))
    lca.best.model <- poLCA(
      analysis.vars,
      analysis.data,
      nclass = bestnoclasses.var,
      maxiter = 5000,
      tol = 1e-5,
      na.rm = TRUE,
      verbose = FALSE,
      calc.se = TRUE,
      probs.start = probs.start.new
    )
    
    ##Create return list
    lcaresults.list$models <- lcmodels.dataframe
    lcaresults.list$best <- lca.best.model
    lcaresults.list$bestnoclasses <- bestnoclasses.var
    writeLines("\nTable of BIC values for model selection\n")
    print(lcaresults.list$models)
    writeLines("\nBest Model From Analysis\n")
    print(lcaresults.list$best)
    return(lcaresults.list)
  }

##Function to plot variable probabilites by latent class
graphsfindbestlca.function <-
  function(lcmodel, str.title, str.fill, str.x) {
    #melt(lcaResults, level = 2)
    lcaresult.horizontal.plot <-
      ggplot(lcmodel, aes(x = L2, y = value, fill = X2)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(X1 ~ .) +
      scale_fill_brewer(type = "seq", palette = "Blues") +
      theme_bw() +
      ggtitle(str.title) +
      labs(fill = str.fill, x = str.x, y = "Share of item-\nresponse categories") +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        panel.grid.major.y = element_blank()
      ) +
      guides(fill = guide_legend(reverse = TRUE))
    print(lcaresult.horizontal.plot)
    lcaresult.vertical.plot <-
      ggplot(lcmodel, aes(x = X1, y = value, fill = X2)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap( ~ L2) +
      scale_x_discrete(str.x, expand = c(0, 0)) +
      scale_y_continuous("Share of item-\nresponse categories", expand = c(0, 0)) +
      scale_fill_brewer(type = "seq", palette = "Blues") +
      theme_bw() +
      ggtitle(str.title) +
      labs(fill = str.fill) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        panel.grid.major.y = element_blank()
      ) +
      guides(fill = guide_legend(reverse = TRUE))
    print(lcaresult.vertical.plot)
  }

##Funciton to plot regresssion lines where LCA has been performed againts predictors
regressiongraph.function <-
  function(analysis.data,
           str.title,
           str.xaxistitle,
           var.range) {
    pidmat <- cbind(1, c(1:var.range))
    exb <- exp(pidmat %*% analysis.data)
    matplot(
      c(1:var.range),
      (cbind(1, exb) / (1 + rowSums(exb))),
      main = str.title,
      xlab = str.xaxistitle,
      ylab = "Probability of latent class membership",
      ylim = c(0, 1),
      type = "l",
      lwd = 3,
      col = 1:3
    )
    legend(
      "topleft",
      inset = .05,
      legend = c(str.c1, str.c2, str.c3),
      pch = 1,
      col = 1:3,
      horiz = TRUE
    )
  }

##Function to plot BIC and relted values
plotbicandentropy.function <- function(models) {
  results.melt <-
    melt(models, id = "NoClasses")
  results.melt
  models.plot <- ggplot(results.melt) +
    geom_point(aes(x = NoClasses, y = value), size = 2) +
    geom_line(aes(NoClasses, value, group = 1)) +
    theme_bw() +
    labs(x = "", y = "", title = "") +
    facet_grid(variable ~ ., scales = "free") +
    theme_bw(base_size = 8, base_family = "") +
    theme(
      panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line(colour = "grey", size = 0.5),
      legend.title = element_text(size = 8, face = 'bold'),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text =  element_text(size = 8),
      axis.line = element_line(colour = "black"),
      strip.text.y = element_text(size = 8)
    )
  print(models.plot)
}

##Function to create regression results table where regresssion has been selected
regressionresults.function <- function(lcamodel, str.regname) {
  var.numberclasses <- length(lcamodel$P)
  var.checkforregression <-
    ifelse(is.na(lcamodel$coeff[1]), 1, nrow(lcamodel$coeff))
  
  if (var.checkforregression > 1) {
    i = 0
    for (r in 2:var.numberclasses) {
      i = i + 1
      regression.table <- data.frame(
        ##reg =  paste(r, "/ 1"),
        coeff = round(lcamodel$coeff[, (r - 1)], 5),
        se = round(lcamodel$coeff.se[, (r - 1)], 5),
        tval = round(lcamodel$coeff[, (r - 1)] / lcamodel$coeff.se[, (r -
                                                                        1)], 3),
        pr = round(1 - (2 * abs(
          pt(lcamodel$coeff[, (r - 1)] / lcamodel$coeff.se[, (r - 1)], lcamodel$resid.df) - 0.5
        )), 3)
      )
      str.filename <-
        paste(getwd(),
              "Regression",
              str.regname,
              r,
              "-1.latex",
              sep = "")
      print(xtable(regression.table), file = str.filename)
      str.filename <-
        paste(getwd(),
              "Regression",
              str.regname,
              r,
              "-1.htm",
              sep = "")
      print(xtable(regression.table),
            type = "html",
            file = str.filename)
    }
    
  }
}

## Function to undertake chisquare analayis and plot graphs of reisduals and contributions
chisquaretest.predictions.function <-
  function(indfactor.data,
           predclass.data,
           noclasses,
           pitem,
           gitem,
           chirows,
           chicols) {
    chisquare.results <- chisq.test(indfactor.data, predclass.data)
    residuals.data <- chisquare.results$residuals
    colnames(residuals.data) <- chicols
    rownames(residuals.data) <- chirows
    title.text <-
      paste(
        "Residuals: chi Square Crosstabulation of\n",
        pitem,
        "and",
        gitem,
        "\n(Chisquare =",
        round(chisquare.results$statistic, 3),
        " p <",
        round(chisquare.results$p.value, 3),
        ")",
        sep = " "
      )
    corrplot(
      residuals.data,
      is.cor = FALSE,
      title = title.text,
      mar = c(0, 0, 4, 0)
    )
    contrib.data <-
      100 * residuals.data ^ 2 / chisquare.results$statistic
    round(contrib.data, 3)
    colnames(contrib.data) <- chicols
    rownames(contrib.data) <- chirows
    title.text <-
      title.text <-
      paste(
        "Contributions: chi Square Crosstabulation of\n",
        pitem,
        "and",
        gitem,
        "\n(Chisquare =",
        round(chisquare.results$statistic, 3),
        " p <",
        round(chisquare.results$p.value, 3),
        ")",
        sep = " "
      )
    corrplot(
      contrib.data,
      is.cor = FALSE,
      title = title.text,
      mar = c(0, 0, 4, 0)
    )
    return(chisquare.results)
  }


##Function to generate Cramers V as needed
cv.test = function(x, y) {
  CV = sqrt(chisq.test(x, y, correct = FALSE)$statistic /
              (length(x) * (min(
                length(unique(x)), length(unique(y))
              ) - 1)))
  print.noquote("Cramér V / Phi:")
  return(as.numeric(CV))
}

##Start of specific analyses of Ofcom data
X2017_data_variables <-
  read_sav("SPSS.sav")

X2017_data_variables$Age[is.na(X2017_data_variables$Age)] <- 8

fit.glmlogit.nonuser <- glm(
  Nonuser ~
    Education +
    Age +
    HealthImpact +
    NumberChildrenByAge +
    NotVeryConfident +
    Income +
    DEP +
    QLOC +
    AWTV3,
  family = binomial(link = "logit"),
  data = X2017_data_variables
)

summary(fit.glmlogit.nonuser)
PseudoR2(fit.glmlogit.nonuser, which = "Nagelkerke")
PseudoR2(fit.glmlogit.nonuser, which = "CoxSnell")
logitgof(fit.glmlogit.nonuser$y,
         fitted(fit.glmlogit.nonuser),
         g = 10)

fit.glmlogit.nonuser0 <- glm(Nonuser ~ 1,
                             family = binomial(link = "logit"), data = X2017_data_variables)

exp(cbind(coef(fit.glmlogit.nonuser),
          confint(fit.glmlogit.nonuser)))
anova(fit.glmlogit.nonuser0,
      fit.glmlogit.nonuser,
      test = "Chisq")

X2017_data_variables$modelglmbin <- 0
X2017_data_variables$modelglmbin[fit.glmlogit.nonuser$linear.predict > 0] <-
  1
CrossTable(X2017_data_variables$Nonuser,
           X2017_data_variables$modelglmbin)
accuracy <-
  table(X2017_data_variables$Nonuser,
        X2017_data_variables$modelglmbin)
sum(diag(accuracy)) / sum(accuracy)

fit.glmlogit.costbarriers <- logistf(
  CostBarriers ~
    Education +
    Age +
    HealthImpact +
    NumberChildrenByAge +
    NotVeryConfident +
    Income +
    DEP +
    QLOC +
    AWTV3,
  family = binomial(link = "logit"),
  data = X2017_data_variables
)

summary(fit.glmlogit.costbarriers)

exp(cbind(
  coef(fit.glmlogit.costbarriers),
  confint(fit.glmlogit.costbarriers)
))

X2017_data_variables$modelglmcost <- 0
X2017_data_variables$modelglmcost[fit.glmlogit.costbarriers$predict > 0.5] <-
  1
CrossTable(X2017_data_variables$CostBarriers,
           X2017_data_variables$modelglmcost)
accuracy <-
  table(X2017_data_variables$CostBarriers,
        X2017_data_variables$modelglmcost)
sum(diag(accuracy)) / sum(accuracy)

fit.glmlogit.fearbarriers <- logistf(
  FearBarriers ~
    Education +
    Age +
    HealthImpact +
    NumberChildrenByAge +
    NotVeryConfident +
    Income +
    DEP +
    QLOC +
    AWTV3,
  family = binomial(link = "logit"),
  data = X2017_data_variables
)

summary(fit.glmlogit.fearbarriers)

exp(cbind(
  coef(fit.glmlogit.fearbarriers),
  confint(fit.glmlogit.fearbarriers)
))

X2017_data_variables$modelglmfear <- 0
X2017_data_variables$modelglmfear[fit.glmlogit.fearbarriers$predict > 0.5] <-
  1
CrossTable(X2017_data_variables$FearBarriers,
           X2017_data_variables$modelglmfear)
accuracy <-
  table(X2017_data_variables$FearBarriers,
        X2017_data_variables$modelglmfear)
sum(diag(accuracy)) / sum(accuracy)

gvarlist = c(
  "Education",
  "Age",
  "HealthImpact",
  "NumberChildrenByAge",
  "NotVeryConfident",
  "Income",
  "DEP",
  "QLOC",
  "AWTV3"
)

pvarlist$Education <- c("16 or under", "17-18", "19-20", "21 or over")
pvarlist$Age <-
  c("16-25",
    "25-35",
    "36-45",
    "46-55",
    "56-65",
    "66-75",
    "76-85",
    "85+")
pvarlist$AgeSimple <-
  c("16-35",
    "36-55",
    "56-75",
    "76+")
pvarlist$HealthImpact <- c("No Impact", "Impact")
pvarlist$Children<- c("No Children", "Children") 
pvarlist$NumberChildrenByAge <-
  c(levels(X2017_data_variables$NumberChildrenByAge))
pvarlist$NotVeryConfident <-
  c("Very Confident", "Not Very Confident")
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
pvarlist$CostBarriers <- c("No Cost", "Cost")
pvarlist$FearBarriers <- c("No Fear", "Fear")
pvarlist$EquipSupportBarriers <- c("No Equip", "Equip")
pvarlist$ComplicationBarriers <- c("No Complexity", "Complexity")

overallchiresults <-
  setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Analysis", "chi", "df", "p", "cv"))

for (gitem in gvarlist) {
  predclass.data <- X2017_data_variables[[gitem]]
  indfactor.data <- X2017_data_variables$FearBarriers
  chicols <- pvarlist[[gitem]]
  chirows <- c("No Fear Barriers", "Fear Barriers")
  chiresults <-
    chisquaretest.predictions.function(
      indfactor.data,
      predclass.data,
      length(chicols),
      "Fear Barriers",
      gitem,
      chirows,
      chicols
    )
  overallchiresults[nrow(overallchiresults) + 1,] = list(
    paste(gitem, " vs ", pitem),
    chiresults$statistic,
    chiresults$parameter,
    chiresults$p.value,
    cv.test(indfactor.data, predclass.data)
  )
}

savefile <- file.path(getwd(), "chiresults.csv")
write.csv(overallchiresults, savefile)

X2017_data_variables$Children <- 0
X2017_data_variables$Children[X2017_data_variables$NumberChildrenByAge > 0] <-
  1
X2017_data_variables$AgeSimple <-
  binVariable(X2017_data_variables$Age, bins = 4, "natural")
# X2017_data_variables$HealthImpact <- X2017_data_variables$HealthImpact+1
# X2017_data_variables$Children <- X2017_data_variables$Children+1
# X2017_data_variables$NotVeryConfident <- X2017_data_variables$NotVeryConfident+1
# X2017_data_variables$Nonuser <- X2017_data_variables$Nonuser+1
# X2017_data_variables$CostBarriers <- X2017_data_variables$CostBarriers+1
# X2017_data_variables$FearBarriers <- X2017_data_variables$FearBarriers+1
# X2017_data_variables$EquipSupportBarriers <- X2017_data_variables$EquipSupportBarriers+1
# X2017_data_variables$ComplicationBarriers <- X2017_data_variables$ComplicationBarriers+1
  
mca.variables <- c(
  "Education",
  #"AgeSimple",
  #"HealthImpact",
  "Children",
  "NotVeryConfident",
  #"Income",
  #"DEP",
  #"QLOC",
  #"AWTV3",
  #"Nonuser",
  "CostBarriers",
  "FearBarriers",
  "EquipSupportBarriers",
  "ComplicationBarriers"
)

lca.variables <- cbind(
  Education,
  AgeSimple,
  HealthImpact,
  Children,
  NotVeryConfident,
  #Income,
  #DEP,
  #QLOC,
  AWTV3,
  #Nonuser,
  CostBarriers,
  FearBarriers,
  EquipSupportBarriers,
  ComplicationBarriers
) ~ 1

mca.data <- X2017_data_variables[, mca.variables]

# for (item in mca.variables) {
#   mca.data[[item]] <- as_label(X2017_data_variables[[item]])
# }

for (item in mca.variables) {
  mca.data[[item]] <- as.factor(X2017_data_variables[[item]])
  print(item)
  print(levels(mca.data[[item]]))
  levels(mca.data[[item]]) <- pvarlist[[item]]
}


cats = apply(mca.data, 2, function(x)
  nlevels(as.factor(x)))

mca.analysis <- MCA(mca.data, graph = FALSE)

# data frame with variable coordinates
mca.vars.df <-
  data.frame(mca.analysis$var$coord, Variable = rep(names(cats), cats))

# data frame with observation coordinates
mca.obs.df <- data.frame(mca.analysis$ind$coord)

# plot of variable categories
ggplot(data = mca.vars.df,
       aes(
         x = Dim.1,
         y = Dim.2,
         label = rownames(mca.vars.df)
       )) +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_text(aes(colour = Variable)) +
  ggtitle("MCA plot of variables using R package FactoMineR")

ggplot(data = mca.obs.df, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  geom_point(colour = "gray50", alpha = 0.7) +
  geom_density2d(colour = "gray80") +
  geom_text(data = mca.vars.df,
            aes(
              x = Dim.1,
              y = Dim.2,
              label = rownames(mca.vars.df),
              colour = Variable
            )) +
  ggtitle("MCA plot of variables using R package FactoMineR") +
  scale_colour_discrete(name = "Variable")

fviz_mca_var(
  mca.analysis,
  col.var = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE,
  # Avoid text overlapping
  ggtheme = theme_minimal()
)

lca.data.users <- X2017_data_variables

lca.var.users.recode <- c(
  "IN13D",
  "IN13J",
  "IN13H",
  "IN13C",
  "IN13E",
  "IN13G",
  "IN13I",
  "IN13F",
  "IN13A",
  "IN21",
  "G1I",
  "IN22H",
  "IN13B",
  "IN22F",
  "IN22I",
  "IN22C",
  "IN19A",
  "T2" 
)

for (item in lca.var.users.recode) {
  lca.data.users[[item]] <-
    recode(lca.data.users[[item]], "0=1;1=2;NA=1")
  print(table(lca.data.users[[item]]))
}

lca.data.users <-
  plyr::rename(
    lca.data.users,
    c("IN13D" = "Pay_online_for_local_service",
"IN13J"	= "Pay_bills_online",
"IN13H"	= "Compare_products",
"IN13C"	= "Complete_Government_processes",
"IN13E"	= "Look_online_for_public_services",
"IN13G"	= "Find_information_for_leisure",
"IN13I"	= "Find_information_online_about_cultural",
"IN13F"	= "Look_online_job_opportunities",
"IN13A"	= "Access_news_websites_politics",
"IN21" = "Use_social_media",
"G1I" = "Play_games",
"IN22H"	= "Twitter",
"IN13B"	= "Sign_an_online_petition",
"IN22F"	= "Snapchat",
"IN22I"	= "WhatsApp",
"IN22C"	= "Instagram",
"IN19A"	= "Watch_videos_Like_YouTube",
"T2" = "On_demand_TV"
))


lca.var.users <- cbind(
  Pay_online_for_local_service,
  Pay_bills_online,
  Compare_products,
  Complete_Government_processes,
  Look_online_for_public_services,
  Find_information_for_leisure,
  Find_information_online_about_cultural,
  Look_online_job_opportunities,
  Access_news_websites_politics,
  Use_social_media,
  Play_games,
  Twitter,
  Sign_an_online_petition,
  Snapchat,
  WhatsApp,
  Instagram,
  Watch_videos_Like_YouTube,
On_demand_TV) ~ 1

lca.users <-
  findbestlca.function(lca.var.users,
                       lca.data.users,
                       2,
                       10)

lca.users.bestprobs <-
  spread(melt(lca.reasons$best$probs, level = 2),
         "L2",
         "value")
savefile <- file.path(getwd(), "lcabestprobs.csv")
write.csv(lca.reasons.bestprobs, savefile)

##Lattice plot of variable proportions by class
title <-
  paste("Reasons variables", "\nNo Regression", sep = "")
str.fill1 <- paste("Measures")
str.x1 <- paste("Film types by platform")
lcmodel2 <-
  melt(lca.users$best$probs, level = 2)
lcmodel2$X2 <-
  revalue(lcmodel2$X2,
          c("Pr(2)" = "Probability viewed", "Pr(1)" = "Probability not viewed"))
graphsfindbestlca.function(lcmodel2,
                           title, str.fill1, str.x1)

##Plot model selection variables
plotbicandentropy.function(lca.users$models)

lca.data.users$usertype<-lca.users$best$predclass
lca.data.users$limiteduser <- 0
lca.data.users$limiteduser[lca.data.users$usertype == 2] <- 1
lca.data.users.nononusers <- subset(lca.data.users, usertype != 3)

fit.glmlogit.limiteduser <- glm(
  limiteduser ~
    Education +
    Age +
    HealthImpact +
    NumberChildrenByAge +
    NotVeryConfident +
    Income +
    DEP +
    QLOC +
    AWTV3,
  family = binomial(link = "logit"),
  data = lca.data.users.nononusers,
  control = list(maxit = 50)
)

summary(fit.glmlogit.limiteduser)
PseudoR2(fit.glmlogit.limiteduser, which = "Nagelkerke")
PseudoR2(fit.glmlogit.limiteduser, which = "CoxSnell")
logitgof(fit.glmlogit.limiteduser$y,
         fitted(fit.glmlogit.limiteduser),
         g = 10)

fit.glmlogit.limiteduser0 <- glm(limiteduser ~ 1,
                             family = binomial(link = "logit"), data = lca.data.users.nononusers, control = list(maxit = 50))

exp(cbind(coef(fit.glmlogit.limiteduser),
          confint(fit.glmlogit.limiteduser)))
anova(fit.glmlogit.limiteduser0,
      fit.glmlogit.limiteduser,
      test = "Chisq")

lca.data.users.nononusers$modelglmlimited <- 0
lca.data.users.nononusers$modelglmlimited[fit.glmlogit.limiteduser$linear.predict > 0.5] <-
  1
CrossTable(lca.data.users.nononusers$limiteduser,
           lca.data.users.nononusers$modelglmlimited)
accuracy <-
  table(lca.data.users.nononusers$limiteduser,
        lca.data.users.nononusers$modelglmlimited)
sum(diag(accuracy)) / sum(accuracy)


dev.off()
sink()