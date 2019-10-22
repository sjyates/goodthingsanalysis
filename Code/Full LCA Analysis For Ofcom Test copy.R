## sink output to ovearll results file
sink("Full LCA", append = FALSE, split = TRUE)

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

## Set up variables and lists
## Set up bic variables
max_II <- -100000
min_bic <- 0

##Set up function for entropy calculation
entropy <- function (var.p) {
  sum(-var.p * log(var.p))
}

##Set up function to find best result from 10 iterations
##Return best LCA results and BIC values table
lcabestmodel.function <-
  function(analysis.vars,
           analysis.data,
           var.minclasses,
           var.maxclasses) {
    ##Set up temporary return, results and models frames
    lcaresults <- list()
    lcmodels <- data.frame(
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
          nrep = 10,
          na.rm = TRUE,
          verbose = TRUE,
          calc.se = TRUE
        )
      ##Check for best model based on BIC criterion
      if (min_bic == 0) {
        min_bic <- lc$bic
      }
      if (lc$bic < min_bic) {
        min_bic <- lc$bic
        LCA_best_model <- lc
      }
      ##Store criteria for comparative table
      lcmodels[i - 1, 1] <- i
      lcmodels[i - 1, 2] <- lc$llik
      lcmodels[i - 1, 3] <- lc$resid.df
      lcmodels[i - 1, 4] <- lc$bic
      lcmodels[i - 1, 5] <- lc$aic
      lcmodels[i - 1, 6] <- lc$Gsq
      lcmodels[i - 1, 7] <- lc$Chisq
      if (i > var.minclasses) {
        error_post <- mean(apply(lc$posterior, 1, entropy), na.rm = TRUE)
        lcmodels[i - 1 , 8] <-
          round(((error_prior - error_post) / error_prior), 3)
        error_prior <- entropy(lc$P) # class proportions model 2
      } else {
        lcmodels[i, 8] <- c("-")
        error_prior <- entropy(lc$P) # class proportions model 1
      }
    }
    
    ##Create return list
    lcaresults$models <- lcmodels
    lcaresults$best <- LCA_best_model
    writeLines("\nTable of BIC values for model selection\n")
    print(lcaresults$models)
    writeLines("\nBest Model From Analysis\n")
    print(lcaresults$best)
    return(lcaresults)
  }

##Function to plot variable probabilites by class
graphsforlca.function <- function(lcaResults, str.title) {
  # str.filename <-
  #   paste(getwd(), "/plots/", str.title, "1.eps", sep = "")
  # postscript(str.filename)
  lcmodel <- melt(lcaResults, level = 2)
  zp1 <- ggplot(lcmodel, aes(x = L2, y = value, fill = Var2))
  zp1 <- zp1 + geom_bar(stat = "identity", position = "stack")
  zp1 <- zp1 + facet_grid(Var1 ~ .)
  zp1 <-
    zp1 + scale_fill_brewer(type = "seq", palette = "Blues") + theme_bw()
  zp1 <-
    zp1 + ggtitle(str.title)
  zp1 <-
    zp1 + labs(x = "Questionnaire items", y = "Share of item-\nresponse categories", fill =
                 "Weekly (6)\nto < once per year (1)")
  zp1 <- zp1 + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major.y = element_blank()
  )
  zp1 <- zp1 + guides(fill = guide_legend(reverse = TRUE))
  print(zp1)
  # dev.off()
  # str.filename <-
  #   paste(getwd(), "/plots/", str.title, "2.eps", sep = "")
  # postscript(str.filename)
  zp2 <- ggplot(lcmodel, aes(x = Var1, y = value, fill = Var2))
  zp2 <- zp2 + geom_bar(stat = "identity", position = "stack")
  zp2 <- zp2 + facet_wrap(~ L2)
  zp2 <-
    zp2 + scale_x_discrete("Questionnaire items", expand = c(0, 0))
  zp2 <-
    zp2 + scale_y_continuous("Share of item-\nresponse categories", expand = c(0, 0))
  zp2 <- zp2 + scale_fill_brewer(type = "seq", palette = "Blues") +
    theme_bw()
  zp2 <-
    zp2 + ggtitle(str.title)
  zp2 <- zp2 + labs(fill = "Weekly (6)\nto < once per year (1)")
  zp2 <- zp2 + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    #legend.justification=c(1,0),
    #legend.position=c(1,0)
  )
  zp2 <- zp2 + guides(fill = guide_legend(reverse = TRUE))
  print(zp2)
  # dev.off()
}

##Funciton to plot regresssion lines
regressiongraph.function <-
  function(analysis.data,
           str.title,
           str.xaxistitle,
           var.range,
           str.c1,
           str.c2,
           str.c3) {
    # str.filename <-
    #   paste(getwd(), "/plots/", str.title, ".eps", sep = "")
    # postscript(str.filename)
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
    # dev.off()
  }

##Plot BIC and relted values
plotbicandentropy.function <- function(models) {
  results.melt <-
    melt(LCAFreqFullResults$models, id = "NoClasses")
  results.melt
  # str.filename <-
  #   paste(getwd(), "/plots/", "Modelresults.png", sep = "")
  # postscript(str.filename)
  ggplot(results.melt) +
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
  # dev.off()
}
##Function to create regression results table
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

## Function to undertake chisquare analayis and plot graphs
chisquaretest.predictions.function <-
  function(indfactor.data,
           predclass.data) {
    chisquare.results <- chisq.test(indfactor.data, predclass.data)
    residuals.data <- chisquare.results$residuals
    corrplot(residuals.data, is.cor = FALSE)
    contrib.data <-
      100 * residuals.data ^ 2 / chisquare.results$statistic
    round(contrib.data, 3)
    corrplot(contrib.data, is.cor = FALSE)
    return(chisquare.results)
  }

## Set up and analyse 2005 data
X2005.uses.tmp <-
  X2005[, c(
    "I11_1",
    "I11_2",
    "I11_3",
    "I11_4",
    "I11_5",
    "I11_6",
    "I11_7",
    "I11_8",
    "I11_9",
    "I11_10",
    "I11_11",
    "I11_12",
    "I11_13",
    "I11_14",
    "I11_15",
    "I11_16",
    "I11_17",
    "I11_18",
    "I11_19",
    "I11_20",
    "I11_21",
    "I11_22",
    "I11_23",
    "I11_24",
    "I11_25",
    "Class",
    "S2"
  )]
colnames(X2005.uses.tmp) <-
  c(
    "E_mails",
    "On_line_chat",
    "Instant_messaging",
    "Sports_scores",
    "Information_work_school",
    "Infomration_leisure",
    "Information_health",
    "Information_sports",
    "Infomration_entertainements",
    "Bookings",
    "Shopping",
    "Banking_paying_bills",
    "Gaming",
    "Gambling",
    "Buying_selling",
    "Downloading_media",
    "Website_blog",
    "Information_public_services",
    "News",
    "Radio",
    "Adult",
    "Renting_DVDs",
    "Politics",
    "Local_community",
    "Other_uses",
    "NRS_social_class",
    "Age"
  )

lcaVars <- c(
  "E_mails",
  "On_line_chat",
  "Instant_messaging",
  "Sports_scores",
  "Information_work_school",
  "Infomration_leisure",
  "Information_health",
  "Information_sports",
  "Infomration_entertainements",
  "Bookings",
  "Shopping",
  "Banking_paying_bills",
  "Gaming",
  "Gambling",
  "Buying_selling",
  "Downloading_media",
  "Website_blog",
  "Information_public_services",
  "News",
  "Radio",
  "Adult",
  "Renting_DVDs",
  "Politics",
  "Local_community",
  "Other_uses"
)

X2005.uses.data <- X2005.uses.tmp
X2005.uses.data[, lcaVars] <-
  sjmisc::rec(X2005.uses.tmp[, lcaVars], rec = "2=4; 4=2; else=copy")
X2005.uses.data[, lcaVars][X2005.uses.data[, lcaVars] < 0] <- 1
X2005.uses.data[, lcaVars][X2005.uses.data[, lcaVars] > 4] <- 1
X2005.uses.data[, lcaVars][X2005.uses.data[, lcaVars] == "."] <- NA
X2005.uses.data[X2005.uses.data == "."] <- NA
X2005.uses.data.naomit <- na.omit(X2005.uses.data)

## Set up varibale binds for each analysis
lcaColNamesFreq2005Bind <-
  cbind(
    E_mails,
    On_line_chat,
    Instant_messaging,
    Sports_scores,
    Information_work_school,
    Infomration_leisure,
    Information_health,
    Information_sports,
    Infomration_entertainements,
    Bookings,
    Shopping,
    Banking_paying_bills,
    Gaming,
    Gambling,
    Buying_selling,
    Downloading_media,
    Website_blog,
    Information_public_services,
    News,
    Radio,
    Adult,
    Renting_DVDs,
    Politics,
    Local_community,
    Other_uses
  ) ~ 1

##Generate results for basic classes
LCAFreqFullResults2005 <- list()
##Call analysis with min 2 and max 10 classes
LCAFreqFullResults2005 <-
  lcabestmodel.function(lcaColNamesFreq2005Bind, X2005.uses.data.naomit, 2, 10)
##Lattice plot of variable proportions by class
graphsforlca.function(LCAFreqFullResults2005$best$probs,
                      "Internet user groups")
##Store and view model selection as latex and htm file
LCAFreqFullResultsModels2005Latex <-
  xtable(LCAFreqFullResults2005$models)
View(LCAFreqFullResultsModels2005Latex)
str.filename <-
  paste(getwd(), "Modelresults.latex", sep = "")
print(LCAFreqFullResultsModels2005Latex, file = str.filename)
str.filename <-
  paste(getwd(), "Modelresults.htm", sep = "")
print(LCAFreqFullResultsModels2005Latex,
      type = "html",
      file = str.filename)
##Plot model selection variables
plotbicandentropy.function(LCAFreqFullResults2005$models)

LCAFreqFullResults$Classchisq <-
  chisquaretest.predictions.function(X2005.uses.data.naomit$NRS_social_class,
                                     LCAFreqFullResults2005$best$predclass)

LCAFreqFullResults$Agechisq <-
  chisquaretest.predictions.function(X2005.uses.data.naomit$Age,
                                     LCAFreqFullResults2005$best$predclass)


X2007.uses.tmp <-
  X2007[, c(
    "IN13_1",
    "IN13_2",
    "IN13_3",
    "IN13_4",
    "IN13_5",
    "IN13_6",
    "IN13_7",
    "IN13_8",
    "IN13_9",
    "IN13_10",
    "IN13_11",
    "IN14_1",
    "IN14_2",
    "IN14_3",
    "IN14_4",
    "IN14_5",
    "IN14_6",
    "IN14_7",
    "IN15_1",
    "IN15_2",
    "IN15_3",
    "WTV3",
    "QX"
  )]
colnames(X2007.uses.tmp) <-
  c(
    "E_mails",
    "On_line_chat",
    "Instant_messaging",
    "Buying_selling",
    "Gaming",
    "Gambling",
    "Banking_paying_bills",
    "Downloading_software",
    "Website_blog",
    "Radio",
    "Social_media",
    "Information_work_school",
    "Bookings",
    "Infomration_leisure",
    "Information_public_services",
    "News",
    "Politics",
    "Adult",
    "Video_clips_pre_YouTube",
    "Film",
    "Music_pre_Itunes",
    "NRS_social_class",
    "Age"
  )

lcaVars <- c(
  "E_mails",
  "On_line_chat",
  "Instant_messaging",
  "Buying_selling",
  "Gaming",
  "Gambling",
  "Banking_paying_bills",
  "Downloading_software",
  "Website_blog",
  "Radio",
  "Social_media",
  "Information_work_school",
  "Bookings",
  "Infomration_leisure",
  "Information_public_services",
  "News",
  "Politics",
  "Adult",
  "Video_clips_pre_YouTube",
  "Film",
  "Music_pre_Itunes"
)

X2007.uses.data <- X2007.uses.tmp
X2007.uses.data[, lcaVars] <-
  sjmisc::rec(X2007.uses.tmp[, lcaVars], rec = "11=1; 12=4; 13=3; 14=2; else=copy")
X2007.uses.data[, lcaVars][X2007.uses.data[, lcaVars] < 0] <- 1
X2007.uses.data[, lcaVars][X2007.uses.data[, lcaVars] > 4] <- 1
X2007.uses.data[, lcaVars][X2007.uses.data[, lcaVars] == "."] <- NA
X2007.uses.data[X2007.uses.data == "."] <- NA
X2007.uses.data.naomit <- na.omit(X2007.uses.data)

lcaColNamesFreqBind <-
  cbind(
    E_mails,
    On_line_chat,
    Instant_messaging,
    Buying_selling,
    Gaming,
    Gambling,
    Banking_paying_bills,
    Downloading_software,
    Website_blog,
    Radio,
    Social_media,
    Information_work_school,
    Bookings,
    Infomration_leisure,
    Information_public_services,
    News,
    Politics,
    Adult,
    Video_clips_pre_YouTube,
    Film,
    Music_pre_Itunes
  ) ~ 1

##Generate results for basic classes
LCAFreqFullResults2007 <- list()
##Call analysis with min 2 and max 10 classes
LCAFreqFullResults2007 <-
  lcabestmodel.function(lcaColNamesFreqBind, X2007.uses.data.naomit, 2, 10)
##Lattice plot of variable proportions by class
graphsforlca.function(LCAFreqFullResults2007$best$probs,
                      "Internet user groups")
##Store and view model selection as latex and htm file
LCAFreqFullResultsModels2007Latex <-
  xtable(LCAFreqFullResults2007$models)
View(LCAFreqFullResultsModels2007Latex)
str.filename <-
  paste(getwd(), "Modelresults.latex", sep = "")
print(LCAFreqFullResultsModels2007Latex, file = str.filename)
str.filename <-
  paste(getwd(), "Modelresults.htm", sep = "")
print(LCAFreqFullResultsModels2007Latex,
      type = "html",
      file = str.filename)
##Plot model selection variables
plotbicandentropy.function(LCAFreqFullResults2007$models)

LCAFreqFullResults$Classchisq <-
  chisquaretest.predictions.function(X2007.uses.data.naomit$NRS_social_class,
                                     LCAFreqFullResults2007$best$predclass)

LCAFreqFullResults$Agechisq <-
  chisquaretest.predictions.function(X2007.uses.data.naomit$Age,
                                     LCAFreqFullResults2007$best$predclass)


sink()