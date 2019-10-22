library("Hmisc")
library("ca")
library("reshape2")
library("dplyr")

tableforplot.mca <- function(results.mjca) {
  results.list <- list()
  nd  <- results.mjca$nd
  if (is.na(nd)) {
    nd <- 2
  } else {
    if (nd > length(results.mjca$sv))
      nd <- length(results.mjca$sv)
  }
  # Eigenvalues:
  Dimension  <- 1:length(results.mjca$sv)
  Value      <- round(results.mjca$sv ^ 2, 6)
  Percentage <-
    paste(as.character(round(100 * results.mjca$inertia.e, 2)), "%", sep = "")
  
  tmp <- rbind(Value = as.character(Value),
               Percentage = as.character(Percentage))
  dimnames(tmp)[[2]] <- Dimension
  Eigenvalues <- tmp
  
  tmp <-
    rbind(
      results.mjca$colmass,
      results.mjca$coldist,
      results.mjca$colinertia,
      t(results.mjca$colcoord[, 1:nd])
    )
  tmpnames <- results.mjca$levelnames
  if (!is.na(results.mjca$colsup[1])) {
    tmpnames[results.mjca$colsup] <-
      paste(tmpnames[results.mjca$colsup], "(*)", sep = "")
  }
  dimnames(tmp)[[2]] <- tmpnames
  dn <- paste("Dim.", 1:nd)
  dimnames(tmp)[[1]] <- c("Mass", "ChiDist", "Inertia", dn)
  Column.profiles <- tmp
  cat("\n Eigenvalues:\n")
  print.table(Eigenvalues, width = 4)
  str.filename <-
    paste(getwd(), "McaresultsEigen.csv", sep = "")
  write.csv(Eigenvalues, file = str.filename)
  results.list$eigenvalues <- melt(Eigenvalues)
  cat("\n\n Columns:\n")
  print(round(Column.profiles, 6))
  results.list$cprofile <- t(Column.profiles)
  str.filename <-
    paste(getwd(), "McaresultsDimensions.csv", sep = "")
  write.csv(results.list$cprofile, file = str.filename)
  return(results.list)
}

X2005.uses.data <- X2005[,c("I11_1","I11_2","I11_3","I11_4","I11_5","I11_6","I11_7","I11_8","I11_9","I11_10","I11_11","I11_12","I11_13","I11_14","I11_15","I11_16","I11_17","I11_18","I11_19","I11_20","I11_21","I11_22","I11_23","I11_24","I11_25")]
colnames(X2005.uses.data) <- c("E-mails","On-line chat","Instant Messaging","Sports scores","Information work/school","Infomration leisure","Information health","Information sports","Infomration entertainements","Bookings","Shopping","Banking and paying bills","Gaming","Gambling","Buying and selling","Downloading media","Website or blog","Information public services","News","Radio","Adult","Renting DVDs","Politics","Local community","Other uses")
X2005.uses.data[X2005.uses.data < 0] <- 1
X2005.uses.data[X2005.uses.data > 4] <- 1
X2005.uses.data[X2005.uses.data == "."] <-NA
X2005.uses.data.naomit <- na.omit(X2005.uses.data)
X2005.uses.mca <- mjca(X2005.uses.data.naomit)
X2005.uses.mctab<-tableforplot.mca(X2005.uses.mca)