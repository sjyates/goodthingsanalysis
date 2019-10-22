# install needed packages
# install.packages(c("ca", "car", "corrplot", "DescTools", "dplyr", "factoextra", "FactoMineR", "forcats", "foreign", "generalhoslem", "ggparallel", "ggplot2", "ggrepel", "gmodels", "haven", "Hmisc", "igraph", "knitr", "lattice", "logistf", "plyr", "poLCA", "RcmdrMisc", "reshape2", "sjlabelled", "tidyr", "xtable"))
# load necessary libraries
library("ca")
library("car")
library("corrplot")
library("DescTools")
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

# Global variables

# Code to load Ofcom SPSS files
location_data_folder <- getwd()
X2005_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2005Report2006/2005_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2005_for_analysis <- X2005_original
X2007_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2007Report2008/2007_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2007_for_analysis_ <- X2007_original
X2009_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2009Report2010/2009_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2009_for_analysis <- X2009_original
X2010_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2010Report2011/2010_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2010_for_analysis <- X2010_original
X2011_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2011Report2012/2011_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2011_for_analysis <- X2011_original
X2012_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2012Report2013/2012_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2012_for_analysis <- X2012_original
X2013_original_ml_only <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2013Report2014/2013_spss_processed_variable_names_ml.sav",
      sep = ""
    )
  )
X2013_for_analysis_ml_only <- X2013_original_ml_only
X2013_original_internet_users_only <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2013Report2014/2013_spss_processed_variable_names_interent_users.sav",
      sep = ""
    )
  )
X2013_for_analysis_internet_users_only <-
  X2013_original_internet_users_only
X2013_original_combined <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2013Report2014/2013_original_data_combined_ml_and_interent_source_ofcom.sav",
      sep = ""
    )
  )
X2013_for_analysis_combined <- X2013_original_combined
X2014_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2014Report2015/2014_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2014_for_analysis <- X2014_original
X2015_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2015Report2016/2015_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2015_for_analysis <- X2015_original
X2016_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2016Report2017/2016_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2016_for_analysis <- X2016_original
X2017_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2017Report2018/2017_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2017_for_analysis <- X2017_original
X2018_original <-
  read_sav(
    paste(
      location_data_folder,
      "/OriginalProcessedDataSPSS/Fieldwork2018Report2019/2018_spss_processed_variable_names.sav",
      sep = ""
    )
  )
X2018_for_analysis <- X2018_original
