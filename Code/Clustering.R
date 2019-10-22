library("pvclust")
library("mclust")
library("Hmisc")
library("factoextra")
library("NbClust")
library("cluster")
library("dendextend")

ahrc_data <- X1807239_SPSS_With_New_NRS_Grades

cluster.results.list<-list()

lcaDataGenreRecode <-
  ahrc_data[, c(
    "Q2r1",
    "Q2r2",
    "Q2r3",
    "Q2r4",
    "Q2r5",
    "Q2r6",
    "Q2r7",
    "Q2r8",
    "Q2r9",
    "Q2r10",
    "Q2r11",
    "Q2r12",
    "Q2r13",
    "Q2r14",
    "Q2r15",
    "Q2r16",
    "Q2r17"
  )]


colnames(lcaDataGenreRecode) <-
  c(
    "Action_adventure",
    "Animation",
    "Art_house",
    "Comedy",
    "Comic_book",
    "Classic",
    "Dcoumentary",
    "Drama",
    "Family",
    "Fantasy",
    "Foreign_language",
    "Horror",
    "Musicals",
    "Romance",
    "Rom_com",
    "Sci_fi",
    "Thriller"
  )


ahrc_dataNoNAN <- na.omit(ahrc_data)
tlahrc_dataNoNAN  <- t(ahrc_dataNoNAN)

## Standard hierarchical cluster
dist.mat <- dist(tlahrc_dataNoNAN, method = "binary")
fviz_dist(dist.mat)
clust.res <- hclust(dist.mat, method = "complete")
fviz_dend(clust.res, cex = 0.5)
cluster.results.list$standard.dist <- dist.mat
cluster.results.list$standard.result <- clust.res

##Check quality of fit
## Compute cophentic distance
res.coph <- cophenetic(clust.res)
cluster.results.list$standard.cophentic<-res.coph
## Correlation between cophenetic distance and the original distance
cor(dist.mat, res.coph)
cluster.results.list$standard.cophcor <- cor(dist.mat, res.coph)

fviz_dend(
  clust.res,
  k = 9,
  # Cut in 13 groups
  cex = 0.5,
  # label size
  k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
  color_labels_by_k = TRUE,
  # color labels by groups
  rect = TRUE # Add rectangle around groups
)

set.seed(123)
res.pv <- pvclust(ahrc_dataNoNAN, method.dist="binary", 
                  method.hclust="complete", nboot = 10)
cluster.results.list$pv.result <- res.pv
# Default plot
plot(res.pv, cex = 0.5)
pvrect(res.pv, alpha=0.99, pv="au", max.only=FALSE)

par(mar=c(6,6,4,1)+.1)
res.pv %>% as.dendrogram %>% 
  set("branches_k_color", k = 9, value = c("purple", "orange", "blue", "red", "green")) %>%
  plot
res.pv %>% text
res.pv %>% pvrect(alpha=0.99, pv="au", max.only=FALSE)