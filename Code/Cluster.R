library("pvclust")
library("mclust")
library("Hmisc")
library("factoextra")
library("NbClust")
library("cluster")
library("dendextend")

#cluster.results.list() <- list()

ClusterData <-
  SPSS[, c(
    'IN8AA',
    'IN8AB',
    'IN8AC',
    'IN8AD',
    'IN8AE',
    'IN8AF',
    'IN8AG',
    'IN8AH'
  )]
colnames(ClusterData) <-
  c(
    "Equipment too expensive",
    "Connection too expensive",
    "Getting on line too comlicated",
    "Using the internet is too complicated",
    "I dont have the right equipment",
    "I dont have the right help",
    "It is just not for me",
    "I dont trust the internet"
  )

ClusterDataNoNAN <- na.omit(ClusterData)
tClusterDataNoNAN  <- t(ClusterDataNoNAN)

## Standard hierarchical cluster
dist.mat <- dist(tClusterDataNoNAN, method = "binary")
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
  k = 8,
  # Cut in 13 groups
  cex = 0.5,
  # label size
  k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
  color_labels_by_k = TRUE,
  # color labels by groups
  rect = TRUE # Add rectangle around groups
)

set.seed(123)
res.pv <- pvclust(ClusterDataNoNAN, method.dist="binary", 
                  method.hclust="complete", nboot = 10)
# cluster.results.list$pv.result <- res.pv
# Default plot
plot(res.pv, cex = 0.5)
pvrect(res.pv, alpha=0.95, pv="au", max.only=FALSE)

par(mar=c(6,6,4,1)+.1)
res.pv %>% as.dendrogram %>% 
  set("branches_k_color", k = 8, value = c("purple", "orange", "blue", "red", "green")) %>%
  plot
res.pv %>% text
res.pv %>% pvrect(alpha=0.95, pv="au", max.only=FALSE)