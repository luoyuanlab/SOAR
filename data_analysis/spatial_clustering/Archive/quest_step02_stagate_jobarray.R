library(mclust)

#### Calculating BIC for each number of clusters ####
args <- commandArgs(trailingOnly=TRUE)

data <- read.csv(paste0(args, "analysis/clustering/obsm_STAGATE.csv"))
spot_ids <- data[,1]
row_size = nrow(data)

print(paste0(">>> ", args, " started calculating BIC <<<"))
bic <- c()
max_cluster = 30
if (nrow(data) < 60) {
  max_cluster = nrow(data)/2
  print(paste0('MAX CLUSTER: ', max_cluster))
} 

for (i in 2:max_cluster) {
  res <- Mclust(data[,-1], i, "EEE")
  if(length(res) == 0) {
    res <- Mclust(data[,-1], i)
  }
  bic <- append(bic, res$BIC[1])
}

#### Visualize using the optimal number of clusters ####
n_clust <- which.min(bic)+1
print(paste0(">>> Optimal number of cluster for", args, ": ", n_clust, " <<<"))
res <- Mclust(data[,-1], n_clust, "EEE")
if(length(res) == 0) {
  res <- Mclust(data[,-1], n_clust)
}
write.table(data.frame(res$classification), file=paste0(args, "analysis/clustering/obsm_STAGATE_cluster_", n_clust, ".csv"), quote=FALSE, sep=",", row.names=spot_ids, col.names=c("cluster"))
