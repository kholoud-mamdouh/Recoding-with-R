#performing k-means clustering on functional profiles using the k-means algorithm
# Load libraries
library(vegan)

# Load input data
inlist <- readRDS("ko_presabs_embed.RDS")
embedding_mds.df = inlist[[4]]

# Perform k-means clustering with k=4
points = embedding_mds.df[c(1,2)]
kk <- kmeans(points, centers = 4)

# Re-assign cluster IDs based on cluster size
kk$cluster[kk$cluster == names(sort(kk$size, decreasing = TRUE)[1])] <- 1
kk$cluster[kk$cluster == names(sort(kk$size, decreasing = TRUE)[2])] <- 2
kk$cluster[kk$cluster == names(sort(kk$size, decreasing = TRUE)[3])] <- 3
kk$cluster[kk$cluster == names(sort(kk$size, decreasing = TRUE)[4])] <- 4

# Save results
k.clusters <- data.frame(group=names(kk$cluster),clusterID=kk$cluster)
saveRDS(list(clusterings=k.clusters, assignments=kk$cluster, centers=kk$centers), "kmeans_ko_presabs.RDS")
