#SPADE!!!

read.csv("Sample 1.csv") -> x
library(igraph)

SPADE <- function(x,k){
  
  #initial clustering and binning 
dist(x, method = "manhattan") -> distx
hclust(distx) -> clus_x
cutree(clus_x, k = k) ->cut_x
datalist = list()
for(i in 1:k){
dat = data.frame(colMeans(x[c(cut_x == i),]))
datalist[[i]] <- dat
}
#cleaning data and assigning to data frame
big_data = do.call(cbind, datalist)
big_data = t(big_data)
clus_num <- c(1:k)
clus_names <- as.character(clus_num)
as.data.frame(big_data, row.names = c(clus_names)) -> cluster_means

#calculating cluster distances and plotting
dist(cluster_means, method = "manhattan") -> distx1
graph.adjacency(as.matrix(distx1),mode="undirected",weighted=TRUE) -> adjgraph
minimum.spanning.tree(adjgraph) -> SPADEgraph
tkplot(SPADEgraph)
}

