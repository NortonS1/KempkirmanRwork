#SPADE!!!

read.csv("Sample 1.csv") -> x
library(igraph)
SPADE <- function(x,k){

dist(x, method = "manhattan") -> distx
hclust(distx) -> clus_x
cutree(clus_x, k = k) ->cut_x


datalist = list()
for(i in 1:k){
dat = data.frame(colMeans(x[c(cut_x == i),]))
datalist[[i]] <- dat
}
big_data = do.call(cbind, datalist)
big_data = t(big_data)
dist(big_data, method = "manhattan") -> distx1
#final graphing
graph.adjacency(as.matrix(distx1),mode="directed",weighted=TRUE) -> adjgraph
minimum.spanning.tree(adjgraph) -> SPADEgraph
plot(SPADEgraph)

return(big_data)
}

