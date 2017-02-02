library(igraph)
library(circlize)
library(ComplexHeatmap)
library(dendextend)
library(ForceAtlas2)
library(cluster)

read.csv("Test cytof data.csv") -> s1
apply(s1,2,mean) -> s1m
apply(s1,2,sd) -> s1s
scale(s1,s1m,s1s) -> s1scaled
as.data.frame(s1scaled) -> s1

clara_clustered <- clara(s1, 50, metric = "manhattan", stand = TRUE,
      samples = 1000, sampsize = (nrow(s1)), rngR = FALSE)

datalist = list()
abundancedatalist = list()

for(i in 1:50){
  dat = data.frame(colMeans(s1[c(clara_clustered$clustering == i),]))
  datalist[[i]] <- dat
  abundancedat = data.frame(clara_clustered$clusinfo[i,1])
  abundancedatalist[[i]] <- abundancedat
}

big_data = do.call(cbind, datalist)
big_data = t(big_data)
clus_num <- c(1:50)
clus_names <- as.character(clus_num)
as.data.frame(big_data, row.names = c(clus_names)) -> cluster_means

abd_data = do.call(cbind, abundancedatalist)
abd_data = t(abd_data)
clus_num <- c(1:50)
clus_names <- as.character(clus_num)
as.data.frame(abd_data, row.names = c(clus_names)) -> cluster_abundance
cluster_abundance[,1] -> cluster_abundance

dist(cluster_means, method = "manhattan") -> distx1
graph.adjacency(as.matrix(distx1),mode="undirected",weighted=TRUE) -> adjgraph

V(adjgraph)$abundance <- cluster_abundance
V(adjgraph)$size <- log10(V(adjgraph)$abundance)*10
cut.off <- mean(E(adjgraph)$weight)+sd(E(adjgraph)$weight)

adjgraph.sp <<- delete_edges(adjgraph,E(adjgraph)[E(adjgraph)$weight < cut.off])

forcedirected<<-layout.forceatlas2(adjgraph.sp, directed = FALSE, iterations = 1000,
                                   linlog = TRUE, pos = NULL, nohubs = FALSE, k = 100, gravity = 0,
                                   ks = 0.2, ksmax = 20, delta = 1, center = NULL,
                                   plotlabels = TRUE )                                                                                                             
V(adjgraph.sp)$expression <- cluster_means[,c("HLA.DR")]
plot(adjgraph.sp, layout = forcedirected, vertex.label.cex = 0.5,vertex.label.color = "black",
     vertex.color = c("navy", "royalblue3", "lightskyblue", "lightsteelblue", "aliceblue","gray95",
                      "mistyrose", "lightpink", "lightcoral", "indianred", "red")
     [1+
       (V(adjgraph.sp)$expression >-0.8)+ 
       (V(adjgraph.sp)$expression >-0.6)+
       (V(adjgraph.sp)$expression >-0.4)+
       (V(adjgraph.sp)$expression >-0.2)+
       (V(adjgraph.sp)$expression >0)+
       (V(adjgraph.sp)$expression >0.2)+
       (V(adjgraph.sp)$expression >0.4)+
       (V(adjgraph.sp)$expression >0.75)+
       (V(adjgraph.sp)$expression >1)+
       (V(adjgraph.sp)$expression >2.5)
     ])
   


