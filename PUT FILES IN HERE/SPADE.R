#SPADE!!!

read.csv("exported.FCS3.csv") -> x
apply(x,2,mean) -> xm
apply(x,2,sd) -> xs
scale(x,xm,xs) -> xscaled
as.data.frame(xscaled) -> x
My_Palette <- colorRampPalette(c("navy","aliceblue","bisque","chocolate1","firebrick"))(256)
library(igraph)
library(circlize)
library(ComplexHeatmap)
library(dendextend)

SPADE <- function(x,k){
  #initial clustering and binning 
dist(x, method = "manhattan") -> distx
hclust(distx) -> clus_x
cutree(clus_x, k = k) ->cut_x
datalist = list()
for(i in 1:k){
dat = data.frame(colMeans(x[c(cut_x == i),]))
datalist[[i]] <- dat}
#cleaning data and assigning to data frame
big_data = do.call(cbind, datalist)
big_data = t(big_data)
clus_num <- c(1:k)
clus_names <- as.character(clus_num)
as.data.frame(big_data, row.names = c(clus_names)) -> cluster_means
# Heatmapping clusters for easy viewing
mypath4 <- file.path("~/Desktop","Lab R work","PUT FILES IN HERE",
                     "Images",paste("Heatmap_","cluster_", ".png", sep = ""))
png(file = mypath4)
heatmap(as.matrix(cluster_means), Colv = NA, col = My_Palette)
dev.off()
# Saving phenotypes as box and whisker graphs
for(i in 1:k){
  mypath3 <- file.path("~/Desktop","Lab R work","PUT FILES IN HERE",
                       "Images",paste("phenotype_","cluster_", i, ".png", sep = ""))
  png(file = mypath3, width = 1700, units = "px")
  phedat = data.frame(x[c(cut_x == i),])
  boxplot.matrix(as.matrix(phedat), cex = 0.5, pch = 20, las = 2,
                 main = paste("cluster",i, sep = " "))
  dev.off()}
#calculating cluster distances and plotting
dist(cluster_means, method = "manhattan") -> distx1
graph.adjacency(as.matrix(distx1),mode="undirected",weighted=TRUE) -> adjgraph
SPADEgraph <-minimum.spanning.tree(adjgraph) %>%
  set_vertex_attr("color", value = "cornflowerblue")
mypath2 <- file.path("~/Desktop","Lab R work","PUT FILES IN HERE",
                     "Images",paste("Network_", ".png", sep = ""))
png(file = mypath2)
plot(SPADEgraph)
dev.off()}

