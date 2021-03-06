#SPADE!!!

read.csv("Sample 1-Sam.csv") -> x
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
abundancedatalist = list()
for(i in 1:k){
dat = data.frame(colMeans(x[c(cut_x == i),]))
datalist[[i]] <- dat
abundancedat = data.frame(dim(x[c(cut_x == i),]))
abundancedatalist[[i]] <- abundancedat}
#cleaning data and assigning to data frame
big_data = do.call(cbind, datalist)
big_data = t(big_data)
clus_num <- c(1:k)
clus_names <- as.character(clus_num)
as.data.frame(big_data, row.names = c(clus_names)) -> cluster_means

abd_data = do.call(cbind, abundancedatalist)
abd_data = t(abd_data)
clus_num <- c(1:k)
clus_names <- as.character(clus_num)
as.data.frame(abd_data, row.names = c(clus_names)) -> cluster_abundance
cluster_abundance[,1] -> cluster_abundance

full_data = data.frame(cluster_means, cluster_abundance)
# View(full_data[,-ncol(full_data)])
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

dist(full_data[,-ncol(full_data)], method = "manhattan") -> distx1
adjgraph <- graph.adjacency(as.matrix(distx1),mode="upper",weighted=TRUE)
SPADEgraph <<-minimum.spanning.tree(adjgraph)
  V(adjgraph)$abundance <- full_data[,ncol(full_data)]
V(adjgraph)$size <- log10(V(adjgraph)$abundance)*10
# E(adjgraph)$edge.width <- (E(adjgraph)$weight)
mypath2 <- file.path("~/Desktop","Lab R work","PUT FILES IN HERE",
                     "Images",paste("Network_", ".png", sep = ""))
cut.off <- mean(E(adjgraph)$weight)+sd(E(adjgraph)$weight)
adjgraph.sp <- delete_edges(adjgraph,E(adjgraph)[E(adjgraph)$weight < cut.off])
print(adjgraph.sp)
print(E(adjgraph)$weight)
png(file = mypath2)
layout.forceatlas2(adjgraph.sp, iterations = 100, linlog = TRUE, k = 100, gravity = 1, ks = 100 ) -> forcedirected
plot(adjgraph.sp, layout = forcedirected)
dev.off()}

