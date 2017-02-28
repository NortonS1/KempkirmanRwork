library(shiny)
library(igraph)
library(circlize)
library(ComplexHeatmap)
library(dendextend)
library(ForceAtlas2)
library(cluster)

outputDir  = "Images"

# assigning palette
My_Palette <-
  colorRampPalette(c("navy", "aliceblue", "bisque", "chocolate1", "firebrick"))(256)

colourBreaks <- 100
# palette <- colorRampPalette(c('blue','red'))

#Reading in data and scaling

s1 <- read.csv("Tumour.csv")
s1 <- as.data.frame(s1)


# gate names
gatingFiles <- lapply(Sys.glob("gating_*.csv"), read.csv)
fileNames<- as.list(dir(pattern = "gating_*"))
gateNames <- list()
for(i in 1:length(fileNames)){
file <- fileNames[[i]]
gateID  <- gsub("gating_","",file)
gateID <- gsub(".csv","",gateID)
gateNames[[i]] <- gateID
}

gateAbundance <- list()

for (i in 1:length(fileNames)){
  rowNum <- nrow(gatingFiles[[i]])
  gateAbundance [[i]] <- rowNum
}
gate_abundance <- do.call(cbind, gateAbundance)
gate_abundance <- t(gate_abundance)
#scaling data together then re-separating

numGates <- length(gatingFiles)
gatelist <- list()
for(i in 1:numGates){
  g <- gatingFiles[[i]]
  gColmean <- colMeans(g)
  gatelist[[i]] <- gColmean
}
gating_data <- do.call(cbind, gatelist)
gating_data <- t(gating_data)
LM_names <- as.character(gateNames)
gating_data <- as.data.frame(gating_data, row.names = c(LM_names))
# pre_scale_all <- rbind(s1, gating_data)

s1m <- apply(s1, 2, mean)
s1s <- apply(s1, 2, sd)
s1scaled <- scale(s1, s1m, s1s)
scaled_all <- as.data.frame(s1scaled)
s1 <- as.data.frame(s1scaled)

g1 <-gating_data
g1m <- apply(g1, 2, mean)
g1s <- apply(g1, 2, sd)
g1scaled <- scale(g1, g1m, g1s)
gating_data <- as.data.frame(g1scaled)

# gating_data <- do.call()

#Clustering function

Cluster <- function(x, k, mkrs, expression, gravity, repel) {
  set.seed(1)
  
  # /////////////////////////////////////////////////////////////////CLARA
  # Landmark clustering with gated files
  
  gating_data <- gating_data[,c(mkrs)]

  
  #raw data clustering
  clara_clustered <-
    clara(
      x[, c(mkrs)],
      k,
      metric = "manhattan",
      stand = TRUE,
      samples = 1000,
      sampsize = (nrow(x)),
      rngR = FALSE
    )
  datalist <- list()
  abundancedatalist <- list()
  
  for (i in 1:k) {
    dat <-
      data.frame(colMeans(x[c(clara_clustered$clustering == i), mkrs]))
    datalist[[i]] <- dat
    abundancedat <- data.frame(clara_clustered$clusinfo[i, ])
    abundancedatalist[[i]] <- abundancedat
  }
  
  # ////////////////////////////////////////////////////////////////CLARA
  
  #cleaning data and assigning to data frame
  big_data <- do.call(cbind, datalist)
  big_data <- t(big_data)
  clus_num <- c(1:k)
  clus_names <- as.character(clus_num)
  cluster_means <-as.data.frame(big_data, row.names = c(clus_names))
  
  c1m <- apply(cluster_means, 2, mean)
  c1s <- apply(cluster_means, 2, sd)
  c1scaled <- scale(cluster_means, c1m, c1s)
  cluster_scaled <- as.data.frame(c1scaled)
  cluster_means <- as.data.frame(cluster_scaled)
 
  
  abd_data <- do.call(cbind, abundancedatalist)
  abd_data <- t(abd_data)
  clus_num <- c(1:k)
  clus_names <- as.character(clus_num)
  cluster_abundance <-as.data.frame(abd_data, row.names = c(clus_names))
  cluster_abundance <- cluster_abundance[, 1]
  full_data <- data.frame(cluster_means, cluster_abundance)
  
  cluster_abundance <- gate_abundance
  full_gate <- cbind(gating_data, cluster_abundance)
  full_data <- rbind(full_data, full_gate)
  
 print("full data complete")
 print("----------------------------------------------")
  # # Switch to the image output directory
  setwd(outputDir)
  
  # Heatmapping clusters for easy viewing
  png(file = paste("Heatmap_", "cluster_", i, ".png"))
  heatmap(as.matrix(cluster_means), Colv = NA, col = My_Palette)
  dev.off()
  
  # Saving phenotypes as box and whisker graphs
  for (i in 1:k) {
    png(
      file = paste("phenotype_", "cluster_", i, ".png"),
      width = 1700,
      units = "px"
    )
    phedat <- data.frame(x[c(clara_clustered$clustering == i), ])
    boxplot.matrix(
      as.matrix(phedat),
      cex = 0.5,
      pch = 20,
      las = 2,
      main = paste0("cluster", i)
    )
    dev.off()
  }
  
  print("box and whisker complete")
  print("----------------------------------------------")
  #calculating cluster distances and assigning attributes
  cluster_means_gated <- rbind(cluster_means, gating_data)
  distxg <- dist(cluster_means_gated, method = "manhattan")
  adjgraph_gated <- graph.adjacency(as.matrix(distxg), mode = "undirected", weighted = TRUE)
  V(adjgraph_gated)$abundance <- full_data[, ncol(full_data)]
  V(adjgraph_gated)$size <- log10(V(adjgraph_gated)$abundance) * 10
  
  print("adjgraph_gated complete")
  print("----------------------------------------------")
  distx1 <- dist(cluster_means, method = "manhattan")
  SPADE_adjgraph <- graph.adjacency(as.matrix(distx1), mode = "undirected", weighted = TRUE)
  SPADEgraph <- minimum.spanning.tree(SPADE_adjgraph)
  
  print("MST made")
  print("----------------------------------------------")
  V(SPADEgraph)$abundance <- full_data[(1:(nrow(full_data)-numGates)), ncol(full_data)]
  V(SPADEgraph)$size <- (log10(V(SPADEgraph)$abundance)) * 10
  
  print("MST V abundance set")
  print("----------------------------------------------")
  #setting edge cuttoffs and min vertex sizes
  
  print("removing none LM edges....")
  print("----------------------------------------------")
  
  i <- 1
  E(adjgraph_gated)$LM <- FALSE
  del_edge_list <- list()
  edgelist <<- get.edgelist(adjgraph_gated)
  for (i in (1:nrow(get.edgelist(adjgraph_gated)))){
      if(nchar(edgelist[i,1]) < 3 && nchar(edgelist[i,2]) < 3){
       
        
        E(adjgraph_gated)$LM[get.edge.ids(adjgraph_gated, vp = c(edgelist[i,1],edgelist[i,2]))] <- FALSE
        
      }else{
        
        E(adjgraph_gated)$LM[get.edge.ids(adjgraph_gated, vp = c(edgelist[i,1],edgelist[i,2]))] <- TRUE
        
      }
      print (i)
  i = i+1
  }
  
  print("set true/false for LM")
  adjgraph_gated <- delete_edges(adjgraph_gated,(which(E(adjgraph_gated)$LM == FALSE)))
  adjgraph_gated <-delete_vertices(adjgraph_gated, V(adjgraph_gated)[V(adjgraph_gated)$abundance < 4])
  
  adjgraph_gated_2 <- delete_edge_attr(adjgraph_gated, "LM")
  print("Deleted LM edge_attr")
  
  E(adjgraph_gated_2)$weight<- scale(E(adjgraph_gated_2)$weight)
  
  print(adjgraph_gated_2)
  
  list <- LM_names
  for (name in list){
    LM_edges <- E(adjgraph_gated_2)[from(name)] 
    LM_edge_weight <- E(adjgraph_gated_2)$weight[LM_edges]
    cutoff <- mean(LM_edge_weight)
    adjgraph_gated_2  <- delete_edges(adjgraph_gated_2, LM_edges[which(LM_edge_weight < cutoff)])
    print(adjgraph_gated_2)
  }
  adjgraph_gated_3 <- adjgraph_gated_2 
  print("removed weak edges and small verteces")
  print("----------------------------------------------")
  #Returning all Values to plot reactively
  
  forcedirected <<-
    layout.forceatlas2(
      adjgraph_gated_3,
      directed = FALSE,
      iterations = 100,
      linlog = TRUE,
      pos = NULL,
      nohubs = TRUE,
      k = 10,
      gravity = 0,
      ks = 5,
      ksmax = 50,
      delta = 2,
      center = NULL,
      plotlabels = TRUE
    )
  
  print("calculated forceatlas co-ords")
  print("----------------------------------------------")
  ceb <<- cluster_edge_betweenness(SPADEgraph)
  SPADEdata <<- SPADEgraph
  all_data <<- full_data[which(full_data[, ncol(full_data)] >= 4),] 
  View(all_data)
  adjgraph_gated_3 <<- adjgraph_gated_3
  setwd("..")
  return(length(ceb))
}

#Cluster data for phenotype plotting

PHESPADE <- function(x, k, clus, mkrs2) {
  #initial clustering and binning
  set.seed(1)
  distx <- dist(x, method = "manhattan")
  clus_x <- hclust(distx)
  cut_x <- cutree(clus_x, k = k)
  # phenotype outputs based on user choice
  phedat <- data.frame(x[c(cut_x == clus), ])
  boxplot.matrix(
    as.matrix(phedat[, c(mkrs2)]),
    cex = 0.5,
    pch = 20,
    las = 2,
    main = paste("cluster", clus, sep = " ")
  )
}

#Begin UI//////////////////////////////////////////////////////////////////////////////////////////////////////

ui <- shinyUI(navbarPage(
  title = "Cluster App",
  tabPanel(title = "Cluster App",
           sidebarLayout(
             sidebarPanel(
               selectInput(
                 "mkrs",
                 "Select which markers to cluster",
                 c(colnames(s1)),
                 multiple = TRUE
               ),
               sliderInput(
                 inputId = "kvalue",
                 label = "How many Clusters",
                 value = 100,
                 min = 0,
                 max = 200
               ),
               actionButton(inputId = "docluster",
                            label = "Cluster"),
               p(""),
               actionButton(inputId = "plotnetwork",
                            label = "Plot MST!"),
               p(""),
               actionButton(inputId = "plotnetworkcoloured",
                            label = "Plot by marker"),
               p(""),
               actionButton(inputId = "plotnetworkforce",
                            label = "Plot force directed"),
               numericInput("gravity", "gravity", value = 0),
               numericInput("repel", "repel", value = 100),
               numericInput("ks", "ks", value = 50),
               numericInput("ksmax", "ksmax", value = 100),
               checkboxInput("linlog", "linlog", value = TRUE),
               checkboxInput("hubs", "hubs", value = FALSE),
               
               p(""),
               selectInput(
                 "expression",
                 "Select which marker expression to colour the plot by",
                 c(colnames(s1)),
                 multiple = FALSE
               ),
               numericInput("clusternumber", "Which cluster phenotype", value = 1),
               selectInput(
                 "mkrs2",
                 "Select which markers to assess",
                 c(colnames(s1)),
                 multiple = TRUE
               ),
               actionButton(inputId = "plotphe",
                            label = "Plot Phenotype")
             ),
             mainPanel(
               headerPanel("Cluster App"),
               p(
                 strong(
                   "Simply choose the markers you wish to use for clustering from the drop down menu."
                 )
               ),
               p(
                 strong("Then use the slide bar to choose a suitable number of clusters.")
               ),
               p(
                 strong("Once you're happy with the settings press cluster and you're away!")
               ),
               p("Depending on your data size, this may take a minute..."),
               p(
                 strong(
                   "When it pops up saying clustering is complete, simply click the graphing buttons and your plot will appear"
                 )
               ),
               p(
                 strong(
                   "Images will be exported to the",
                   strong("images"),
                   "folder of the working directory"
                 )
               ),
               textOutput("working"),
               textOutput("done"),
               textOutput("clusterdone"),
               plotOutput("Network"),
               p(""),
               plotOutput("Networkcoloured"),
               p(""),
               plotOutput("Networkforce"),
               p(""),
               plotOutput("Phenotype")
             )
           )),
  
  #CSS Themeing
  includeCSS("cyborg-theme2.css")
))

#END UI ///////////////////////////////////////////////////////////////////////////////////////////////

#BEGIN SERVER /////////////////////////////////////////////////////////////////////////////////////////

server <- shinyServer(function(input, output) {
  #Do Clustering
  
  observeEvent(input$docluster, {
    print("Clustering...")
    Cluster(s1,
            input$kvalue,
            input$mkrs,
            expression = NULL,
            input$gravity,
            input$repel)
    output$clusterdone <- renderText({
      print("Clustering Complete")
    })
  })
  
  #Plot SPADE with ceb colouring
  
  observeEvent(input$plotnetwork, {
    output$working <- renderText({
      print("Rendering plot...")
    })
    output$Network <- renderPlot({
      set.seed(1)
      plot(ceb, SPADEdata)
    })
    output$done <- renderText({
      print("Complete!")
    })
  })
  
  #plot SPADE with user input expression colouring
  
  observeEvent(input$plotnetworkcoloured, {
    output$working <- renderText({
      print("Rendering plot...")
    })
    output$Networkcoloured <- renderPlot({
      set.seed(1)
      V(SPADEdata)$expression <- all_data[, c(input$expression)]
      plot(
        SPADEdata,
        vertex.label.cex = 0.5,
        vertex.label.color = "black",
        vertex.color = c(
          "royalblue3",
          "dodgerblue1",
          "darkorchid1",
          "chocolate1",
          "brown1"
        )[1 + (V(SPADEdata)$expression > (mean(V(SPADEdata)$expression)) - (sd(V(SPADEdata)$expression))) +
            (V(SPADEdata)$expression > mean(V(SPADEdata)$expression)) +
            (V(SPADEdata)$expression > (mean(V(SPADEdata)$expression)) + (sd(V(SPADEdata)$expression)))]
      )
    })
    output$done <- renderText({
      print("Complete!")
    })
  })
  
  #plot forceAtlas clustering with user input coulouring
  
  observeEvent(input$plotnetworkforce, {
    output$working <- renderText({
      print("Rendering plot...")
    })
    output$Networkforce <- renderPlot({
      set.seed(1)
      
      
      V(adjgraph_gated_3)$expression <<- all_data[, c(input$expression)]
      print(V(adjgraph_gated_3)$expression)
      V(adjgraph_gated_3)$color <- colorRampPalette(c("bisque","firebrick"))(100)[as.numeric (c((cut(V(adjgraph_gated_3)$expression, breaks= colourBreaks))))]
      
      
      plot(
        adjgraph_gated_3,
        layout = forcedirected,
        vertex.label.cex = 0.8,
        vertex.label.color = "black")
      
    })
    output$done <- renderText({
      print("Complete!")
    })
  })
  
  #plot phenotype of given cluster
  
  observeEvent(input$plotphe, {
    output$Phenotype <- renderPlot({
      PHESPADE(s1, input$kvalue, input$clusternumber, input$mkrs2)
    })
  })
})
# END SERVER  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

shinyApp(ui = ui, server = server)
