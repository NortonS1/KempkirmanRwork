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

#Reading in data and scaling

s1 <- read.csv("Test cytof data.csv")
s1m <- apply(s1, 2, mean)
s1s <- apply(s1, 2, sd)
s1scaled <- scale(s1, s1m, s1s)
s1 <- as.data.frame(s1scaled)

#Clustering function

Cluster <- function(x, k, mkrs, expression, gravity, repel) {
  set.seed(1)
  
  # /////////////////////////////////////////////////////////////////CLARA
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
  cluster_means <-
    as.data.frame(big_data, row.names = c(clus_names))
  
  abd_data <- do.call(cbind, abundancedatalist)
  abd_data <- t(abd_data)
  clus_num <- c(1:k)
  clus_names <- as.character(clus_num)
  cluster_abundance <-
    as.data.frame(abd_data, row.names = c(clus_names))
  cluster_abundance <- cluster_abundance[, 1]
  full_data <- data.frame(cluster_means, cluster_abundance)
  
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
  
  #calculating cluster distances and assigning attributes
  
  distx1 <- dist(cluster_means, method = "manhattan")
  adjgraph <-
    graph.adjacency(as.matrix(distx1), mode = "undirected", weighted = TRUE)
  SPADEgraph <- minimum.spanning.tree(adjgraph)
  
  V(SPADEgraph)$abundance <- full_data[, ncol(full_data)]
  V(SPADEgraph)$size <- (log10(V(SPADEgraph)$abundance)) * 10
  
  V(adjgraph)$abundance <- full_data[, ncol(full_data)]
  V(adjgraph)$size <- log10(V(adjgraph)$abundance) * 10
  
  #setting edge cuttoffs and min vertex sizes
  
  cut.off <- mean(E(adjgraph)$weight) + (sd(E(adjgraph)$weight))
  adjgraph.sp <-
    delete_edges(adjgraph, E(adjgraph)[E(adjgraph)$weight < cut.off])
  adjgraph.sp <<-
    delete_vertices(adjgraph.sp, V(adjgraph.sp)[V(adjgraph.sp)$size < 1])
  
  
  #Returning all Values to plot reactively
  
  forcedirected <<-
    layout.forceatlas2(
      adjgraph.sp,
      directed = FALSE,
      iterations = 1000,
      linlog = TRUE,
      pos = NULL,
      nohubs = FALSE,
      k = 100,
      gravity = 0,
      ks = 50,
      ksmax = 100,
      delta = 1,
      center = NULL,
      plotlabels = TRUE
    )
  
  ceb <<- cluster_edge_betweenness(SPADEgraph)
  cebforce <<- cluster_edge_betweenness(adjgraph.sp)
  SPADEdata <<- SPADEgraph
  all_data <<- full_data
  
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
                 value = 50,
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
      
      
      V(adjgraph.sp)$expression <- all_data[, c(input$expression)]
      plot(
        adjgraph.sp,
        layout = forcedirected,
        vertex.label.cex = 0.5,
        vertex.label.color = "black",
        vertex.color = c(
          "navy",
          "royalblue3",
          "lightskyblue",
          "lightsteelblue",
          "aliceblue",
          "gray95",
          "mistyrose",
          "lightpink",
          "lightcoral",
          "indianred",
          "red"
        )
        [1 +
          (V(adjgraph.sp)$expression > -0.8) +
          (V(adjgraph.sp)$expression > -0.6) +
          (V(adjgraph.sp)$expression > -0.4) +
          (V(adjgraph.sp)$expression > -0.2) +
          (V(adjgraph.sp)$expression > 0) +
          (V(adjgraph.sp)$expression > 0.2) +
          (V(adjgraph.sp)$expression > 0.4) +
          (V(adjgraph.sp)$expression > 0.75) +
          (V(adjgraph.sp)$expression > 1) +
          (V(adjgraph.sp)$expression > 2.5)]
      )
      
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
