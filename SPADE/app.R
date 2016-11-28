library(shiny)
library(igraph)
library(circlize)
library(ComplexHeatmap)
library(dendextend)
setwd("~/Desktop/R work/PUT FILES IN HERE")
read.csv("Sample 1.csv") -> s1
apply(s1,2,mean) -> s1m
apply(s1,2,sd) -> s1s
scale(s1,s1m,s1s) -> s1scaled
as.data.frame(s1scaled) -> s1

SPADE <- function(x,k,mkrs){
  #initial clustering and binning 
  dist(x[,c(mkrs)], method = "manhattan") -> distx
  hclust(distx) -> clus_x
  cutree(clus_x, k = k) ->cut_x
  datalist = list()
  for(i in 1:k){
    dat = data.frame(colMeans(x[c(cut_x == i),mkrs]))
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
  dev.off()
  plot(SPADEgraph)
  }

PHESPADE <- function(x,k,clus,mkrs){
  #initial clustering and binning 
  dist(x, method = "manhattan") -> distx
  hclust(distx) -> clus_x
  cutree(clus_x, k = k) ->cut_x
  
  phedat = data.frame(x[c(cut_x == clus),])
   boxplot.matrix(as.matrix(phedat[,c(mkrs)]), cex = 0.5, pch = 20, las = 2,
                   main = paste("cluster",clus, sep = " "))
  
  }
  

My_Palette <- colorRampPalette(c("navy","aliceblue","bisque","chocolate1","firebrick"))(256)


ui <- shinyUI(navbarPage(title = "SPADE",
                         tabPanel(title = "SPADE",
                                  sidebarLayout(
                                    sidebarPanel(
                                      selectInput("mkrs", "Select which markers to assess", c(colnames(s1)), multiple = TRUE ),
                                      sliderInput(inputId = "kvalue", 
                                                  label = "How many Clusters", 
                                                  value = 50, min = 0, max = 200),
                                      actionButton(inputId = "plotnetwork", 
                                                   label = "Plot!"),
                                      numericInput("clusternumber", "Which cluster phenotype", value = 1),
                                      actionButton(inputId = "plotphe", 
                                                   label = "Plot Phenotype")
                                    ),
                                    mainPanel(
                                      headerPanel("SPADE"),
                                      
                                      p("Simply click on either Tumour or NTB to plot either tissue."),
                                      
                                      p("Use the slide bar to choose a suitable number of clusters."),
                                      p("Images will be exported to the", strong("images"), "folder of the working directory"),
                                      
                                      plotOutput("Network"),
                                      plotOutput("Phenotype")
                                    )
                                  )
    ),
    includeCSS("cyborg-theme2.css")
  )
)


server <- shinyServer(function(input, output) {
   
  observeEvent(input$plotnetwork, {
    output$Network <- renderPlot({
      SPADE(s1, input$kvalue, input$mkrs)
        })
  })
  
  observeEvent(input$plotphe, {
    output$Phenotype <- renderPlot({
      PHESPADE(s1, input$kvalue, input$clusternumber, input$mkrs)
    })
  })
  
})


shinyApp(ui = ui, server = server)

