library(shiny)
library(ComplexHeatmap)
library(dendextend)
library(circlize)
setwd("~/Desktop/R work/PUT FILES IN HERE")
ht <- read.csv("Sample 1.csv")
ht2 <- read.csv("Sample 2.csv")
gt <- read.csv("Sample 1.csv")
zt <- read.csv("Tumourz.csv")
zn <- read.csv("NTBz.csv")

gt <- as.matrix(ht)
apply(gt,2,mean) -> gtm
apply(gt,2,sd) -> gts
scale(gt,gtm,gts) -> gtscaled

as.matrix(ht) -> htm
# log10(htm) -> htl
htm[!rowSums(!is.finite(htm)),] -> htm
apply(htm,2,scale) -> hts

as.matrix(ht2) -> htm2
log10(htm2) -> htl2
htl2[!rowSums(!is.finite(htl2)),] -> htl2
apply(htl2,2,scale) -> hts2

gt <- as.matrix(gt)
apply(gt,2,mean) -> gtm
apply(gt,2,sd) -> gts
scale(gt,gtm,gts) -> gtscaled

ui <- shinyUI(navbarPage(title = "How to:",
  
  tabPanel(title = "Introduction",
  
  headerPanel("Alternative ways to analyse cytometry data"),
  
  wellPanel(
  
  p("This app provides different ways to assess cytometry data"),
  
  p("This app uses flow cytometry data in the form of '.csv' files."), 
    
  p("To create these files, simply export your desired fcs file as a .csv file direct from flowjo. Then, put your files in the R work folder as '.csv' files "),
  
  p( "You can copy and paste the produced graphs directly from the browser!"))
  
  ),
  
  tabPanel(title = "Scaffold Contour",
           
  sidebarLayout(
    sidebarPanel(
  
  p(strong("Choose which sample")),    
      
  actionButton(inputId = "Tum", 
               label = "Sample1"),
  
  actionButton(inputId = "NTB", 
               label = "Sample2"),
    
  sliderInput(inputId = "theta", 
              label = "rotate left/right", 
              value = 90, min = 0, max = 360),
  
  sliderInput(inputId = "phi", 
              label = "rotate up/down", 
              value = 30, min = 0, max = 360),
  
  textInput(inputId = "s1title", label = "Sample 1 title" ,value = "Insert your title for sample 1 here"),
  
  textInput(inputId = "s2title", label = "sample 2 title" ,value = "Insert your title for sample 2 here")
  
  ),
  
mainPanel(
  
  headerPanel("Immune landscapes of Scaffold plots"),
  
  p("Simply click on either Tumour or NTB to plot either landscape."),
  
  p("Use the slide bars to rotate the image in all directions."),
  
  plotOutput("contour")))),

tabPanel(title = "Heatmap",
         
         sidebarLayout(
        sidebarPanel(
          selectInput("markers", "Select which markers to assess", c(colnames(hts)), multiple = TRUE ),
          
          actionButton(inputId = "heat", 
                       label = "Sample 1"),
          
          textInput("heat1tit", "Insert title here",
                    value = "Delete this text and write your title"),
          
          actionButton(inputId = "heat2", 
                       label = "Sample 2"),
          
          textInput("heat2tit", "Insert title here",
                    value = "Delete this text and write your title"),
          
          checkboxInput("ccol", "Cluster by Columns", value = TRUE),
          
          checkboxInput("crow", "Cluster by Rows", value = TRUE)
          
        ),
        mainPanel(
          
          headerPanel("Single cell Heatmaps"),
          
          p("Simply select which markers from the drop down panel you wish to assess and press 'Plot Heatmap' "),
          
          p("It may take a while to process ( up to a few minutes) depending on the size of your file, so only click the button ONCE!"),
          
          plotOutput("heatmap")
          
        )
         )),

tabPanel( title = "Cluster plotting over gates",
          
          
          sidebarLayout(
            sidebarPanel(
              
              selectInput("markersx1", "x-axis", c(colnames(gt)), multiple = FALSE ),
              
              selectInput("markersy1", "y-axis", c(colnames(gt)), multiple = FALSE ),
              
              actionButton(inputId = "plotscatter", 
                           label = "Plot Scatter Graph"),
              
              textInput("scattertitle", "Insert title here",
                        value = "Delete this text and write your title"),
              
              numericInput("kvalue", "How many clusters", value = 1),
              
              actionButton(inputId = "phenotypes", 
                           label = "Show cluster relative phenotypes"),
              
              numericInput("clusternumber", "Which cluster phenotype", value = 1)
              
            ),
            
            
            mainPanel(
              
              plotOutput("scatter"),
              
              plotOutput("phenotypes")
              
              
              
            )
          )
          
),

includeCSS("cyborg-theme2.css")
))



as.matrix(zt) -> zt
log10(zt) -> ztl
yt <-1:ncol(zt)
xt <- 1:nrow(zt)

as.matrix(zn) -> zn
log10(zn) -> znl
yn <-1:ncol(zn)
xn <- 1:nrow(zn)

server <- shinyServer(
  function(input, output) {
  
observeEvent(input$Tum, {
  
  output$contour <- renderPlot({
    
    persp(xt,yt,ztl,
        theta = input$theta, phi = input$phi, scale = FALSE,
        main = input$s1title,
        border = "cyan",shade = 0.2,col = "azure4")
    
    })
  
})  

    observeEvent(input$NTB, {
      
      output$contour <- renderPlot({
        
        persp(xn,yn,znl,
              theta = input$theta, phi = input$phi, scale = FALSE,
              main = input$s2title,
              border = "cyan",shade = 0.2,col = "azure4")
        
      })
      
    })
    
    observeEvent(input$heat, {
      
      output$heatmap <- renderPlot({
        
      Heatmap(hts[,c(input$markers)], col = colorRamp2(c(-4, -2, 0, 2, 4),c("darkblue","cornflowerblue","white","chocolate1","firebrick")), input$heat1tit,
              cluster_columns = input$ccol, cluster_rows = input$crow ,row_dend_width = unit(3,"cm"))
      
      })
      
    })
    
    observeEvent(input$heat2, {
      
      output$heatmap <- renderPlot({
       
        Heatmap(hts2[,c(input$markers)], col = colorRamp2(c(-4, -2, 0, 2, 4),c("darkblue","cornflowerblue","white","chocolate1","firebrick")), input$heat2tit,
                cluster_columns = input$ccol, cluster_rows = input$crow , row_dend_width = unit(3,"cm"))
        
      })
      
    })
    
    observeEvent(input$plotscatter, {
      
      output$scatter <- renderPlot({
        
        kgt <- reactive(kmeans(gt, input$kvalue))
        
        plot(log10(gt[,c(input$markersx1)]),log10(gt[,c(input$markersy1)]),
             xlab = input$markersx,
             ylab = input$markersy,
             pch = 20,
             main = input$scattertitle,
             col = kgt()$cluster)
        
        
      })
      
    }) 
    
    
    observeEvent(input$phenotypes, {
      
      output$phenotypes <- renderPlot({
        
        kgt <- reactive(kmeans(gt, input$kvalue))
        
        barplot((as.matrix(kgt()$centers))[input$clusternumber,])
        
      })
      
    })
    
  }

)

shinyApp(ui = ui, server = server)

