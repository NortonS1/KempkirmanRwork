# Learning Network analysis

> dfr1 <- data.frame( ID=1:4,
                                           
                               FirstName=c("John","Jim","Jane","Jill"),
                                         
                              Female=c(F,F,T,T), 
                                         
                               Age=c(22,33,44,55), stringsAsFactors = F )

# formatting edges and verteces
plot(g4, edge.arrow.size=.5, vertex.color="gold", vertex.size=15, 
     
     vertex.frame.color="gray", vertex.label.color="black", 
     
     vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2) 

# Add atributes to vertices/edges

V(g4)$gender <- c("male", "male", "male", "male", "female", "female", "male")

E(g4)$type <- "email" # Edge attribute, assign "email" to all edges

E(g4)$weight <- 10    # Edge weight, setting all existing edges to 10, so can control the weight

plot(g4, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     
     vertex.color=c( "pink", "skyblue")[1+(V(g4)$gender=="male")] ) 


