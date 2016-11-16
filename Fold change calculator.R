# ddCT fold change calculator
foldchange <- function(cnCT,tnCT,ctCT,ttCT){
  difc <- cnCT-ctCT
  dift <- tnCT-ttCT
  ddCT <- difc-dift
  print(paste("ddCT = ",ddCT))
  Fold <- 2^-ddCT
  print(paste("Fold-change =",Fold))
  
}

filefoldchange <- function(file,rows){
  row = 0
  while(row<rows){
    row = row+1
    print(row)
    foldchange(file[row,1],file[row,2],file[row,3],file[row,4])
  }
}