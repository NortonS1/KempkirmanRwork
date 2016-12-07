#' For a whole file

filefoldchange <- function(file,rows){
  row = 0
  while(row<rows){
    row = row+1
    print(row)
    foldchange(file[row,1],file[row,2],file[row,3],file[row,4])
  }
}