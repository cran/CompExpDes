MAC<-function(matrix){
  matrix<-as.matrix(matrix)
  value<-sort(unique(abs(c(cor(matrix)))),decreasing = TRUE)[2]
  return((value))
}


