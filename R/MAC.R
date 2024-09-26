MAC<-function(matrix){
  matrix=as.matrix(matrix)
  mat<-cor(matrix)
  upper_tri_values<-abs(mat[upper.tri(mat)])
  return(max(upper_tri_values))
}

