max_coincidence_number<-function(matrix){
  matrix=as.matrix(matrix)
  max_coincidence<-NULL
  for(i in 1:((nrow(matrix)-1))){
    for(j in setdiff((1:nrow(matrix)),i)){
      sub<-matrix[i,]-matrix[j,]
      max_coincidence<-c(max_coincidence,length(sub[sub==0]))
    }
  }
  return(max(max_coincidence))
}
