
Maxpro_Measure<-function(Design){
 matrix<-as.matrix(Design)
  sumsq<-NULL
  for(i in 1:(nrow(matrix)-1)){
    j=i+1
    while(j<=nrow(matrix)){
      sumsq=c(sumsq,(1/prod(((matrix[i,]-matrix[j,])^2))))
      j=j+1
    }
  }
  return((sum(sumsq)/choose(nrow(matrix),2))^(1/ncol(matrix)))
}


