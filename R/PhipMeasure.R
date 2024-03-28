PhipMeasure<-function(design,p=15,q=2){
  design=as.matrix(design)
    EucliDist <- function(vec1, vec2) {
      return(sqrt(sum((vec1 - vec2)^2)))
      }
    ############
    phip=c()
    for(i in 1:(nrow(design)-1)){
      phip=c(phip,(EucliDist(design[i,],design[i+1,]))^(-p))
    }
    ##########
    phip=(sum(phip))^(1/p)
    return(phip)
  }