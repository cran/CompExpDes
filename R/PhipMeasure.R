PhipMeasure<-function (design, p = 15, q = 1){
  design = as.matrix(design)
  GeneralDist <- function(vec1, vec2, q) {
    if (q == 1) {
      return(sum(abs(vec1 - vec2)))  # Manhattan Distance for q = 1
    } else {
      return(sum((abs(vec1 - vec2))^q)^(1/q))  # Minkowski Distance for q > 1
    }
  }
  phip = c()
  for (i in 1:(nrow(design) - 1)) {
    for (j in (i + 1):nrow(design)) {  
      phip = c(phip, (GeneralDist(design[i, ], design[j, ], q))^(-p))
    }
  }
  phip = (sum(phip))^(1/p)
  return(phip)
}