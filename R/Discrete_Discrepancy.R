Discrete_Discrepancy<-function(Design,a,b){
  matrix<-as.matrix(Design)
  s<-ncol(matrix)
  q<-max(matrix)
  n<-nrow(matrix)
  collect_for_sum<-c()
  for(i in 1:(nrow(matrix)-1)){
    for(j in (i+1):nrow(matrix)){
      substract_row<-matrix[i,]-matrix[j,]
      now<-length(substract_row[substract_row!=0])
      hamming<-length(matrix[i,])-now
      collect_for_sum<-c(collect_for_sum,(a/b)^hamming)
    }
  }
  DD<--((a+((q-1)*b))/q)^s+((a^s)/n)+((2*(b^s))/(n^2))*sum(collect_for_sum)
  ############LDD
  psy<-(s*(n-q)/(q*(n-1)))
  gamma<-as.integer(psy)
  LDD<--((a+((q-1)*b))/q)^s+((a^s)/n)+(((n-1)*(b*(1-psy)+a*psy)*b^s)/(n*b))*(a/b)^gamma
  list1=list("Discrete Discrepancy Measure"=DD,"Lower Bound of Discrete Discrepancy Measure"=LDD)
return(list1)
  }
