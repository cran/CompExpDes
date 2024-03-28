UDesigns_1<-function(p,q,type){
  if(p<=2){
    return(message("Please enter p>2 and q>=2, such that v>=6."))
  }
  if(q<2){
    return(message("Please enter p>2 and q>=2, such that v>=6."))
  }
  #########
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
    list1=list(DD,LDD)
    return(list1)
  }
  ################################
  v=p*q
  scheme<-matrix(1:v,p,q,byrow=TRUE)
  if(type==2){
    rectangular_scheme<-list()
    for(i in 1:v){
      row_pos<-which(scheme==i,arr.ind = TRUE)[1,1]
      col_pos<-which(scheme==i,arr.ind = TRUE)[1,2]
      third_AS<-t(t(c(scheme[-row_pos,-col_pos])))
      rectangular_scheme<-append(rectangular_scheme,list(third_AS))
    }
    now_the_scheme<-rectangular_scheme
  }
  ###################
  #Scheme="GD"
  if(type==1){
    GD_scheme<-list()
    for(i in 1:v){
      row_pos<-which(scheme==i,arr.ind = TRUE)[1,1]
      #col_pos<-which(scheme==i,arr.ind = TRUE)[1,2]
      second_AS<-t(t(c(scheme[-row_pos,])))
      GD_scheme<-append(GD_scheme,list(second_AS))
    }
    now_the_scheme<-GD_scheme
  }
  ###########
  list_maker <- function(element, length) {
    result <- list()
    
    for (i in 1:length) {
      column_vector <- rep(i, element)
      result[[i]] <- t(t(column_vector))
    }
    
    return(result)
  }
  list2<-list_maker(length(now_the_scheme[[1]]),length(now_the_scheme))
  for(i in 1:100000){
    design<-cbind(unlist(list2),unlist(now_the_scheme[c(sample(1:v,v,FALSE))]))
    minimum_abs_cor<-(min(cor(design)))
    if(minimum_abs_cor==0){
      colnames(design)<-c( "Regions","Treatments")
      final_list<-list("Uniform_Design"=design,"Number of Factors"=2,"Number of Levels"=v,"Number of Runs"=nrow(design),"Maximum Absolute Correlation"=abs(minimum_abs_cor),"Discrete Discrepancy Measure"=Discrete_Discrepancy(design,1,0.3)[[1]],
                       "Lower Bound of Discrete Discrepancy"=Discrete_Discrepancy(design,1,0.3)[[2]])
      return(final_list)
    }
  }
}

