UDesigns_I<-function(levels,type){
  t0=Sys.time()
  factors<-NULL
  for(i in 3:as.integer(levels/2)){
    if(levels%%i==0){
      factors<-c(i,levels/i)
      p<-factors[1]
      q<-factors[2]
      break
    }
  }
  #########
  ################################
  v=p*q
  scheme<-matrix(1:v,p,q,byrow=TRUE)
  if(type=="Good"){
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
  if(type=="Excellent"){
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
  repeat{
    design<-cbind(unlist(list2),unlist(now_the_scheme[c(sample(1:v,v,FALSE))]))
    minimum_abs_cor<-(min(cor(design)))
    if(minimum_abs_cor==0){
      colnames(design)<-NULL
      t1=Sys.time()
      total_time=t1-t0
      final_list<-list("Uniform design"=design,"Number of factors"=2,"Number of levels"=v,"Number of runs"=nrow(design),"Maximum absolute correlation"=abs(minimum_abs_cor),"Discrete discrepancy measure"=Discrete_Discrepancy(design)$`Discrete Discrepancy Measure`,
                       "Lower bound of discrete discrepancy Measure"=Discrete_Discrepancy(design)$`Lower Bound of Discrete Discrepancy Measure`,"Total system time"=total_time)
      return(final_list)
    }
  }
}

##########end UDesigns_I

