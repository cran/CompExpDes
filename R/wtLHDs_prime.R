wtLHDs_prime<-function(levels,factors,weights=c(0.2,0.4,0.4),iterations=300){
  t0<-Sys.time()
  if(sum(weights)!=1){
    return(message("Please enter weights such that sum of weights = 1."))
  }
  is.prime <- function(value) {
    if (value <= 1) return(FALSE)
    if (value <= 3) return(TRUE)
    if (value %% 2 == 0 || value %% 3 == 0) return(FALSE)
    i <- 5
    while (i * i <= value) {
      if (value %% i == 0 || value %% (i + 2) == 0) return(FALSE)
      i <- i + 6
    }
    return(TRUE)
  }
  
  #####error messsage
  if(levels>factors^2 || levels<factors ||is.prime(factors)==FALSE ||factors<3){
    return(message("Factors, F (>2) is a prime number and Levels, L should be in the range from F to (F^2)"))
  }
  ###############
  s=factors
  row_wise_add=c()
  for(i in 1:(s-1)){
    row_wise_add=c(row_wise_add,s*i+1)
  }
  ###
  column_wise_add=((s^2+1)-row_wise_add)
  #########
  basic=matrix(1:s^2,nrow=s,ncol=s,byrow=TRUE)
  ########
  list_of_first_rows<-list()
  for(k in 1:(s-1)){
    base_row=c(1)
    for(l in 1:(s-1)){
      base_row=c(base_row,base_row[length(base_row)]+row_wise_add[k])
    }
    base_row<-base_row%%(s^2)
    base_row[base_row==0]<-(s^2)
    list_of_first_rows<-append(list_of_first_rows,list(base_row))
  }
  ##############################Develop columns
  list_of_arrays=list()
  for(m in 1:(s-1)){
    array=t(as.matrix(list_of_first_rows[[m]]))
    for(n in 1:(s-1)){
      array=rbind(array,array[nrow(array),]+column_wise_add[m])
    }
    list_of_arrays=append(list_of_arrays,list(array)) 
  }
  ###########################apply mod
  ###########
  list_of_arrays1=lapply(list_of_arrays,function(mat)mat%%s^2)
  replace_element=function(matrix,value){
    matrix[matrix==0]<-value
    return(matrix)
  }
  list_of_arrays2=lapply(list_of_arrays1,function(mat)replace_element(mat,s^2))
  ####Final list of arrays
  list_of_arrays3=append(list(basic),list_of_arrays2)
  ###### Final design
  new=NULL
  a=c(1:s)
  for(i in 0:(s-1)){
    new=rbind(new,a+i)
  }
  new=new%%s
  new[new==0]<-s
  ############  
  Final_list<-list()
  for(i in 1:s){
    Final_list=append(Final_list,list(list_of_arrays3[[i]][,new[,i]]))
  }
  ############# LHD
  LHD=NULL
  for(i in 1:length(Final_list)){
    LHD=rbind(LHD,Final_list[[i]])
  }
  #}
  ####################################Design Generation above
  design<-LHD
  ######
  ############### 
  des<-LHD
  simu<-iterations
  n=levels
  k=factors
  store<-NULL
  store_des<-list()
  for(s in 1:simu){
    des1<-as.matrix(des[(1:n),]) #needed
    des1[des1 > n] <- NA 
    rest<-list()
    for(r in 1:ncol(des1)){
      rest<-append(rest,list(setdiff(1:n,des1[,r])))
    }
    mysample<-function(vec){
      if(length(vec)==1){
        return(vec)
      }else{
        return(sample(vec))
      }
    }
    rest<-lapply(rest,function(vec)mysample(vec))
    for(m in 1:ncol(des1)) {
      des1[is.na(des1[,m]), m] <- rest[[m]]
    }
    ##############
    store_des<-append(store_des,list(des1))
    store<-rbind(store,c(MAC(des1),PhipMeasure(des1,q=2),Maxpro_Measure(des1)))
  }
  weightage<-store%*%diag(weights)
  row_sum<-apply(weightage,1,sum)
  pos<-which.min(row_sum)
  desired_des<-store_des[[pos]]
  ##############
  LHD<-desired_des
  t1=Sys.time()
  time_req=t1-t0
  lm=list("LHD"=LHD,"Maxpro_Measure"=Maxpro_Measure(LHD),"Phi_p Measure"=PhipMeasure(LHD),"Maximum Absolute Correlation"=MAC(LHD),"Total System Time Requires"=time_req) 
  return(lm)
}

