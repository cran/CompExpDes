LHDs_I<-function(levels,factors,weight=c(0.3,0.3,0.4),iterations=400){
  if(levels>factors^2 || levels<=factors){
    return(message("Levels should be in the range from (factos+2) to (factors^2)"))
  }
  t0<-Sys.time()
  ############
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
  if(is.prime(factors)==FALSE){
    return(message("Please enter number of factors interms of prime numbers"))
  }
  #is.prime(factors)
  ###############
  if(is.prime(factors)==FALSE || factors<3){
    return(message("Please enter a prime number (>2)"))
  }else{
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
  }
  ####################################Design Generation above
  design<-LHD
  ######
  if(levels==factors^2){
    t1=Sys.time()
    time_req=t1-t0
    lm=list("LHD"=LHD,"Maxpro_Measure"=Maxpro_Measure(LHD),"Phi_p Measure"=PhipMeasure(LHD),"Maximum Absolute Correlation"=MAC(LHD),"Total System Time Requires"=time_req) 
    return(lm)
  }
  
 ###############@@@@@@@@@@@@@@@@@@@@@@@@@@@
  final<-LHD
  mod_final<-final[1:levels,]
  mod_final<-apply(mod_final,2,rank)
  ####################
  p1<-mod_final[1:factors,]
  p2<-mod_final[(1+factors):(nrow(mod_final)),]
  partitions_list<-list(p1,p2)
  store_mp_phip_mac<-NULL
  store_des<-list()
  for(i in 1:iterations){
    des<-lapply(partitions_list,function(mat)apply(mat,2,sample))
    des<-do.call(rbind,des)
    #des<-rbind(p1,apply(p2,2,sample))
    store_des<-append(store_des,list(des))
    store_mp_phip_mac<-rbind(store_mp_phip_mac,c(Maxpro_Measure(des),PhipMeasure(des),MAC(des)))
    if(store_mp_phip_mac[i,1]<0.1 && store_mp_phip_mac[i,2]<0.1 && store_mp_phip_mac[i,3]<0.2){
      final_des<-des
      t1=Sys.time()
      time_req<-t1-t0
      lm=list("LHD"=final_des,"Maxpro_Measure"=store_mp_phip_mac[i,1],"Phi_p_measure"=store_mp_phip_mac[i,2],"MAC"=store_mp_phip_mac[i,3],"Total System Time Requires"=time_req)
      return(lm)
    }
    if(i==iterations){
      measure_mat<-store_mp_phip_mac
      measure_mat[,1]<-store_mp_phip_mac[,1]*weight[1]
      measure_mat[,2]<-store_mp_phip_mac[,2]*weight[2]
      measure_mat[,3]<-store_mp_phip_mac[,3]*weight[3]
      wtsum<-apply(measure_mat,1,sum)
      min_value<-which(wtsum==min(wtsum))
      final_des<-store_des[[min_value]]
      t1=Sys.time()
      time_req<-t1-t0
      lm=list("LHD"=final_des,"Maxpro_Measure"=store_mp_phip_mac[min_value,1],"Phi_p_measure"=store_mp_phip_mac[min_value,2],"MAC"=store_mp_phip_mac[min_value,3],"Total System Time Requires"=time_req)
      return(lm)
    }
  }
}

