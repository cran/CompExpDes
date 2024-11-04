
LHDs_II<-function(levels,factors,weight=c(0.3,0.3,0.4),iterations=400){
  Levels=levels
  factor=factors
  iteration=iterations
  n=2*factor+1
  v=choose(n,2)
  if(levels>v || levels<=factors){
    return(message("Levels, L should be in the range from (F+2) to sC2, where s = 2F+1"))
  }
  
  
  seq1=c(1:(n-1))
  seq2=c((n-1):2)
  list<-list()
  ######################
  t0=Sys.time()
  ##############
  j=1
  for(j in seq1){
    a=c(j)
    b<-c(j)
    i=1
    while(i <= length(seq2)){
      b=c(b,a+seq2[i])
      a<-b[length(b)]
      i=i+1
    }
    list<-append(list,list(b))
    seq2=seq2[-length(seq2)]
  }
  ###
  Tri1<-list
  Tri2<-Tri1[length(Tri1):1]
  #############
  array1<-NULL
  for(i in 1:length(Tri1)){
    array1<-rbind(array1,c(Tri1[[i]],Tri2[[i]]))
  }
  ################
  final_array=t(array1)
  ######################
  partition1<-final_array[,1:((n-1)/2)]
  partition2<-final_array[,(1+((n-1)/2)):(n-1)]
  #################################################
  ###########################################
  lhd1<-partition1
  lhd2<-partition2
  lhd3<-cbind(lhd1,lhd2)
  tr1=lhd1
  tr2=lhd2
  if(((n-1)/2)%%2!=0){
    middle<-median(1:((n-1)/2))
    tr22<-tr2
    tr2[,middle]<-tr22[,(middle+1)]
    tr2[,middle+1]<-tr22[,middle]
    rank1<-rank(apply(tr1,2,sum))
    rank2<-rank(apply(tr2,2,sum))
  }
  #write.table(lhd2,"clipboard",sep="\t")
  tr=rbind(tr1,tr2)
  ################
  rest_mat<-NULL
  for(i in 1:ncol(tr)){
    setdif<-setdiff((1:choose(n,2)),tr[,i])
    rest_mat<-cbind(rest_mat,setdif)
  }
  ############
  ####################
  final<-rbind(tr,rest_mat)
  mod_final<-final[1:levels,]
  mod_final<-apply(mod_final,2,rank)
  ####################
  p1<-mod_final[1:factors,]
  p2<-mod_final[(1+factors):levels,]
  partitions_list<-list(p1,p2)
  store_mp_phip_mac<-NULL
  store_des<-list()
  for(i in 1:iterations){
    des<-lapply(partitions_list,function(mat)apply(mat,2,sample))
    des<-do.call(rbind,des)
    store_des<-append(store_des,list(des))
    store_mp_phip_mac<-rbind(store_mp_phip_mac,c(Maxpro_Measure(des),PhipMeasure(des),MAC(des)))
    if(store_mp_phip_mac[i,1]<0.1 && store_mp_phip_mac[i,2]<0.1 && store_mp_phip_mac[i,3]<0.2){
      final_des<-des
      t1=Sys.time()
      time_req=t1-t0
      colnames(final_des)<-NULL
      lm=list("Latin hypercube design"=final_des,"Number of factors"=ncol(final_des),"Number of levels"=nrow(final_des),"Maxpro measure"=store_mp_phip_mac[i,1],"Phi_p measure"=store_mp_phip_mac[i,2],"Maximum absolute correlation"=store_mp_phip_mac[i,3],"Total system time"=time_req)
      return(lm)
    }
    if(i==iterations){
      measure_mat<-store_mp_phip_mac
      measure_mat[,1]<-store_mp_phip_mac[,1]*weight[1]
      measure_mat[,2]<-store_mp_phip_mac[,2]*weight[2]
      measure_mat[,3]<-store_mp_phip_mac[,3]*weight[3]
      wtsum<-apply(measure_mat,1,sum)
      min_value<-which(wtsum==min(wtsum))[1]
      final_des<-store_des[[min_value]]
      t1<-Sys.time()
      time_req<-t1-t0
      colnames(final_des)<-NULL
      lm=list("Latin hypercube design"=final_des,"Number of factors"=ncol(final_des),"Number of levels"=nrow(final_des),"Maxpro measure"=store_mp_phip_mac[min_value,1],"Phi_p measure"=store_mp_phip_mac[min_value,2],"Maximum absolute correlation"=store_mp_phip_mac[min_value,3],"Total system time"=time_req)
      return(lm)
    }
  }
}