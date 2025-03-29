wtLHDs<-function(levels,factors,weights=c(0.2,0.4,0.4),iterations=300,population=100){
  k=factors
  n=levels
  p=2*k+1
  if(levels>choose(p,2) || levels<factors ||factors<2){
    return(message("Factors, F (>2) is a prime number and Levels, L should be in the range from F to (F^2)"))
  }
  rotation_mat<-function(matrix){
    a=as.matrix(matrix)
    ini=c(1:ncol(a))
    final=NULL
    for(i in 1:(ncol(a)-1)){
      x=ini-i
      final=rbind(final,x)
    }
    final=rbind(ini,final)
    final=final%%ncol(a)
    final[final==0]<-ncol(a)
    
    fmat=NULL
    for(j in 1:ncol(a)){
      output=NULL
      for(k in 1:ncol(a)){
        output=cbind(output,(a[,final[j,k]]))
      }
      fmat=rbind(fmat,output)
    }
    return(fmat) 
  }
  t0=Sys.time()
  v=choose(p,2)
  seq1=c(1:(p-1))
  seq2=c((p-1):2)
  list<-list()
  ######################
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
  partition1<-final_array[,1:((p-1)/2)]
  #partition2<-final_array[,(1+((p-1)/2)):(p-1)]
  #################################################
  des<-rotation_mat(partition1)
  ################ store regions
  reg_store<-list()
  for(a in 1:(nrow(des)/p)){
    reg_store<-append(reg_store,list(des[((a*p)-(p-1)):(a*p),])) 
  }
  #######Apply GA algorithm
  generate_permuted_matrices <- function(input_matrix, times) {
    permuted_list <- lapply(1:times, function(i) {
      apply(input_matrix, 2, sample)  # Shuffle each column independently
    })
    return(permuted_list)
  }
  ###############
  finallist<-list()
  for(b in 1:length(reg_store)){
    all_perm<-generate_permuted_matrices(reg_store[[b]],population)
    cal_mac<-lapply(all_perm,function(mat)MAC(mat))
    min_mac_pos<-which.min(cal_mac)
    finallist<-append(finallist,list(all_perm[[min_mac_pos]]))
  }
  ##########reconstruct design using reduced MAC
  des<-do.call(rbind,finallist)
  
  ######## delete as per number of levels
  #n=15
  ########
  store<-NULL
  store_des<-list()
  for(sim in 1:iterations){
    if(n==v){
      des1<-des
    }else{
      random_rows<-sample(1:v,v-n)
      des1<-des[-c(random_rows),]
      ##########random replacement code
      des2<-NULL
      for(c in 1:ncol(des1)){
        des2<-cbind(des2,c(setdiff(1:v,des1[,c])))
      }
      #now replace elements
      for(i in 1:ncol(des1)){
        now_ele<-c(des1[,i])
        del_ele<-c(des2[,i])
        need_replacement <-now_ele[now_ele>length(now_ele)]
        replacing_ele<-del_ele[del_ele<=length(now_ele)]
        if(length(replacing_ele)!=0){
          # Find positions
          positions <- which(des1[,i] %in% c(need_replacement))
          if(length(replacing_ele)!=1){
            des1[,i][c(positions)]<-sample(replacing_ele)
          }else{
            des1[,i][c(positions)]<-(replacing_ele)
          }
        }
      }
    }
    store_des<-append(store_des,list(des1))
    store<-rbind(store,c(MAC(des1),PhipMeasure(des1,q=1),Maxpro_Measure(des1)))
  }
  weight_given<-store%*%diag(weights)
  row_sum<-apply(weight_given,1,sum)
  find_pos<-which.min(row_sum)
  desired_des<-store_des[[find_pos]]
  LHD<-desired_des
  t1=Sys.time()
  time_req=t1-t0
  lm=list("LHD"=LHD,"Maxpro_Measure"=Maxpro_Measure(LHD),"Phi_p Measure"=PhipMeasure(LHD),"Maximum Absolute Correlation"=MAC(LHD),"Total System Time Requires"=time_req) 
  return(lm)
}

