###UDesigns_2 start
UDesigns_II<-function(factors){
  t0=Sys.time()
  n=factors+1
  ########################
  # if(n%%2==0 || n<5){
  #   return(message("Please provide an odd value of factors (>=5)."))
  # }
  ###################################3
  v=choose(n,2)
  tri<-matrix(NA,n,n)
  k=1
  i=1
  while(i <n){
    for(j in 1:(n-i)){
      tri[i,j+i]<-k
      k=k+1
    }
    i=i+1
  }
  for(l in 1:n){
    tri[,l]<-tri[l,]
  }
  ############
  required_array=t(apply(tri,1,sort))
  ########################
  reg_col=ncol(required_array)/2
  reg_row=n
  first_region=required_array[,1:reg_col]
  second_region=required_array[,(1+reg_col):(n-1)]
  blank_reg=matrix(NA,n,reg_col)
  #####################
  for(j in 1:reg_row){
    if(length(setdiff(first_region[j,],blank_reg[1:j,]))==reg_col){
      blank_reg[j,]<-first_region[j,]
    }else{
      #k=1
      for(k in first_region[j,]){
        if(k %in% blank_reg[1:j,]){
          pos_repeated_element_in_first_region=which(first_region[j,]==k)
          same_pos_in_second_reg=pos_repeated_element_in_first_region
          now<-first_region[j,pos_repeated_element_in_first_region]
          next1<-second_region[j,pos_repeated_element_in_first_region]
          first_region[j,pos_repeated_element_in_first_region]<-next1
          second_region[j,pos_repeated_element_in_first_region]<-now
        }
      }
      blank_reg[1:j,]<-first_region[1:j,]
    }
  }
  ##################
  final=cbind(first_region,second_region)
  ###############
  ration_seq=c()
  seq=c(1:reg_col)
  for(i in 0:(reg_col-1)){
    ration_seq<-c(ration_seq,(seq+i)%%reg_col)
  }
  ration_seq[ration_seq==0]<-reg_col
  ########################
  #ration_seq=ration_seq[-c(1:reg_col)]
  ration_seq=matrix(ration_seq,(reg_col),reg_col,byrow=TRUE)
  rotation_seq_we_need=((ration_seq[2:nrow(ration_seq),]))
  if(n==5){
    rotation_seq_we_need=(t(rotation_seq_we_need))
  }
  ###############
  f1<-first_region
  for(i in 1:nrow((rotation_seq_we_need))){
    f1<-rbind(f1,first_region[,rotation_seq_we_need[i,]])
  }
  #######
  f2<-second_region
  for(i in 1:nrow(rotation_seq_we_need)){
    f2<-rbind(f2,second_region[,rotation_seq_we_need[i,]])
  }
  #################
  final1=cbind(f1,f2)
  #######################
  ration_seq1=c(ncol(final1):1)
  ######################final task
  final1<-rbind(final1,final1[,ration_seq1])
  DD<-Discrete_Discrepancy(final1)
  t1=Sys.time()
  total_time<-t1-t0
  final_list<-list("Uniform design"=final1,"Number of factors"=n-1,"Number of levels"=n*(n-1)/2,"Number of runs"=nrow(final1),"Maximum absolute correlation"=MAC(final1),"Discrete discrepancy measure"=DD[[1]],
                   "Lower bound of discrete discrepancy deasure"=DD[[2]],"Total system time"=total_time)
  return(final_list)
}

###Udesigns_II End
