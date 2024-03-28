UDesigns_2<-function(n){
  ########################
  if(n%%2==0 || n<5){
    return(message("Please provide an odd value of n (>=5)."))
  }
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
  final_list<-list("Uniform_Design"=final1,"Number of Factors"=n-1,"Number of Levels"=n*(n-1)/2,"Number of Runs"=nrow(final1),"Discrete Discrepancy Measure"=Discrete_Discrepancy(final1,1,0.3)[[1]],
                   "Lower Bound of Discrete Discrepancy"=Discrete_Discrepancy(final1,1,0.3)[[2]])
  return(final_list)
}

