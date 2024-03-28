MaxproLHD_1<-function(prime_number){
  
 is.prime<-function(value){
  for(i in 2:(value/2)){ 
    k="green"
  if(value%%i==0){
    k="red"
    return(FALSE)
  }
    if(k=="green"){
      return(TRUE)
    }
  }
 }
###############
 if(is.prime(prime_number)==FALSE || prime_number<3){
   return(message("Please enter a prime number (>2)"))
 }else{
  s=prime_number
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
  #######
  Maxpro_Measure<-function(Design){
    matrix<-as.matrix(Design)
    sumsq<-NULL
    for(i in 1:(nrow(matrix)-1)){
      j=i+1
      while(j<=nrow(matrix)){
        sumsq=c(sumsq,(1/prod(((matrix[i,]-matrix[j,])^2))))
        j=j+1
      }
    }
    return((sum(sumsq)/choose(nrow(matrix),2))^(1/ncol(matrix)))
  }
  ##########
  PhipMeasure<-function(design,p=15,q=2){
    design=as.matrix(design)
    EucliDist <- function(vec1, vec2) {
      return(sqrt(sum((vec1 - vec2)^2)))
    }
    ############
    phip=c()
    for(i in 1:(nrow(design)-1)){
      phip=c(phip,(EucliDist(design[i,],design[i+1,]))^(-p))
    }
    ##########
    phip=(sum(phip))^(1/p)
    return(phip)
  }
##########
  MAC<-function(matrix){
    matrix<-as.matrix(matrix)
    value<-sort(unique(abs(c(cor(matrix)))),decreasing = TRUE)[2]
    return((value))
  }
  
  #####################
  list=list("Maxpro_LHD"=LHD,"Maxpro Measure"=Maxpro_Measure(LHD),"Phi_p Measure"=PhipMeasure(LHD),"Maximum Absolute Correlation"=MAC(LHD))
return(list)
 }
}


