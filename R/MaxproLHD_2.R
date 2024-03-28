
MaxproLHD_2<-function(s){
#########
  if(s%%2==0 || s<5){
    return(message("Please provide an odd value of s (>=5)."))
  }
  n=s
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
  
  
  
  ##############################################
v=choose(n,2)
seq1=c(1:(n-1))
seq2=c((n-1):2)
list<-list()
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
###############
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
p1<-rotation_mat(partition1)
p2<-rotation_mat(partition2)
final<-cbind(p1,p2)
########################
final1=p1
final1p1=final1[1:n,]
final1p2<-final1[(nrow(final1p1)+1):v,]
if(n>5){
for(i in 1:50000){
  deg<-rbind(final1p1,apply(final1p2,2,sample))
  if( Maxpro_Measure(deg)<0.1 && MAC(deg)<0.2){
    break
  }
}
}else{
  for(i in 1:50000){
    deg<-rbind(final1p1,apply(final1p2,2,sample))
    if(MAC(deg)<0.1){
      break
    }
  }
}
###########2nd time
if(n>5){
for(i in 1:50000){
  deg11<-rbind(apply(deg[1:n,],2,sample),deg[(nrow(final1p1)+1):v,])
  if( Maxpro_Measure(deg11)<0.1 && MAC(deg11)<MAC(deg)){
    #final1<-deg
    break
  } 
}
}else{
  deg11=deg
}
############
final2=p2
final2p1=final2[1:n,]
final2p2<-final2[(nrow(final2p1)+1):v,]
if(n>5){
for(i in 1:50000){
  deg2<-rbind(final2p1,apply(final2p2,2,sample))
  if( Maxpro_Measure(deg2)<0.1 && MAC(deg2)<0.2){
    break
  }
}
}else{
  for(i in 1:50000){
    deg2<-rbind(final2p1,apply(final2p2,2,sample))
    if(MAC(deg2)<0.1){
      break
    }
  } 
}
#########
###########2nd time
if(n>5){
for(j in 1:50000){
  deg22<-rbind(apply(deg2[1:n,],2,sample),deg2[(nrow(final1p1)+1):v,])
  if( Maxpro_Measure(deg22)<0.1 && MAC(deg22)<MAC(deg2)){
    #final1<-deg
    break
  } 
}
}else{
  deg22=deg2
}
###############
#new_deg=cbind(deg11,deg22)
new_deg=final
pt1=new_deg[1:n,]
pt2=new_deg[(nrow(pt1)+1):nrow(new_deg),]
if(n>5){
for(i in 1:20000){
  f<-rbind(pt1,apply(pt2,2,sample))
  if(Maxpro_Measure(f)<=0.1){
    break
  }
}
}else{
  for(i in 1:20000){
  f<-rbind(pt1,apply(pt2,2,sample))
  if(MAC(f)<0.1){
    break
  }
  }
}
###################################################################
  list12=list("Maxpro_LHD1"=deg11,"Maxpro Measure"=Maxpro_Measure(deg11),
                      "Phi_p Measure"=PhipMeasure(deg11),"Maximum Absolute Correlation"=MAC(deg11),
              "Maxpro_LHD2"=deg22,"Maxpro Measure"=Maxpro_Measure(deg22),
              "Phi_p Measure"=PhipMeasure(deg22),"Maximum Absolute Correlation"=MAC(deg22),
              "Maxpro_LHD3"=f,"Maxpro Measure"=Maxpro_Measure(f), "Phi_p Measure"=PhipMeasure(f),"Maximum Absolute Correlation"=MAC(f))
  return(list12)
    
  }

#########

