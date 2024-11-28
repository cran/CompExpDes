OLHDs_2F<-function(levels){
  t0<-Sys.time()
  v=levels+1
  h_mat<-HadamardR::Hadamard_Matrix(v)[-1,-1]
  ###pos and neg
  for(i in 1:ncol(h_mat)){
    plus_one<-which(h_mat[,i]==1)
    neg_one<-which(h_mat[,i]==-1)
  }
  converted_mat<-NULL
  for(i in 1:ncol(h_mat)){
    plus_one<-which(h_mat[,i]==1)
    neg_one<-which(h_mat[,i]==-1)
    converted_mat<-cbind(converted_mat,c(plus_one,neg_one))
  }
  #write.table(converted_mat,"clipboard",sep="\t")
  pos<-(max(converted_mat)-1)/2
  neg<-(max(converted_mat)+1)/2
  #############
  col<-ncol(converted_mat)
  d<-converted_mat
  comb<-matrix(1:(col-1),ncol=2,byrow=TRUE)
  olhds<-list()
  for(i in 1:nrow(comb)){
    repeat{
      des<-rbind(apply(d[1:pos,c(comb[i,])],2,sample),apply(d[(pos+1):(pos+neg),c(comb[i,])],2,sample))
      if(MAC(des)==0 && Maxpro_Measure(des)<0.2){
        olhds<-append(olhds,list(des))
        break
      }
    }
  }
  maxpro<-lapply(olhds,function(mat)Maxpro_Measure(mat))
  index<-which(maxpro==min(unlist(maxpro)))
  #####################
  t1<-Sys.time()
  list<-list("OLHD"=olhds[[index]],"Factors"=2,"Runs"=v-1,"Required Time"=t1-t0)
  return(list)
}