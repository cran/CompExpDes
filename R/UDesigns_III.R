###Gerechte designs for even numbers
UDesigns_III<-function(v){
  m=v/2
  if(v%%2!=0 || m<3){
    return(message("Please enter an even number, v >=6"))
  }
seq=c(1,v:2)
sudoku_mat<-NULL
for(i in seq){
  if(i%%2!=0){
    upto=i+(v-1)
    odd_seq<-(c(i,(i+1):upto))
    sudoku_mat<-cbind(sudoku_mat,odd_seq)
  }else{
    upto=i-(v-1)
    even_seq<-c(i,(i-1):upto)
    sudoku_mat<-(cbind(sudoku_mat,even_seq))
  }
}
sudoku_mat=sudoku_mat%%v
sudoku_mat[sudoku_mat==0]<-v
colnames(sudoku_mat)<-NULL
#############################UD
UD<-NULL
for(i in 1:(v/2)){
 UD<-rbind(UD,sudoku_mat[,((2*i)-(2-1)):(2*i)]) 
}
############# UD second
odd_numbers_rows<-seq(1,(v-1),by=2)
sudoku_mat_odd_rows<-sudoku_mat[odd_numbers_rows,]
################
UD1<-NULL
for(i in 1:(v/2)){
  UD1<-rbind(UD1,sudoku_mat_odd_rows[,((2*i)-(2-1)):(2*i)]) 
}
##########
unique1<-sort(unique(c(UD1[,1])))
unique2<-sort(unique(c(UD1[,2])))
ranks<-c(1:length(unique1))

for(i in 1:length(unique1)){
  UD1[,1][c(which(UD1[,1]==unique1[i]))]<-ranks[i]
}
for(i in 1:length(unique2)){
  UD1[,2][c(which(UD1[,2]==unique2[i]))]<-ranks[i]
}
################
third_col<-(UD1[,1]+UD1[,2])%%max(UD1)
third_col[third_col==0]<-max(UD1)
fourth_col<-(UD1[,1]+third_col)%%max(UD1)
fourth_col[fourth_col==0]<-max(UD1)
OA<-cbind(UD1,third_col,fourth_col)
colnames(OA)<-NULL
if(v>=10 && m%%2!=0){
  list1=list("Uniform Design 1"=UD,"Number of Factors"=2,"Number of Levels"=v,"Number of Runs"=(v*v)/2,"MAC"=MAC(UD),Discrete_Discrepancy(UD))
  list2=list("Uniform Design 2"=UD1,"Number of Factors"=2,"Number of Levels"=v/2,"Number of Runs"=(v*v)/4,"MAC"=MAC(UD1),Discrete_Discrepancy(UD1))
  list3=list("Uniform Design 3"=OA,"Number of Factors"=4,"Number of Levels"=v/2,"Number of Runs"=(v*v)/4,"MAC"=MAC(OA),Discrete_Discrepancy(OA))
}else{
  list1=list("Uniform Design 1"=UD,"Number of Factors"=2,"Number of Levels"=v,"Number of Runs"=(v*v)/2,"MAC"=MAC(UD))
  list2=list("Uniform Design 2"=UD1,"Number of Factors"=2,"Number of Levels"=v/2,"Number of Runs"=(v*v)/4,"MAC"=MAC(UD1))
  list3=list("Uniform Design 3"=OA,"Number of Factors"=4,"Number of Levels"=v/2,"Number of Runs"=(v*v)/4,"MAC"=MAC(OA)) 
}
all_list<-c(list1,list2,list3)
return(all_list)
}

