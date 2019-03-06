
add="/SinaMc/University/codeMy/CodePaper1-UsingMatrixcompletion/data/ConvertingGeraci-For-MATLAB/exp-700-e0.2-c10-h0.4/read-exp-700-e0.2-c10/"
  
  
setwd(add)
R=matrix(0,1990*100,700) #c=8 > 1550  c=5>1000 c=10 1990
lengAll=matrix(0,100,1)
rowIDX=0;
#start.time <- Sys.time()
for (ii in 1:100 ) {
  print(ii)
  data_R=read.delim(paste("exp", ii-1,".m.AT", sep=""),header = F)
  leng=dim(data_R)[1]
  lengAll[ii]=leng
  for (i in 1:leng ) {
    rowIDX=rowIDX+1
    row=as.character(data_R[i,])
    row1=strsplit(row,'')
    b2=sub("a","1",row1[[1]])
    b3=sub("-","0",b2)
    b4=sub("t","-1",b3)
    b5=as.numeric(b4) 
    R[rowIDX,] =b5
  }
}


#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken
R1=as.data.frame(R[1:rowIDX,],col.names=F)
colnames(R1)=NULL

write.csv(R1,  file ="ReadMatrix-e0.2-c10.txt", row.names = F)
write.csv(lengAll,  file ="NumberRowsofMatrix.txt")
