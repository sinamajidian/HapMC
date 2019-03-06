add="/SinaMc/University/codeMy/CodePaper1-UsingMatrixcompletion/data/ConvertingGeraci-For-MATLAB/exp-700-e0.2-c10-h0.4/Haplo-exp-700-e0.2-c10"
  

setwd(add)
Hall=matrix(NA,2*100,700)
for (ii in 1:100 ) {
  a=read.delim(paste("exp", ii-1,".h.AT", sep=""),header = F)
  h1=as.character(a[1,])
  h2=as.character(a[2,])
  for (i in 1:700 ) {
    if (substr(h1,i,i)=='t') {
      Hall[2*ii-1,i]=-1}
    if (substr(h1,i,i)=='a') {
      Hall[2*ii-1,i]=1}
    
    if (substr(h2,i,i)=='t') {
      Hall[2*ii,i]=-1}
    if (substr(h2,i,i)=='a') {
      Hall[2*ii,i]=1}
    }
}


H1=as.data.frame(Hall,col.names=F)
colnames(H1)=NULL

write.csv(H1,  file ="Haplotypes-e0.2-c8.txt", row.names = F)

