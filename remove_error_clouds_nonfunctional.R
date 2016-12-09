args = commandArgs(trailingOnly=TRUE)
fileF = args[1]
visit = strsplit(fileF,"_")[[1]][2]
fileT = args[2]

package_path = '/local10G/debourcy/tools/Rpackages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
library('stringdist', lib.loc=package_path)

extract_Vgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  s=sapply(w,function(x) {paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
  return(as.vector(s))
}
extract_Jgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  return(as.vector(w))
}


dataF = readChangeoDb(fileF)
dataT = readChangeoDb(fileT)
dataT = subset(dataT,VISIT==visit)

### group by V-J-length combination:
dataF$V_SEGMENT = extract_Vgene(dataF$V_CALL_GENOTYPED)
dataF$J_SEGMENT = sapply(dataF$J_CALL,function(s) strsplit(s," ")[[1]][2])
dataF$J_SEGMENT = extract_Jgene(dataF$J_SEGMENT)
dataF$VJL_GROUP = paste(dataF$V_SEGMENT,dataF$J_SEGMENT,dataF$JUNCTION_LENGTH,sep="_")

dataT$V_SEGMENT = extract_Vgene(dataT$V_CALL_GENOTYPED)
dataT$J_SEGMENT = sapply(dataT$J_CALL,function(s) strsplit(s," ")[[1]][2])
dataT$J_SEGMENT = extract_Jgene(dataT$J_SEGMENT)
dataT$VJL_GROUP = paste(dataT$V_SEGMENT,dataT$J_SEGMENT,dataT$JUNCTION_LENGTH,sep="_")

data = rbind(dataF,dataT[,names(dataF)])

### Remove sequences that are 1 bp substitution removed from highly abundant sequence within the same V-J-Length group:
highly_abundant_seqs = subset(data,DUPCOUNT>=20)
for (i in 1:nrow(dataF)) {
ref_seqs = subset(highly_abundant_seqs,VJL_GROUP==dataF$VJL_GROUP[i])
if (!(dataF$SEQUENCE_VDJ[i] %in% ref_seqs$SEQUENCE_VDJ)) {
  D = stringdist(dataF$SEQUENCE_VDJ[i],ref_seqs$SEQUENCE_VDJ,method="hamming")
  if (any(D==1)) dataF[i,"DUPCOUNT"]=0
}
}
dataF = dataF[dataF$DUPCOUNT!=0,]

### Write result file:
outfile = gsub(".tab","_remove-error-clouds-pass.tab",fileF)
writeChangeoDb(dataF,outfile)
