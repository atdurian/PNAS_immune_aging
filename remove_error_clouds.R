args = commandArgs(trailingOnly=TRUE)
file = args

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

data0 = readChangeoDb(file)

### Group by V-J-length combination:
data0$V_SEGMENT = extract_Vgene(data0$V_CALL_GENOTYPED)
data0$J_SEGMENT = sapply(data0$J_CALL,function(s) strsplit(s," ")[[1]][2])
data0$J_SEGMENT = extract_Jgene(data0$J_SEGMENT)
data0$VJL_GROUP = paste(data0$V_SEGMENT,data0$J_SEGMENT,data0$JUNCTION_LENGTH,sep="_")

### Remove sequences that are 1 bp substitution removed from highly abundant sequence within the same V-J-Length group:
data_out = data.frame()
for (visit in unique(data0$VISIT)) {
  data = subset(data0,VISIT==visit)
  highly_abundant_seqs = subset(data,DUPCOUNT>=20)
  for (i in 1:nrow(data)) {
  ref_seqs = subset(highly_abundant_seqs,VJL_GROUP==data$VJL_GROUP[i])
  if (!(data$SEQUENCE_VDJ[i] %in% ref_seqs$SEQUENCE_VDJ)) {
    D = stringdist(data$SEQUENCE_VDJ[i],ref_seqs$SEQUENCE_VDJ,method="hamming")
    if (any(D==1)) data[i,"DUPCOUNT"]=0
  }
  }
  data = data[data$DUPCOUNT!=0,]
  data_out = rbind(data_out,data)
}

### Write result file:
outfile = gsub(".tab","_remove-error-clouds-pass.tab",file)
writeChangeoDb(data_out,outfile)
