args = commandArgs(trailingOnly=TRUE)
file = args

package_path = '/local10G/debourcy/tools/Rpackages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
library('data.table',lib.loc=package_path)
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

CDR3_frac_dissim_thresholds = c(0.1,0.15,0.2)
data = readChangeoDb(file)
    
### group by V-J-length combination:
data$V_SEGMENT = extract_Vgene(data$V_CALL_GENOTYPED)
data$J_SEGMENT = sapply(data$J_CALL,function(s) strsplit(s," ")[[1]][2])
data$J_SEGMENT = extract_Jgene(data$J_SEGMENT)
data$VJL_GROUP = paste(data$V_SEGMENT,data$J_SEGMENT,data$JUNCTION_LENGTH,sep="_")

### cluster CDR3 sequences in each V-J-length group:
for (comb in unique(data$VJL_GROUP)) {    
  CDR3seqs = unique(data$JUNCTION[data$VJL_GROUP==comb]) 
  CDR3length = unique(data$JUNCTION_LENGTH[data$VJL_GROUP==comb])
  CDR3seqs = as.vector(CDR3seqs)
  hamming_d = as.matrix(stringdistmatrix(CDR3seqs,method="hamming"))
  rownames(hamming_d)=colnames(hamming_d)=CDR3seqs
  for (frac_threshold in CDR3_frac_dissim_thresholds) {
    if (length(CDR3seqs)==1) {
      clus = 1; names(clus)=CDR3seqs
    } else {
      hc = hclust(as.dist(hamming_d),method="single")
      clus = cutree(hc,h=frac_threshold*CDR3length)       
    }
    data[data$VJL_GROUP==comb,paste("CLUSTER",frac_threshold,sep="_dissim")] = 
      sapply(data$JUNCTION[data$VJL_GROUP==comb],
             function(s) paste(comb,clus[s],sep="_"))
  }
}

### change cluster names to simple numbers:
for (frac_threshold in CDR3_frac_dissim_thresholds) { 
  cluster_column = paste("CLUSTER",frac_threshold,sep="_dissim")
  lineage_column = paste("LINEAGE",frac_threshold,sep="_dissim")
  data[,cluster_column] = as.factor(data[,cluster_column])
  levels(data[,cluster_column]) = 1:length(levels(data[,cluster_column]))
  data[,cluster_column] = as.numeric(data[,cluster_column])
  setnames(data,old=cluster_column,new=lineage_column)
  data = data[order(data[,lineage_column]),]
}

### write result file:
outfile = gsub(".tab","_lineage-pass.tab",file)
writeChangeoDb(data,outfile)
 


