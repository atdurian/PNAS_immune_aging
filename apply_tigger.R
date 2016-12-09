package_path = "/local10G/debourcy/tools/Rpackages/"
library("ggplot2"); library("plyr"); library("reshape2"); library("scales")
library("lazyeval",lib.loc=package_path)
library("igraph",lib.loc=package_path)
library("alakazam",lib.loc=package_path)
library("data.table",lib.loc=package_path)
library("snow",lib.loc=package_path)
library("iterators",lib.loc=package_path)
library("foreach",lib.loc=package_path)
library("doSNOW",lib.loc=package_path)
library("SDMTools",lib.loc=package_path)
library("ade4", lib.loc=package_path)
library("seqinr",lib.loc=package_path)
library("shm",lib.loc=package_path)
library("Rcpp",lib.loc=package_path)
library("dplyr",lib.loc=package_path)
library("tigger",lib.loc=package_path)


args = commandArgs(trailingOnly=TRUE)
file = args[1]
pat_id = strsplit(file,"_")[[1]][1]
germline_file = args[2]
nproc = as.numeric(args[3])

data = readChangeoDb(file)


germline_ighv = read.fasta(germline_file,as.string=T)
germline_ighv = sapply(germline_ighv,function(s) toupper(s))

geno=inferGenotype(data,find_unmutated=F,germline_db=germline_ighv)
genotype_seqs = genotypeFasta(geno, germline_ighv)


########### Collapse germlines that are indistinguishable within sequenced region
a=genotype_seqs
a=gsub("-",".",a)
a=sapply(1:length(a),function(s) {
  allele = names(a)[s]
  sequence = a[s]
  if (grepl("IGHV1",allele)) {start=276}       # start position of sequence given the primers that were used
  if (grepl("IGHV2",allele)) {start=280}
  if (grepl("IGHV3",allele)) {start=284}
  if (grepl("IGHV4",allele)) {start=281}
  if (grepl("IGHV5",allele)) {start=279}
  if (grepl("IGHV6",allele)) {start=281}
  if (grepl("IGHV7",allele)) {start=288}
  return(substr(sequence,start,nchar(sequence)))  
})
map = character()
for (i in 1:length(a)) {
  allele = names(a)[i]
  seq = a[allele]
  map[allele] = sort(names(a)[a==seq])[1]
}
names(genotype_seqs)=map[names(genotype_seqs)]
###########

write.fasta(sequences=as.list(genotype_seqs), names=names(genotype_seqs), as.string=T,
            file.out=paste(pat_id,"_personal_hIGHV.fasta",sep=""))

V_CALL_GENOTYPED = reassignAlleles(data, genotype_seqs)
data = bind_cols(data, V_CALL_GENOTYPED)

data$V_CALL_GENOTYPED = sapply(strsplit(as.character(data$V_CALL_GENOTYPED),","),
                               function(s) paste(unique(map[s]),collapse=","))
outfile = paste(strsplit(file,".tab")[[1]][1],"_distinct-alleles-pass.tab",sep="")
writeChangeoDb(data,file=outfile)