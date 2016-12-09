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
library('ape', lib.loc=package_path)

args = commandArgs(trailingOnly=TRUE)
file = args[1]
pat_id = strsplit(file,"_")[[1]][1]
germline_file = args[2]

data = readChangeoDb(file)
genotype_seqs = read.FASTA(paste(pat_id,"personal_hIGHV.fasta",sep="_"))
genotype_seqs = as.character(genotype_seqs)
genotype_seqs = sapply(genotype_seqs,function(s) paste(s,collapse=""))

V_CALL_GENOTYPED = reassignAlleles(data, genotype_seqs)
data = bind_cols(data, V_CALL_GENOTYPED)
data$V_CALL_GENOTYPED = sapply(strsplit(as.character(data$V_CALL_GENOTYPED),","),
                               function(s) paste(unique(s),collapse=","))

outfile = paste(strsplit(file,".tab")[[1]][1],"_V-CALL-GENOTYPED-pass.tab",sep="")
writeChangeoDb(data,file=outfile)