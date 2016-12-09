args = commandArgs(trailingOnly=TRUE)
file = args[2]
name = args[1]


package_path = "/local10G/debourcy/tools/Rpackages/"
library("ade4",lib.loc=package_path)
library("seqinr",lib.loc=package_path)

   
data = read.fasta(file)

# Trim isotype designations:
names(data) = gsub("IgE_1","IgE",names(data)) 
names(data) = gsub("IgE_2","IgE",names(data))
data = lapply(data,function(s) replace(s,s=="-","n"))

# Write results:
n_seqs = length(data)
max_size = 5e5 # maximum number of sequences allowed per file uploaded to IMGT/hIGHV-QUEST
if (n_seqs > max_size) {
 n_chunks = ceiling(n_seqs/max_size)
 for (i in 1:n_chunks) {
   start = (i-1)*max_size+1
   end = min(i*max_size,n_seqs)
   part_data = data[start:end]
   write.fasta(sequences=part_data, names=names(part_data), nbchar=80,
               file.out=paste(name,"_minCONSCOUNT2_part",i,".fasta",sep=""))
 }
} else {
 write.fasta(sequences=data, names=names(data), nbchar=80,
               file.out=paste(name,"_minCONSCOUNT2.fasta",sep=""))
}