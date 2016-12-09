args = commandArgs(trailingOnly=TRUE)
wdir = args[1]
file = args[2]

package_path = '/local10G/debourcy/tools/Rpackages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
removeBrackets = function(entries) as.numeric(sapply(entries,function(s) strsplit(s," ")[[1]][1]))

data = readChangeoDb(file)

### Compute mutations using the inferred germline for the V and J segments:
data$NMUT_GERMLINE = sapply(1:nrow(data), function(i) {
seq=strsplit(data$SEQUENCE_IMGT[i],"")[[1]]
germ=strsplit(data$GERMLINE_IMGT_D_MASK[i],"")[[1]]
sites = ((seq %in% c("A","C","T","G")) & (germ %in% c("A","C","T","G")))
cond = (sites & (seq!=germ))
return(sum(cond))
})

data$NMUT_PER_SITE_GERMLINE = sapply(1:nrow(data), function(i) {
seq=strsplit(data$SEQUENCE_IMGT[i],"")[[1]]
germ=strsplit(data$GERMLINE_IMGT_D_MASK[i],"")[[1]]
sites = ((seq %in% c("A","C","T","G")) & (germ %in% c("A","C","T","G")))
cond = (sites & (seq!=germ))
return(sum(cond)/sum(sites))
})


### Count mutations in the V segment reported by IMGT:
# Read IMGT output:
patient = strsplit(file,"_")[[1]][1]
visit = strsplit(file,"_")[[1]][2]
patient_visit = paste(patient,visit,sep="_")
IMGTdata = data.frame()
IMGTfolders = list.files(wdir,patient_visit)
IMGTfolders = IMGTfolders[grepl("_minCONSCOUNT2",IMGTfolders) & 
							!grepl("txz",IMGTfolders) & 
							!grepl("tab",IMGTfolders) & 
							!grepl("fasta",IMGTfolders) &
							!grepl("extracted",IMGTfolders)]
for (IMGTfolder in IMGTfolders) {
	IMGTfile = paste(wdir,"/",IMGTfolder,"/8_V-REGION-nt-mutation-statistics.txt",sep="")
	newIMGTdata = read.delim(IMGTfile,stringsAsFactors=F)
	newIMGTdata$SEQUENCE_ID = sapply(newIMGTdata$Sequence.ID, function(s) strsplit(s,"_")[[1]][1])
	newIMGTdata$SEQUENCE_VISIT_ID = paste(newIMGTdata$SEQUENCE_ID,visit,sep=";")
	IMGTdata = rbind(IMGTdata,newIMGTdata)
}

# Count mutations:
IMGTdata$NMUT_IMGT_V = removeBrackets(IMGTdata$V.REGION.a.g) + removeBrackets(IMGTdata$V.REGION.g.a) + 
removeBrackets(IMGTdata$V.REGION.c.t) + removeBrackets(IMGTdata$V.REGION.t.c) + 
removeBrackets(IMGTdata$V.REGION.a.c) + removeBrackets(IMGTdata$V.REGION.c.a) + 
removeBrackets(IMGTdata$V.REGION.a.t) + removeBrackets(IMGTdata$V.REGION.t.a) + 
removeBrackets(IMGTdata$V.REGION.g.c) + removeBrackets(IMGTdata$V.REGION.c.g) + 
removeBrackets(IMGTdata$V.REGION.g.t) + removeBrackets(IMGTdata$V.REGION.t.g)
IMGTdata$NMUT_PER_SITE_IMGT_V = IMGTdata$NMUT_IMGT_V/removeBrackets(IMGTdata$V.REGION.Nb.of.nucleotides)
IMGTdata = IMGTdata[,c("SEQUENCE_VISIT_ID","NMUT_IMGT_V","NMUT_PER_SITE_IMGT_V")]
# Bind mutation counts to the changeo database:
data$SEQUENCE_VISIT_ID = paste(data$SEQUENCE_ID,visit,sep=";")
data = merge(data,IMGTdata,by="SEQUENCE_VISIT_ID")

### Write result file:
outfile = gsub(".tab","_mutations-pass-corrected.tab",file)
writeChangeoDb(data,outfile)





