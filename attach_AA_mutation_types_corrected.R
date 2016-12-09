args = commandArgs(trailingOnly=TRUE)
wdir = args[1]
file = args[2]

package_path = '/local10G/debourcy/tools/Rpackages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
library('data.table',lib.loc=package_path)
removeBrackets = function(entries) as.numeric(sapply(entries,function(s) strsplit(s," ")[[1]][1]))
parseAAchangeStats = function(entry) {
  numbers = strsplit(entry,"\\+")
  numbers = lapply(numbers,function(n) as.numeric(sapply(n,function(s) strsplit(s," ")[[1]][1])))
  result = lapply(numbers,sum)
  return(unlist(result))
}

data = readChangeoDb(file)
								
## Read IMGT output:
patient = strsplit(file,"_")[[1]][1]
IMGTdata = data.frame()
for (vis in c("V1","V2","V3")) {
	IMGTfolders = list.files(wdir,paste(patient,vis,sep="_"))
	IMGTfolders = IMGTfolders[grepl("_minCONSCOUNT2",IMGTfolders) & 
								!grepl("txz",IMGTfolders) & 
								!grepl("tab",IMGTfolders) & 
								!grepl("fasta",IMGTfolders) &
								!grepl("extracted",IMGTfolders)]
	for (IMGTfolder in IMGTfolders) {
		IMGTfile = paste(wdir,"/",IMGTfolder,"/9_V-REGION-AA-change-statistics.txt",sep="")
		newIMGTdata = read.delim(IMGTfile,stringsAsFactors=F)
		newIMGTdata$SEQUENCE_ID = sapply(newIMGTdata$Sequence.ID, function(s) strsplit(s,"_")[[1]][1])
		newIMGTdata$SEQUENCE_VISIT_ID = paste(newIMGTdata$SEQUENCE_ID,vis,sep=";")
		IMGTdata = rbind(IMGTdata,newIMGTdata)
	}
}

## Get desired fields:
IMGTdata$VERY_SIMILAR = removeBrackets(IMGTdata$V.REGION.Very.similar)
IMGTdata$SIMILAR = parseAAchangeStats(IMGTdata$V.REGION.Similar)
IMGTdata$DISSIMILAR = parseAAchangeStats(IMGTdata$V.REGION.Dissimilar)
IMGTdata$VERY_DISSIMILAR = removeBrackets(IMGTdata$V.REGION.Very.dissimilar)    

IMGTdata$ppp = removeBrackets(IMGTdata$V.REGION....)
IMGTdata$ppm = removeBrackets(IMGTdata$V.REGION.....1)
IMGTdata$pmp = removeBrackets(IMGTdata$V.REGION.....2)
IMGTdata$pmm = removeBrackets(IMGTdata$V.REGION.....3)    
IMGTdata$mpm = removeBrackets(IMGTdata$V.REGION.....4)    
IMGTdata$mmp = removeBrackets(IMGTdata$V.REGION.....5)    
IMGTdata$mmm = removeBrackets(IMGTdata$V.REGION.....6)    

IMGTdata = IMGTdata[,c("SEQUENCE_VISIT_ID","VERY_SIMILAR","SIMILAR","DISSIMILAR","VERY_DISSIMILAR",
					 "ppp","ppm","pmp","pmm","mpm","mmp","mmm")]

## Bind fields to the changeo database:
data$SEQUENCE_VISIT_ID = paste(data$SEQUENCE_ID,data$VISIT,sep=";")
data = merge(data,IMGTdata,by="SEQUENCE_VISIT_ID")

## Compute number of Xs in AA_SEQUENCE_VDJ field:
data$X_COUNT_AA_SEQUENCE_VDJ = nchar(data$AA_SEQUENCE_VDJ)-nchar(gsub("X","",data$AA_SEQUENCE_VDJ))

## Write result file:
outfile = gsub(".tab","_attach-AA-mutationtypes-pass-corr.tab",file)
writeChangeoDb(data,outfile)


