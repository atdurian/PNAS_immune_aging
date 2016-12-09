args = commandArgs(trailingOnly=TRUE)
files = args

package_path = '/local10G/debourcy/tools/Rpackages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)

data = data.frame()
for (file in files) {
visit = strsplit(file,'_')[[1]][2]
data1 = readChangeoDb(file)
data1$VISIT = visit
data = rbind(data,data1)
}


data = data[order(data$VISIT),]
outfile = gsub(paste("_",visit,sep=""),"",file)
outfile = gsub(".tab","_visits-pooled.tab",outfile)

writeChangeoDb(data,outfile)