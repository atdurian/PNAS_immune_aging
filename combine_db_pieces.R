args = commandArgs(trailingOnly=TRUE)
files = args

package_path = "/local10G/debourcy/tools/Rpackages/"
library("ggplot2"); library("plyr"); library("reshape2"); library("scales")
library("lazyeval",lib.loc=package_path)
library("igraph",lib.loc=package_path)
library("alakazam",lib.loc=package_path)

data = data.frame()
for (file in files) {
data1 = readChangeoDb(file)
data = rbind(data,data1)
}

part = grep("part", strsplit(file,"_")[[1]], value = TRUE)
if (length(part)==0) {
  outfile = file
} else {
  outfile = gsub(paste("_",part,sep=""),"",file)
}
outfile = gsub("_db-pass","_all-chunks-db-pass",outfile)

writeChangeoDb(data,outfile)