library(beeswarm)
library(data.table)
library(alakazam)
library(RColorBrewer)
library(stringdist)
calculate = T
save_eps = T
fig_path = "/home/charles/Ellison_paper_figures/"

### Functions
group_color = function(patient) {
  colorpalette = brewer.pal(n=8,name="Paired")
  if (patient %in% pats_ellison_young_cmvn) return(colorpalette[4])
  if (patient %in% pats_ellison_young_cmvp) return(colorpalette[2])
  if (patient %in% pats_ellison_old_cmvn) return(colorpalette[8])
  if (patient %in% pats_ellison_old_cmvp) return(colorpalette[6])
}


### Organize patient information from a file "patients.tab"
pats_ellison_old_cmvp = ##### vector of the 5 relevant raw study participant IDs
pats_ellison_old_cmvn = ##### vector of the 5 relevant raw study participant IDs
pats_ellison_young_cmvp = ##### vector of the 5 relevant raw study participant IDs
pats_ellison_young_cmvn = ##### vector of the 5 relevant raw study participant IDs
pats_ellison = c(pats_ellison_old_cmvp,pats_ellison_old_cmvn,pats_ellison_young_cmvp,pats_ellison_young_cmvn)

path = "/media/charles/New_Volume/Ellison/Analysis/"
pat_data = read.table(paste(path,"patients.tab",sep=""),header=T,check.names=F)
pat_data = subset(pat_data,ID %in% pats_ellison)
pat_data$group = factor(sapply(pat_data$ID,function(s) {
  if (s %in% pats_ellison_young_cmvn) "young, CMV-"
  else if (s %in% pats_ellison_young_cmvp) "young, CMV+"
  else if (s %in% pats_ellison_old_cmvn) "elderly, CMV-"
  else if (s %in% pats_ellison_old_cmvp) "elderly, CMV+"
}),
levels = c("young, CMV-","young, CMV+","elderly, CMV-","elderly, CMV+"))
pat_data = pat_data[order(pat_data$group,pat_data$Age),]
pat_data$DummyID[pat_data$group == "young, CMV-"] = paste("YN",1:5,sep="")
pat_data$DummyID[pat_data$group == "young, CMV+"] = paste("YP",1:5,sep="")
pat_data$DummyID[pat_data$group == "elderly, CMV-"] = paste("EN",1:5,sep="")
pat_data$DummyID[pat_data$group == "elderly, CMV+"] = paste("EP",1:5,sep="")

id_to_dummy = pat_data$DummyID; names(id_to_dummy) = pat_data$ID
dummy_to_id = pat_data$ID; names(dummy_to_id) = pat_data$DummyID
ages_ellison = pat_data$Age; names(ages_ellison) = pat_data$DummyID


### Compute:
path = "/media/charles/New_Volume/Ellison_analysis_output/"
lineage_field = "LINEAGE_dissim0.1"
visitA = "V1"
visitB = "V3"
if (calculate) {
radius_increase = list()
for (pat in pats_ellison) {
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab",sep="")
  data0 = readChangeoDb(file)
  setnames(data0,old=lineage_field,new="LINEAGE")
  dataA = subset(data0,VISIT==visitA)
  dataB = subset(data0,VISIT==visitB)
  shared_lineages = intersect(dataA$LINEAGE,dataB$LINEAGE)
  dataA = subset(dataA, LINEAGE %in% shared_lineages)
  dataB = subset(dataB, LINEAGE %in% shared_lineages)
  radius_increase[[pat]] = sapply(shared_lineages,function(x) {
    seqsA = unique(subset(dataA,LINEAGE==x,JUNCTION,drop=T))
    distA = as.matrix(stringdistmatrix(seqsA,method="hamming"))
    radiusA = max(distA)
    seqsB = unique(subset(dataB,LINEAGE==x,JUNCTION,drop=T))
    seqsAB = union(seqsA,seqsB)
    distAB = as.matrix(stringdistmatrix(seqsAB,method="hamming"))
    radiusAB = max(distAB)
    return(radiusAB-radiusA)})
}
saveRDS(radius_increase,paste(fig_path,"radius_increase_corrected",sep=""))
}



### Make Fig. 3F
radius_increase = readRDS(paste(fig_path,"radius_increase_corrected",sep=""))
if (save_eps) {
  pdf(file=paste(fig_path,"Figure3F.eps",sep=""), 
      width=3.5,height=3.1,paper="a4",
      colormodel="rgb",pointsize = 10)
}
xrange=range(unlist(radius_increase))
for (i in 1:length(pats_ellison)) {
  pat = pats_ellison[i]
  h = hist(radius_increase[[pat]],breaks=seq(-0.5,xrange[2]+0.5,1),
           plot=F)
  if (i==1) {
    plot(h$mids,h$counts/sum(h$counts),
         xlim=c(0,10),ylim=c(0,0.85),
         type='l',col=group_color(pat),
         ann=F,axes=F,xaxt="n",yaxt="n",
         xlab="Lineage radius increase",
         ylab="Frequency")
  } else {
    lines(h$mids,h$counts/sum(h$counts),type='l',col=group_color(pat))
  }
}
mtext(side=1,line=1.5,"Lineage radius increase (day 0 to day 28)")
mtext(side=2,line=1.75,"Frequency")
axis(side=1,at=axTicks(1),mgp=c(3, .5, 0))
axis(side=2,at=axTicks(2),mgp=c(3, .5, 0))
box()
if (save_eps) dev.off()

