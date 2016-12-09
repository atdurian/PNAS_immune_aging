library(beeswarm)
library(data.table)
library(alakazam)
library(RColorBrewer)
calculate = T
save_eps = T
fig_path = "/home/charles/Ellison_paper_figures/"


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


### Compute fraction of persistent lineages that become more abundant from V1 to V2:
path = "/media/charles/New_Volume/Ellison_analysis_output/"

if (calculate) {
n_responding_V1V2 = n_persistent_V1V2 = numeric()
lineage_field = "LINEAGE_dissim0.1"
for (pat in pats_ellison) {  
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass.tab",sep="")
  data = readChangeoDb(file)
  setnames(data,old=lineage_field,new="LINEAGE")
  data = aggregate(DUPCOUNT ~ LINEAGE+VISIT,data,sum) 
  
  df1 = subset(data,VISIT =="V1")
  df2 = subset(data,VISIT =="V2")
  df1 = df1[df1$DUPCOUNT < 0.05*sum(df1$DUPCOUNT),] # exclude the huge lineages
  df2 = df2[df2$DUPCOUNT < 0.05*sum(df2$DUPCOUNT),] # exclude the huge lineages
  
  cl1 = df1$LINEAGE # V1 lineages
  cl2 = df2$LINEAGE # V2 lineages
  cl = intersect(cl1,cl2) # persistent lineages
  
  df1 = subset(df1,LINEAGE %in% cl)
  rownames(df1) = paste("lineage",df1$LINEAGE,sep="_")
  df2 = subset(df2,LINEAGE %in% cl)
  rownames(df2) = paste("lineage",df2$LINEAGE,sep="_")
  
  cl_names = paste("lineage",cl,sep="_")
  n_responding_V1V2[pat] = sum(df2[cl_names,"DUPCOUNT"] > df1[cl_names,"DUPCOUNT"])
  n_persistent_V1V2[pat] = length(cl)   
}

### Save plot data:
saveRDS(n_responding_V1V2,paste(fig_path,"n_responding_V1V2",sep=""))
saveRDS(n_persistent_V1V2,paste(fig_path,"n_persistent_V1V2",sep=""))
}


### Make Fig. 3E
n_responding_V1V2 = readRDS(paste(fig_path,"n_responding_V1V2",sep=""))
n_persistent_V1V2 = readRDS(paste(fig_path,"n_persistent_V1V2",sep=""))
fraction_responding = n_responding_V1V2/n_persistent_V1V2

if (save_eps) {
  pdf(file=paste(fig_path,"Figure3E.eps",sep=""), 
      width=3.5,height=2.55,paper="a4",
      colormodel="rgb",pointsize = 10)
}
par(mar=c(4,6,2,2))
data_list = list(100*fraction_responding[pats_ellison_young_cmvn],
                 100*fraction_responding[pats_ellison_young_cmvp],
                 100*fraction_responding[pats_ellison_old_cmvn],
                 100*fraction_responding[pats_ellison_old_cmvp])
color = brewer.pal(n=6,name="Set1")[4]
boxplot(data_list,
        names=rep("",4),
        ylab="",
        border=color,
        outline=F,range=Inf,boxwex=0.7)
beeswarm(data_list,"swarm",add=T,pch=20,cex=1.25)
mtext(side=2,line=3, "% persistent lineages\nparticipating in response")
mtext(side=1,line=1.5,at=1:4,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
test = wilcox.test(100*fraction_responding[c(pats_ellison_young_cmvn,pats_ellison_old_cmvn)],
                   100*fraction_responding[c(pats_ellison_young_cmvp,pats_ellison_old_cmvp)])
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Fig3E_legend.eps",sep=""), 
      width=3.5,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("top","day 7",text.col=color,bty='n')
if (save_eps) dev.off()


