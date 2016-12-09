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


### Compute:
path = "/media/charles/New_Volume/Ellison_analysis_output/"
visits = c("V1","V2","V3")
types = c("AA_SEQUENCE_VDJ","AA_CDR3")
isotypes = c("IgD","IgM","IgA","IgG","IgE")
if (calculate) {
  stop_in_CDR3 = array(dim=c(length(pats_ellison),length(visits),length(isotypes),length(types)),
                       dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,types=types))
  for (pat in pats_ellison) {  
    for (vis in visits) {
      file = paste(path,pat,"_",vis,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab",sep="")
      data = readChangeoDb(file)
      data = subset(data,STOP==T & IN_FRAME==T)
      
      data$STOPCOUNT_AA_SEQUENCE_VDJ = nchar(data$AA_SEQUENCE_VDJ)-nchar(gsub("\\*","",data$AA_SEQUENCE_VDJ))
      data$STOPCOUNT_AA_CDR3 = nchar(data$AA_CDR3)-nchar(gsub("\\*","",data$AA_CDR3))
      data$STOP_INCIDENCE_AA_CDR3 =  as.numeric(data$STOPCOUNT_AA_CDR3 > 0)
      data$STOP_INCIDENCE_AA_SEQUENCE_VDJ =  as.numeric(data$STOPCOUNT_AA_SEQUENCE_VDJ > 0)
      
      summary_data = aggregate(cbind(STOP_INCIDENCE_AA_SEQUENCE_VDJ,STOP_INCIDENCE_AA_CDR3,STOPCOUNT_AA_SEQUENCE_VDJ,STOPCOUNT_AA_CDR3) ~ PRCONS,data,sum)
      
      for (iso in isotypes) {
        isodata = subset(summary_data,PRCONS==iso)
        if (nrow(isodata)==0) {
          stop_in_CDR3[pat,vis,iso,"AA_SEQUENCE_VDJ"] = 0
          stop_in_CDR3[pat,vis,iso,"AA_CDR3"] = 0
        } else {
          stop_in_CDR3[pat,vis,iso,"AA_SEQUENCE_VDJ"] = isodata$STOPCOUNT_AA_SEQUENCE_VDJ
          stop_in_CDR3[pat,vis,iso,"AA_CDR3"] = isodata$STOPCOUNT_AA_CDR3
        }
      }
    }
  }
  saveRDS(stop_in_CDR3,paste(fig_path,"stop_in_CDR3_inframe",sep=""))
}

### Make Fig. 4C
stop_in_CDR3 = readRDS(paste(fig_path,"stop_in_CDR3_inframe",sep=""))
visit = "V1"
y = 100*apply(stop_in_CDR3[pats_ellison,visit,isotypes,"AA_CDR3"],"patients",sum)/
  apply(stop_in_CDR3[pats_ellison,visit,isotypes,"AA_SEQUENCE_VDJ"],"patients",sum)

wilcox.test(y[c(pats_ellison_young_cmvn,pats_ellison_young_cmvp)],
            y[c(pats_ellison_old_cmvn,pats_ellison_old_cmvp)])

data_list = list(y[pats_ellison_young_cmvn],
                 y[pats_ellison_young_cmvp],
                 y[pats_ellison_old_cmvn],
                 y[pats_ellison_old_cmvp])

if (save_eps) {
  pdf(file=paste(fig_path,"Figure4C.eps",sep=""), 
      width=3.5,height=2.55,paper="a4",
      colormodel="rgb",pointsize = 10)
}
par(mar=c(4,6,2,2))
color=brewer.pal(n=6,name="Set1")[3]
boxplot(data_list,
        names=rep("",4),
        border=color,
        ylab="",boxwex=0.7,
        outline=F,range=Inf)
beeswarm(data_list,"swarm",add=T,pch=20,cex=1.25)
mtext(side=2,line=3, "% premature stop\ncodons in CDR3")
mtext(side=1,line=1.5,at=1:4,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure4C_legend.eps",sep=""), 
      width=3.5,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("top","day 0",text.col=color,bty='n')
if (save_eps) dev.off()

