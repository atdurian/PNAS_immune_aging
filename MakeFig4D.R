library(beeswarm)
library(data.table)
library(alakazam)
library(RColorBrewer)
library(plotrix)
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
types = c("VERY_SIMILAR","SIMILAR","DISSIMILAR","VERY_DISSIMILAR")
isotypes = c("IgM","IgD","IgG","IgA","IgE")
if (calculate) {
  mutation_types = array(dim=c(length(pats_ellison),length(visits),length(isotypes),length(types)),
                         dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,types=types))
  for (pat in pats_ellison) {
    file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
    data = readChangeoDb(file)
	
    #*#*#*#*# Ns lead to Xs in AA sequence, which are interpreted as "VERY DISSIMILAR" mutations by IMGT!
    data$VERY_DISSIMILAR = data$VERY_DISSIMILAR - data$X_COUNT_AA_SEQUENCE_VDJ
    #*#*#*#*#
    
	summary_data = aggregate(cbind(VERY_SIMILAR,SIMILAR,DISSIMILAR,VERY_DISSIMILAR) ~ VISIT+PRCONS,data,sum)
    for (vis in visits) {
      for (type in types) {
        for (iso in isotypes) {
          mutation_types[pat,vis,iso,type] = subset(summary_data,VISIT==vis & PRCONS==iso,type,drop=T)
        }
      }
    }
  }
  saveRDS(mutation_types,paste(fig_path,"mutation_types_byIsotype",sep=""))
}


### Figure 4D (radical mutations)
mutation_types = readRDS(paste(fig_path,"mutation_types_byIsotype",sep=""))
vis = "V2"
pats = c(pats_ellison_young_cmvn[order(id_to_dummy[pats_ellison_young_cmvn])],
         pats_ellison_young_cmvp[order(id_to_dummy[pats_ellison_young_cmvp])],
         pats_ellison_old_cmvn[order(id_to_dummy[pats_ellison_old_cmvn])],
         pats_ellison_old_cmvp[order(id_to_dummy[pats_ellison_old_cmvp])])
types = c("VERY_SIMILAR","SIMILAR","DISSIMILAR","VERY_DISSIMILAR")

isos = c("IgA","IgG")
M = t(as.matrix(apply(mutation_types[pats,vis,isos,types],c("patients","types"),sum)))
M = apply(M,2,function(s) s/sum(s))

select_types = c("VERY_DISSIMILAR","DISSIMILAR")
M_collapseType = apply(M[select_types,],"patients",sum)

wilcox.test(M_collapseType[c(pats_ellison_young_cmvn,pats_ellison_young_cmvp)],
            M_collapseType[c(pats_ellison_old_cmvn,pats_ellison_old_cmvp)])

data_list = list(100*M_collapseType[pats_ellison_young_cmvn],
                 100*M_collapseType[pats_ellison_young_cmvp],
                 100*M_collapseType[pats_ellison_old_cmvn],
                 100*M_collapseType[pats_ellison_old_cmvp])


if (save_eps) {
  pdf(file=paste(fig_path,"Figure4D.eps",sep=""), 
      width=3.5,height=2.55,paper="a4",
      colormodel="rgb",pointsize = 10)
}
par(mar=c(4,6,2,2))
color = brewer.pal(n=6,name="Set1")[4]
boxplot(data_list,
        names=rep("",4),
        ylab="",
        border=color,
        outline=F,range=Inf,boxwex=0.7)
beeswarm(data_list,"swarm",add=T,pch=20,cex=1.25)
mtext(side=2,line=3, "% radical mutations\n(isotype-switched sequences)")
mtext(side=1,line=1.5,at=1:4,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure4D_legend.eps",sep=""), 
      width=3.5,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("top","day 7",text.col=color,bty='n')
if (save_eps) dev.off()


