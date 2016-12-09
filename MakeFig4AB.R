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
types = c("OUTFRAME","PREMATURE_STOP")
isotypes = c("IgD","IgM","IgA","IgG")
if (calculate) {
frac_unprod = array(dim=c(length(pats_ellison),length(visits),length(types)),
                            dimnames=list(patients=pats_ellison,visits=visits,types=types))
unprod_mutation_loads = array(dim=c(length(pats_ellison),length(visits),length(isotypes),length(types)),
                                    dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,types=types))
prod_mutation_loads = array(dim=c(length(pats_ellison),length(visits),length(isotypes)),
                            dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes))
total_unprod = array(dim=c(length(pats_ellison),length(visits),length(isotypes),length(types)),
                     dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,types=types))
total_prod = array(dim=c(length(pats_ellison),length(visits),length(isotypes)),
                     dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes))
mutation_field = "NMUT_IMGT_V" 
for (pat in pats_ellison) {  
  file_prod = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab",sep="")
  data_prod0 = readChangeoDb(file_prod)
  for (vis in visits) {
    data_prod = subset(data_prod0,VISIT==vis)
    file_unprod = paste(path,pat,"_",vis,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected_remove-error-clouds-pass.tab",sep="")
    data_unprod = readChangeoDb(file_unprod)
    data_unprod$OUTFRAME = !(data_unprod$IN_FRAME)
    data_unprod$PREMATURE_STOP = (data_unprod$IN_FRAME & data_unprod$STOP)
    total = nrow(data_prod)+nrow(data_unprod)
    frac_unprod[pat,vis,"OUTFRAME"] = sum(data_unprod$OUTFRAME)/total
    frac_unprod[pat,vis,"PREMATURE_STOP"] = sum(data_unprod$PREMATURE_STOP)/total
    for (iso in isotypes) { 
      unprod_mutation_loads[pat,vis,iso,"OUTFRAME"] = sum(subset(data_unprod,OUTFRAME==T & PRCONS==iso,mutation_field,drop=T))
      unprod_mutation_loads[pat,vis,iso,"PREMATURE_STOP"] = sum(subset(data_unprod,PREMATURE_STOP==T & PRCONS==iso,mutation_field,drop=T))
      prod_mutation_loads[pat,vis,iso] = sum(subset(data_prod,PRCONS==iso,mutation_field,drop=T),
                                             na.rm=T)
      total_unprod[pat,vis,iso,"OUTFRAME"] = nrow(subset(data_unprod,OUTFRAME==T & PRCONS==iso))
      total_unprod[pat,vis,iso,"PREMATURE_STOP"] = nrow(subset(data_unprod,PREMATURE_STOP==T & PRCONS==iso))
      total_prod[pat,vis,iso] = length(na.omit(subset(data_prod,PRCONS==iso,mutation_field,drop=T)))
    }
  }
}
saveRDS(frac_unprod,paste(fig_path,"frac_unprod_corrected",sep=""))
saveRDS(unprod_mutation_loads,paste(fig_path,"unprod_mutation_loads_corrected",sep=""))
saveRDS(prod_mutation_loads,paste(fig_path,"prod_mutation_loads_corrected",sep=""))
saveRDS(total_unprod,paste(fig_path,"total_unprod_corrected",sep=""))
saveRDS(total_prod,paste(fig_path,"total_prod_corrected",sep=""))
}

  
## Make Fig. 4A
frac_unprod = readRDS(paste(fig_path,"frac_unprod_corrected",sep=""))
if (save_eps) {
  pdf(file=paste(fig_path,"Figure4A.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
y = 100*frac_unprod[,,"PREMATURE_STOP"]
data_list = list(y[pats_ellison_young_cmvn,"V1"],
                 y[pats_ellison_young_cmvn,"V2"],
                 y[pats_ellison_young_cmvn,"V3"],
                 y[pats_ellison_young_cmvp,"V1"],
                 y[pats_ellison_young_cmvp,"V2"],
                 y[pats_ellison_young_cmvp,"V3"],
                 y[pats_ellison_old_cmvn,"V1"],
                 y[pats_ellison_old_cmvn,"V2"],
                 y[pats_ellison_old_cmvn,"V3"],
                 y[pats_ellison_old_cmvp,"V1"],
                 y[pats_ellison_old_cmvp,"V2"],
                 y[pats_ellison_old_cmvp,"V3"])
wilcox.test(unlist(data_list[c(1,4)]),unlist(data_list[c(7,10)]))
wilcox.test(unlist(data_list[c(2,5)]),unlist(data_list[c(8,11)]))
wilcox.test(unlist(data_list[c(3,6)]),unlist(data_list[c(9,12)]))
wilcox.test(unlist(data_list[1:6]),unlist(data_list[7:12]))

par(mar=c(8,5,2,2))
delta1 = 0.5
delta2 = 1
start = 2
box_pos = start+delta2*rep(0:3,each=3)+delta1*c(0:2,2:4,4:6,6:8)
name_pos = c(box_pos[2],box_pos[5],box_pos[8],box_pos[11])
colors = brewer.pal(n=6,name="Set1")[3:5]
boxplot(list(c(),c(),c(),c()),
        at=name_pos,
        las=1,
        names = rep("",4),
        ylab = "",axes=F,
        ylim=range(data_list))
boxplot(data_list,add=T,
        border = rep(colors,4),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",12),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=3,"% sequences with\npremature stop codons")
mtext(side=1,line=2,at=name_pos,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
axis(side=1,at=name_pos,label=F)
axis(side=2,at=axTicks(2),las=2,label=axTicks(2))
box()
if (save_eps) dev.off()


## Make Fig. 4B
total_unprod = readRDS(paste(fig_path,"total_unprod_corrected",sep=""))
unprod_mutation_loads = readRDS(paste(fig_path,"unprod_mutation_loads_corrected",sep=""))
total_prod = readRDS(paste(fig_path,"total_prod_corrected",sep=""))
prod_mutation_loads = readRDS(paste(fig_path,"prod_mutation_loads_corrected",sep=""))
if (save_eps) {
  pdf(file=paste(fig_path,"Figure4B.eps",sep=""), 
      width=3,height=3.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
visit = "V1"
y = (unprod_mutation_loads[pats_ellison,visit,"IgA","PREMATURE_STOP"]+
       unprod_mutation_loads[pats_ellison,visit,"IgG","PREMATURE_STOP"])/
  (total_unprod[pats_ellison,visit,"IgA","PREMATURE_STOP"]+
     total_unprod[pats_ellison,visit,"IgG","PREMATURE_STOP"])
x = (prod_mutation_loads[pats_ellison,visit,"IgA"]+
       prod_mutation_loads[pats_ellison,visit,"IgG"])/
  (total_prod[pats_ellison,visit,"IgA"]+
     total_prod[pats_ellison,visit,"IgG"])
data_list = list(y[c(pats_ellison_young_cmvn,pats_ellison_young_cmvp)],
                 x[c(pats_ellison_young_cmvn,pats_ellison_young_cmvp)],
                 y[c(pats_ellison_old_cmvn,pats_ellison_old_cmvp)],
                 x[c(pats_ellison_old_cmvn,pats_ellison_old_cmvp)])
testp = wilcox.test(y[c(pats_ellison_young_cmvn,pats_ellison_young_cmvp)],
            y[c(pats_ellison_old_cmvn,pats_ellison_old_cmvp)])$p.value
par(mar=c(8,5,2,2))
delta1 = 0.5
delta2 = 0.75
start = 1
box_pos = start+delta2*rep(0:1,each=2)+delta1*c(0:1,1:2)
name_pos = 0.2+c(mean(box_pos[1],box_pos[2]),
             mean(box_pos[3],box_pos[4]))
colors = brewer.pal(n=8,name="Set1")[1:2]
boxplot(list(c(),c()),
        at=name_pos,
        las=1,
        names = rep("",2),
        ylab = "",axes=F,
        ylim=c(range(data_list)[1],1.5+range(data_list)[2]))
boxplot(data_list,add=T,
        border = rep(colors,2),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",4),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=3,"Average number of mutations\n(isotype-switched sequences, day 0)")
mtext(side=1,line=1,at=name_pos,c("young","elderly"))
axis(side=1,at=name_pos,label=F)
axis(side=2,at=axTicks(2),las=2,label=axTicks(2))
box()
lines(c(box_pos[1],box_pos[3]),c(0.5+range(data_list)[2],0.5+range(data_list)[2]),col=colors[1])
text(mean(c(box_pos[1],box_pos[3])),c(1+range(data_list)[2],1+range(data_list)[2]),
     paste("p = ",round(testp,3),sep=""),col=colors[1])
if (save_eps) dev.off()
  
if (save_eps) {
  pdf(file=paste(fig_path,"Figure4B_legend.eps",sep=""), 
      width=10,height=10,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0.5,0.5,c("sequences\nwith\npremature\nstop codons\n",
                 "productive\nsequences"),
       text.col=colors,bty='n')
if (save_eps) dev.off()
