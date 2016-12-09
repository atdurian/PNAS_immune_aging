library(beeswarm)
library(data.table)
library(alakazam)
library(RColorBrewer)
library(plotrix)
save_eps = T
calculate = T
fig_path = "/home/charles/Ellison_paper_figures/"

### Functions
group_color = function(patient) {
  colorpalette = brewer.pal(n=8,name="Paired")
  if (patient %in% pats_ellison_young_cmvn) return(colorpalette[4])
  if (patient %in% pats_ellison_young_cmvp) return(colorpalette[2])
  if (patient %in% pats_ellison_old_cmvn) return(colorpalette[8])
  if (patient %in% pats_ellison_old_cmvp) return(colorpalette[6])
}

oligo_color = function(patient) {
  colorpalette = brewer.pal(n=8,name="Paired")
  if (patient %in% young) return(colorpalette[4])
  if (patient %in% nonoligoclonal_elderly) return(colorpalette[8])
  if (patient %in% oligoclonal_elderly) return(colorpalette[6])
}

visit_to_day = function(visit) {
  if (visit=="V1") return("day 0")
  if (visit=="V2") return("day 7")
  if (visit=="V3") return("day 28")
  if (visit=="all") return("all time points")
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

### Compute mutation loads by lineage and isotype compartment
path = "/media/charles/New_Volume/Ellison_analysis_output/"
lineage_field = "LINEAGE_dissim0.1"
visits = c("V1","V2","V3","all")
isotypes = c("IgD","IgM","IgA","IgG","IgE")
variables = c("SEQCOUNT","DUPCOUNT","NMUT_PER_SITE_GERMLINE")
if (calculate) {
results = array(dim=c(length(pats_ellison),length(visits),length(isotypes),2,length(variables)),
                 dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,lineages=c("top-lineage",'rest-repertoire'),variables=variables))
frac_IgD_also_IgM = array(dim=c(length(pats_ellison),length(visits),2),
                          dimnames=list(patients=pats_ellison,visits=visits,lineages=c("top-lineage",'rest-repertoire')))
for (pat in pats_ellison) {
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab",sep="")
  data0 = readChangeoDb(file)
  setnames(data0,old=lineage_field,new="LINEAGE")
  data0$SEQCOUNT = 1
  formula1 = as.formula(paste("cbind(",paste(variables,collapse=","),") ~ LINEAGE+VISIT+PRCONS",sep=""))
  summary_data = aggregate(formula1, data=data0, FUN=sum) 
  
  # determine top lineages
  lineage_abundances = aggregate(DUPCOUNT ~ LINEAGE+VISIT, data=summary_data, FUN=sum) 
  lineage_abundances = lineage_abundances[order(lineage_abundances$VISIT,-lineage_abundances$DUPCOUNT),]
  lineage_abundances_overall = aggregate(DUPCOUNT ~ LINEAGE, data=lineage_abundances, FUN=sum) 
  lineage_abundances_overall = lineage_abundances_overall[order(-lineage_abundances_overall$DUPCOUNT),]
  top_lineage_V1 = lineage_abundances$LINEAGE[lineage_abundances$VISIT=="V1"][1]
  top_lineage_V2 = lineage_abundances$LINEAGE[lineage_abundances$VISIT=="V2"][1]
  top_lineage_V3 = lineage_abundances$LINEAGE[lineage_abundances$VISIT=="V3"][1]
  top_lineage_overall = lineage_abundances_overall$LINEAGE[1]
  
  # summarize data by "top lineage" versus "rest of repertoire"
  summary_data_overall = summary_data
  summary_data_overall$TOP_LINEAGE_OVERALL = (summary_data_overall$LINEAGE==top_lineage_overall)
  formula3 = as.formula(paste("cbind(",paste(variables,collapse=","),") ~ PRCONS+TOP_LINEAGE_OVERALL",sep=""))
  summary_data_overall = aggregate(formula3, data=summary_data_overall, FUN=sum) 
  
  summary_data$TOP_LINEAGE_AT_VISIT = (summary_data$VISIT=="V1" & summary_data$LINEAGE==top_lineage_V1) |
                                      (summary_data$VISIT=="V2" & summary_data$LINEAGE==top_lineage_V2) |
                                      (summary_data$VISIT=="V3" & summary_data$LINEAGE==top_lineage_V3)
  formula2 = as.formula(paste("cbind(",paste(variables,collapse=","),") ~ VISIT+PRCONS+TOP_LINEAGE_AT_VISIT",sep=""))
  summary_data = aggregate(formula2, data=summary_data, FUN=sum) 
  
  # make sure missing isotypes in the top lineage are represented with the appropriate value
  summary_data$VISIT = as.factor(summary_data$VISIT)
  summary_data$PRCONS = as.factor(summary_data$PRCONS)
  summary_data$TOP_LINEAGE_AT_VISIT = as.factor(summary_data$TOP_LINEAGE_AT_VISIT)
  all_factor_combinations = with(summary_data,expand.grid(VISIT=levels(VISIT),PRCONS=levels(PRCONS),TOP_LINEAGE_AT_VISIT=levels(TOP_LINEAGE_AT_VISIT)))
  summary_data = merge(summary_data,all_factor_combinations,all.y=T)
  summary_data$SEQCOUNT[is.na(summary_data$SEQCOUNT)] = 0
  summary_data$DUPCOUNT[is.na(summary_data$DUPCOUNT)] = 0
  
  summary_data_overall$PRCONS = as.factor(summary_data_overall$PRCONS)
  summary_data_overall$TOP_LINEAGE_OVERALL = as.factor(summary_data_overall$TOP_LINEAGE_OVERALL)
  all_factor_combinations = with(summary_data_overall,expand.grid(PRCONS=levels(PRCONS),TOP_LINEAGE_OVERALL=levels(TOP_LINEAGE_OVERALL)))
  summary_data_overall = merge(summary_data_overall,all_factor_combinations,all.y=T)
  summary_data_overall$SEQCOUNT[is.na(summary_data_overall$SEQCOUNT)] = 0
  summary_data_overall$DUPCOUNT[is.na(summary_data_overall$DUPCOUNT)] = 0

  # record results 
  for (vis in visits) {
    for (iso in isotypes) {
      for (var in variables) {
        if (vis=="all") {
          results[pat,vis,iso,"top-lineage",var] = subset(summary_data_overall,PRCONS==iso & TOP_LINEAGE_OVERALL==T,var,drop=T)
          results[pat,vis,iso,"rest-repertoire",var] = subset(summary_data_overall,PRCONS==iso & TOP_LINEAGE_OVERALL==F,var,drop=T)
        } else {
          results[pat,vis,iso,"top-lineage",var] = subset(summary_data,VISIT==vis & PRCONS==iso & TOP_LINEAGE_AT_VISIT==T,var,drop=T)
          results[pat,vis,iso,"rest-repertoire",var] = subset(summary_data,VISIT==vis & PRCONS==iso & TOP_LINEAGE_AT_VISIT==F,var,drop=T)
        }
      }
    }
  }
  
  # compute fraction of IgD sequences also found in IgM
  for (vis in visits) {
    if (vis=="all") {
      top_lineage = top_lineage_overall
      seqs_D_top = subset(data0,LINEAGE==top_lineage & PRCONS=="IgD",SEQUENCE_VDJ,drop=T)
      seqs_M_top = subset(data0,LINEAGE==top_lineage & PRCONS=="IgM",SEQUENCE_VDJ,drop=T)
      seqs_D_rest = subset(data0,LINEAGE!=top_lineage & PRCONS=="IgD",SEQUENCE_VDJ,drop=T)
      seqs_M_rest = subset(data0,LINEAGE!=top_lineage & PRCONS=="IgM",SEQUENCE_VDJ,drop=T)
      frac_IgD_also_IgM[pat,vis,"top-lineage"] = length(intersect(seqs_D_top,seqs_M_top))/length(seqs_D_top)
      frac_IgD_also_IgM[pat,vis,"rest-repertoire"] = length(intersect(seqs_D_rest,seqs_M_rest))/length(seqs_D_rest)
    } else {
      top_lineage = lineage_abundances$LINEAGE[lineage_abundances$VISIT==vis][1]
      data = subset(data0,VISIT==vis)
      seqs_D_top = subset(data,LINEAGE==top_lineage & PRCONS=="IgD",SEQUENCE_VDJ,drop=T)
      seqs_M_top = subset(data,LINEAGE==top_lineage & PRCONS=="IgM",SEQUENCE_VDJ,drop=T)
      seqs_D_rest = subset(data,LINEAGE!=top_lineage & PRCONS=="IgD",SEQUENCE_VDJ,drop=T)
      seqs_M_rest = subset(data,LINEAGE!=top_lineage & PRCONS=="IgM",SEQUENCE_VDJ,drop=T)
      frac_IgD_also_IgM[pat,vis,"top-lineage"] = length(intersect(seqs_D_top,seqs_M_top))/length(seqs_D_top)
      frac_IgD_also_IgM[pat,vis,"rest-repertoire"] = length(intersect(seqs_D_rest,seqs_M_rest))/length(seqs_D_rest)
    }
  }
  
}
saveRDS(results,paste(fig_path,"Fig2results_corrected",sep=""))
saveRDS(frac_IgD_also_IgM,paste(fig_path,"frac_IgD_also_IgM_corrected",sep=""))
}

### Plot results
results = readRDS(paste(fig_path,"Fig2results_corrected",sep=""))
frac_IgD_also_IgM = readRDS(paste(fig_path,"frac_IgD_also_IgM_corrected",sep=""))
nmut_var = "NMUT_PER_SITE_GERMLINE"
oligoclonal_elderly = ##### vector of the relevant raw study participant IDs
nonoligoclonal_elderly = ##### vector of the relevant raw study participant IDs
young = c(pats_ellison_young_cmvp,pats_ellison_young_cmvn)

# Make Fig. 2B
if (save_eps) {
  pdf(file=paste(fig_path,"Fig2B.eps",sep=""), 
      width=4.25,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
visit = "V1"
pats = c(oligoclonal_elderly[order(id_to_dummy[oligoclonal_elderly])],
         nonoligoclonal_elderly[order(id_to_dummy[nonoligoclonal_elderly])],
         young[order(id_to_dummy[young])])
isotypes = c("IgE","IgA","IgG","IgD","IgM")
colors = brewer.pal(n=7,name="Dark2")[3:7]
M = as.matrix(apply(results[pats,visit,isotypes,"top-lineage","SEQCOUNT"],
                    1,function(s) s/sum(s)))
par(mar=c(4,6,4,2))
bp = barplot(M,
             ylab="",
             border=NA,las=2,
             ylim=c(0,1),xpd=F,
             space=c(1,rep(0.1,length(oligoclonal_elderly)-1),
                     1,rep(0.1,length(nonoligoclonal_elderly)-1),
                     1,rep(0.1,length(young)-1)),
             col=colors,
             names.arg=id_to_dummy[pats])
text(x=bp[c(3,8,16)],y=1,pos=3,xpd=T,
     labels=c("oligocl.\nelderly","non-oligocl.\nelderly","young"))
mtext(side=2,line=2.5,
      paste("Proportion of sequences\nin dominant lineage (",visit_to_day(visit),")",sep=""))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Fig2B.eps",sep=""), 
      width=3.5,height=5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0.5,0.5,legend=rev(isotypes),bty='n',pt.cex=2,
       pch=15,col=rev(colors),text.col=rev(colors))
if (save_eps) dev.off()


# Make Fig. 2C
if (save_eps) {
  pdf(file=paste(fig_path,"Figure2C.eps",sep=""), 
      width=3.4,height=3.4,paper="a4",
      colormodel="rgb",pointsize = 10)
}
visit = "V1"
x = (results[,visit,"IgD",,nmut_var]+results[,visit,"IgM",,nmut_var])/
  (results[,visit,"IgD",,"SEQCOUNT"]+results[,visit,"IgM",,"SEQCOUNT"])
y = (results[,visit,"IgA",,nmut_var]+results[,visit,"IgG",,nmut_var])/
  (results[,visit,"IgA",,"SEQCOUNT"]+results[,visit,"IgG",,"SEQCOUNT"])
plot(x[oligoclonal_elderly,"top-lineage"],y[oligoclonal_elderly,"top-lineage"],
     pch=17,col=oligo_color(oligoclonal_elderly[1]),
     xlim=range(x,na.rm=T),ylim=range(y,na.rm=T),
     xlab="",
     ylab="")
points(x[oligoclonal_elderly,"rest-repertoire"],y[oligoclonal_elderly,"rest-repertoire"],
       pch=2,col=oligo_color(oligoclonal_elderly[1]))
points(x[nonoligoclonal_elderly,"top-lineage"],y[nonoligoclonal_elderly,"top-lineage"],
       pch=15,col=oligo_color(nonoligoclonal_elderly[1]))
points(x[nonoligoclonal_elderly,"rest-repertoire"],y[nonoligoclonal_elderly,"rest-repertoire"],
       pch=0,col=oligo_color(nonoligoclonal_elderly[1]))
points(x[young,"top-lineage"],y[young,"top-lineage"],
       pch=16,col=oligo_color(young[1]))
points(x[young,"rest-repertoire"],y[young,"rest-repertoire"],
       pch=1,col=oligo_color(young[1]))
mtext(side=1,line=2,paste("Mutations per site in IgM & IgD (",visit_to_day(visit),")",sep=""))
mtext(side=2,line=2,paste("Mutations per site in IgA & IgG (",visit_to_day(visit),")",sep=""))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure2C_legend1.eps",sep=""), 
      width=10,height=3.6,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legendg(0.5,0.5,bty='n',
        c("oligocl. elderly","non-oligocl. elderly","young"),
        pch=list(c(17,2),c(15,0),c(16,1)),
        pt.space=1,
        text.col = c(oligo_color(oligoclonal_elderly[1]),
                     oligo_color(nonoligoclonal_elderly[1]),
                     oligo_color(young[1])),
        col=list(rep(oligo_color(oligoclonal_elderly[1]),2),
                 rep(oligo_color(nonoligoclonal_elderly[1]),2),
                 rep(oligo_color(young[1]),2)))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure2C_legend2.eps",sep=""), 
      width=10,height=3.6,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
leg = legendg(0.5,0.5,bty='n',
        "dominant lineage",
        pch=list(c(17,15,16)),
        pt.space=0.7,
        text.col = "black",
        col=list(rep('black',3)))
legendg(0.5,leg$text[['y']]-0.005,bty='n',
        "rest of repertoire",
        pch=list(c(2,0,1)),
        pt.space=2*0.7,
        text.col = "black",
        col=list(rep('black',3)))
if (save_eps) dev.off()

# Make Fig. 2D
if (save_eps) {
  pdf(file=paste(fig_path,"Figure2D.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
visit = "V1"
isotype = "IgD"
y = results[,visit,isotype,,nmut_var]/
  results[,visit,isotype,,"SEQCOUNT"]
data_list = list(y[oligoclonal_elderly,"top-lineage"],
                 y[oligoclonal_elderly,"rest-repertoire"],
                 y[nonoligoclonal_elderly,"top-lineage"],
                 y[nonoligoclonal_elderly,"rest-repertoire"],
                 y[young,"top-lineage"],
                 y[young,"rest-repertoire"])
delta1 = 0.5
delta2 = 1
start = 2
box_pos = start+c(0,delta1,delta1+delta2,2*delta1+delta2,2*delta1+2*delta2,3*delta1+2*delta2)
name_pos = 0.2+c(mean(box_pos[1],box_pos[2]),
                 mean(box_pos[3],box_pos[4]),
                 mean(box_pos[5],box_pos[6]))
colors = brewer.pal(n=6,name="Set1")[3:4]
par(mar=c(4,8,2,2))
boxplot(list(c(),c(),c()),
        at=name_pos,
        las=1,
        names = rep("",3),
        ylab = "",
        ylim=range(data_list,na.rm=T))
boxplot(data_list,add=T,
        border = rep(colors,3),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",6),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=3,
      paste("Mutations per site in ",isotype," (",visit_to_day(visit),")",sep=""))
mtext(side=1,line=1.5,at=name_pos,
      c("oligocl.\nelderly","non-oligocl.\nelderly","young\n"))
legend("topright",bty='n',c('dominant lineage','rest of repertoire'),
       text.col=colors)
p1 = wilcox.test(y[oligoclonal_elderly,"top-lineage"],
                 y[nonoligoclonal_elderly,"top-lineage"])$p.value
p2 = wilcox.test(y[oligoclonal_elderly,"top-lineage"],
                 y[young,"top-lineage"])$p.value
wilcox.test(y[oligoclonal_elderly,"top-lineage"],
            y[c(nonoligoclonal_elderly,young),"top-lineage"])$p.value
if (save_eps) dev.off()



# Make Fig. 2E
if (save_eps) {
  pdf(file=paste(fig_path,"Figure2E.eps",sep=""), 
      width=4,height=3,paper="a4",#bg="white",fonts="Helvetica",
      colormodel="rgb",pointsize = 10)
}
visit = "V1"
y = 100*(1-frac_IgD_also_IgM[,visit,])
data_list = list(y[oligoclonal_elderly,"top-lineage"],
                 y[oligoclonal_elderly,"rest-repertoire"],
                 y[nonoligoclonal_elderly,"top-lineage"],
                 y[nonoligoclonal_elderly,"rest-repertoire"],
                 y[young,"top-lineage"],
                 y[young,"rest-repertoire"])
delta1 = 0.5
delta2 = 1
start = 2
box_pos = start+c(0,delta1,delta1+delta2,2*delta1+delta2,2*delta1+2*delta2,3*delta1+2*delta2)
name_pos = 0.2+c(mean(box_pos[1],box_pos[2]),
                 mean(box_pos[3],box_pos[4]),
                 mean(box_pos[5],box_pos[6]))
colors = brewer.pal(n=6,name="Set1")[3:4]
par(mar=c(4,8,2,2))
boxplot(list(c(),c(),c()),
        at=name_pos,
        las=1,
        names = rep("",3),
        ylim=c(0,100),
        ylab = "",
        ylim=range(data_list,na.rm=T))
boxplot(data_list,add=T,
        border = rep(colors,3),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",6),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=3,
      paste("% IgD sequences\nabsent from IgM (",
            visit_to_day(visit),")",sep=""))
mtext(side=1,line=1.5,at=name_pos,
      c("oligocl.\nelderly","non-oligocl.\nelderly","young\n"))
legend("topright",bty='n',c('dominant lineage','rest of repertoire'),
       text.col=colors)
p1 = wilcox.test(y[oligoclonal_elderly,"top-lineage"],
                 y[nonoligoclonal_elderly,"top-lineage"])$p.value
p2 = wilcox.test(y[oligoclonal_elderly,"top-lineage"],
                 y[young,"top-lineage"])$p.value
wilcox.test(y[oligoclonal_elderly,"top-lineage"],
            y[c(nonoligoclonal_elderly,young),"top-lineage"])$p.value
if (save_eps) dev.off()




