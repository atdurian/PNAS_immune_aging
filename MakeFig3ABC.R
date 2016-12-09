library(beeswarm)
library(data.table)
library(alakazam)
library(RColorBrewer)
calculate = T
save_eps = T
fig_path = "/home/charles/Ellison_paper_figures/"

### Functions
visit_to_day = function(visit) {
  if (visit=="V1") return("day 0")
  if (visit=="V2") return("day 7")
  if (visit=="V3") return("day 28")
  if (visit=="all") return("all time points")
}

chao2006 = function(seqs1,seqs2,abunds1,abunds2) {
  shared_seqs = intersect(seqs1,seqs2)  
  names(abunds1) = seqs1
  names(abunds2) = seqs2
  shared_abunds1 = abunds1[shared_seqs]
  shared_abunds2 = abunds2[shared_seqs]
  D12 = length(shared_seqs)
  f11 = sum(shared_abunds1==1 & shared_abunds2==1,na.rm=T)
  f1p = sum(shared_abunds1==1,na.rm=T)
  fp1 = sum(shared_abunds2==1,na.rm=T)    
  f2p = sum(shared_abunds1==2,na.rm=T)
  fp2 = sum(shared_abunds2==2,na.rm=T)   
  if (f2p!=0 & fp2!=0) {
    S = D12+(f11/4)*(f1p/f2p)*(fp1/fp2)+f1p^2/(2*f2p)+fp1^2/(2*fp2)
  } else {
    S = D12+(f11/4)*(f1p/(f2p+1))*(fp1/(fp2+1))+f1p*(f1p-1)/(2*(f2p+1))+fp1*(fp1-1)/(2*(fp2+1))
  }
  return(S)
}

chao1 = function(abunds) {
  f1 = sum(abunds==1)
  f2 = sum(abunds==2)
  D = length(abunds) + f1^2/(2*f2)
  return(D)
}

sciNotation <- function(x, digits = 1) {
  if (length(x) > 1) {
    return(append(sciNotation(x[1]), sciNotation(x[-1])))
  }
  if (!x) return(0)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
    substitute(base %.% 10^exponent,
                   list(base = base, exponent = exponent))
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


### Compute sequence count estimates
path = "/media/charles/New_Volume/Ellison_analysis_output/"
visits = c("V1","V2","V3","all")
isotypes = c("IgD","IgM","IgA","IgG","IgE")

if (calculate) {
richnesses0 = array(dim=c(length(pats_ellison),length(visits),length(isotypes),length(c("DUPCOUNT","SEQCOUNT","CHAO1"))),
                     dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,variables=c("DUPCOUNT","SEQCOUNT","CHAO1")))
richnesses = array(dim=c(length(pats_ellison),length(visits),length(isotypes),2,length(c("DUPCOUNT","SEQCOUNT","CHAO1"))),
                dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,mutations=c("mutated","unmutated"),variables=c("DUPCOUNT","SEQCOUNT","CHAO1")))
shared_richnesses = array(dim=c(length(pats_ellison),length(visits),length(isotypes),length(isotypes),length(c("SEQCOUNT","CHAO"))),
                          dimnames=list(patients=pats_ellison,visits=visits,isotypes=isotypes,isotypes=isotypes,variables=c("SEQCOUNT","CHAO")))
for (pat in pats_ellison) {
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab",sep="")
  data0 = readChangeoDb(file)
  data0$SEQCOUNT = 1
  data0$MUTATED = (data0$NMUT_GERMLINE>0)
  
  # compute richnesses
  nmolec_data0 = aggregate(DUPCOUNT ~ VISIT+PRCONS, data=data0, FUN=sum) 
  nmolec_data_allvisits0 = aggregate(DUPCOUNT ~ PRCONS, data=data0, FUN=sum) 
  nseq_data0 = aggregate(SEQCOUNT ~ VISIT+PRCONS, data=data0, FUN=sum) 
  nseq_data_allvisits0 = aggregate(SEQCOUNT ~ PRCONS, data=data0, FUN=sum) 
  chao_data0 = aggregate(DUPCOUNT ~ VISIT+PRCONS, data=data0, FUN=chao1) 
  chao_data_allvisits0 = aggregate(DUPCOUNT ~ PRCONS, data=data0, FUN=chao1) 
  
  nmolec_data = aggregate(DUPCOUNT ~ VISIT+PRCONS+MUTATED, data=data0, FUN=sum) 
  nmolec_data_allvisits = aggregate(DUPCOUNT ~ PRCONS+MUTATED, data=data0, FUN=sum) 
  nseq_data = aggregate(SEQCOUNT ~ VISIT+PRCONS+MUTATED, data=data0, FUN=sum) 
  nseq_data_allvisits = aggregate(SEQCOUNT ~ PRCONS+MUTATED, data=data0, FUN=sum) 
  chao_data = aggregate(DUPCOUNT ~ VISIT+PRCONS+MUTATED, data=data0, FUN=chao1) 
  chao_data_allvisits = aggregate(DUPCOUNT ~ PRCONS+MUTATED, data=data0, FUN=chao1) 

  for (vis in visits) {
    for (iso in isotypes) {
        if (vis=="all") {
          richnesses0[pat,vis,iso,"DUPCOUNT"] = subset(nmolec_data_allvisits0,PRCONS==iso,DUPCOUNT,drop=T)
          richnesses0[pat,vis,iso,"SEQCOUNT"] = subset(nseq_data_allvisits0,PRCONS==iso,SEQCOUNT,drop=T)
          richnesses0[pat,vis,iso,"CHAO1"] = subset(chao_data_allvisits0,PRCONS==iso,DUPCOUNT,drop=T)

          richnesses[pat,vis,iso,"mutated","DUPCOUNT"] = subset(nmolec_data_allvisits,PRCONS==iso & MUTATED==T,DUPCOUNT,drop=T)
          richnesses[pat,vis,iso,"unmutated","DUPCOUNT"] = subset(nmolec_data_allvisits,PRCONS==iso & MUTATED==F,DUPCOUNT,drop=T)
          richnesses[pat,vis,iso,"mutated","SEQCOUNT"] = subset(nseq_data_allvisits,PRCONS==iso & MUTATED==T,SEQCOUNT,drop=T)
          richnesses[pat,vis,iso,"unmutated","SEQCOUNT"] = subset(nseq_data_allvisits,PRCONS==iso & MUTATED==F,SEQCOUNT,drop=T)
          richnesses[pat,vis,iso,"mutated","CHAO1"] = subset(chao_data_allvisits,PRCONS==iso & MUTATED==T,DUPCOUNT,drop=T)
          richnesses[pat,vis,iso,"unmutated","CHAO1"] = subset(chao_data_allvisits,PRCONS==iso & MUTATED==F,DUPCOUNT,drop=T)
        } else {
          richnesses0[pat,vis,iso,"DUPCOUNT"] = subset(nmolec_data0,VISIT==vis & PRCONS==iso,DUPCOUNT,drop=T)
          richnesses0[pat,vis,iso,"SEQCOUNT"] = subset(nseq_data0,VISIT==vis & PRCONS==iso,SEQCOUNT,drop=T)
          richnesses0[pat,vis,iso,"CHAO1"] = subset(chao_data0,VISIT==vis & PRCONS==iso,DUPCOUNT,drop=T)

          richnesses[pat,vis,iso,"mutated","DUPCOUNT"] = subset(nmolec_data,VISIT==vis & PRCONS==iso & MUTATED==T,DUPCOUNT,drop=T)
          richnesses[pat,vis,iso,"unmutated","DUPCOUNT"] = subset(nmolec_data,VISIT==vis & PRCONS==iso & MUTATED==F,DUPCOUNT,drop=T)
          richnesses[pat,vis,iso,"mutated","SEQCOUNT"] = subset(nseq_data,VISIT==vis & PRCONS==iso & MUTATED==T,SEQCOUNT,drop=T)
          richnesses[pat,vis,iso,"unmutated","SEQCOUNT"] = subset(nseq_data,VISIT==vis & PRCONS==iso & MUTATED==F,SEQCOUNT,drop=T)
          richnesses[pat,vis,iso,"mutated","CHAO1"] = subset(chao_data,VISIT==vis & PRCONS==iso & MUTATED==T,DUPCOUNT,drop=T)
          richnesses[pat,vis,iso,"unmutated","CHAO1"] = subset(chao_data,VISIT==vis & PRCONS==iso & MUTATED==F,DUPCOUNT,drop=T)
        }
      }
    }
  
  # compute shared richnesses
  for (vis in visits) {
    for (iso1 in isotypes) {
      for (iso2 in isotypes) {
        if (vis=="all") {
          seqs_iso1 = subset(data0,PRCONS==iso1,SEQUENCE_VDJ,drop=T)
          abunds_iso1 = subset(data0,PRCONS==iso1,DUPCOUNT,drop=T)
          seqs_iso2 = subset(data0,PRCONS==iso2,SEQUENCE_VDJ,drop=T)
          abunds_iso2 = subset(data0,PRCONS==iso2,DUPCOUNT,drop=T)
          shared_richnesses[pat,vis,iso1,iso2,"SEQCOUNT"] = length(intersect(seqs_iso1,seqs_iso2)) 
          shared_richnesses[pat,vis,iso1,iso2,"CHAO"] = chao2006(seqs_iso1,seqs_iso2,abunds_iso1,abunds_iso2) 
        } else {
          data = subset(data0,VISIT==vis)
          seqs_iso1 = subset(data,PRCONS==iso1,SEQUENCE_VDJ,drop=T)
          abunds_iso1 = subset(data,PRCONS==iso1,DUPCOUNT,drop=T)
          seqs_iso2 = subset(data,PRCONS==iso2,SEQUENCE_VDJ,drop=T)
          abunds_iso2 = subset(data,PRCONS==iso2,DUPCOUNT,drop=T)
          shared_richnesses[pat,vis,iso1,iso2,"SEQCOUNT"] = length(intersect(seqs_iso1,seqs_iso2)) 
          shared_richnesses[pat,vis,iso1,iso2,"CHAO"] = chao2006(seqs_iso1,seqs_iso2,abunds_iso1,abunds_iso2) 
        }
      }
    }
  }
  
}
saveRDS(richnesses0,paste(fig_path,"richnesses0_corrected",sep=""))
saveRDS(richnesses,paste(fig_path,"richnesses_corrected",sep=""))
saveRDS(shared_richnesses,paste(fig_path,"shared_richnesses_corrected",sep=""))
}



### Plot results
richnesses0 = readRDS(paste(fig_path,"richnesses0_corrected",sep=""))
richnesses = readRDS(paste(fig_path,"richnesses_corrected",sep=""))
shared_richnesses = readRDS(paste(fig_path,"shared_richnesses_corrected",sep=""))


# Make Fig. 3A, upper panel
if (save_eps) {
  pdf(file=paste(fig_path,"Fig3A_upper.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
par(mar=c(4,6,2,2))
pats = c(pats_ellison_young_cmvn[order(id_to_dummy[pats_ellison_young_cmvn])],
         pats_ellison_young_cmvp[order(id_to_dummy[pats_ellison_young_cmvp])],
         pats_ellison_old_cmvn[order(id_to_dummy[pats_ellison_old_cmvn])],
         pats_ellison_old_cmvp[order(id_to_dummy[pats_ellison_old_cmvp])])
visit = "V1"
isotypes = c("IgE","IgA","IgG","IgD","IgM")
colors = brewer.pal(n=7,name="Dark2")[3:7]
M = t(as.matrix(richnesses0[pats,visit,isotypes,"CHAO1"]))
bp = barplot(M,
             axes=F,ann=F,
             border=NA,las=2,xpd=F,
             space=c(1,rep(0.1,length(pats_ellison_young_cmvn)-1),
                     1,rep(0.1,length(pats_ellison_young_cmvp)-1),
                     1,rep(0.1,length(pats_ellison_old_cmvn)-1),
                     1,rep(0.1,length(pats_ellison_old_cmvp)-1)),
             col=colors,
             names.arg=rep("",length(pats)))
axis(side=2,at=axTicks(2),las=2,label=as.expression(sciNotation(axTicks(2), 1)))
mtext(side=2,line=3.5,paste("Sequence diversity (",visit_to_day(visit),")",sep=""))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Fig3A_upper_legend.eps",sep=""), 
      width=3.5,height=5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0.5,0.5,legend=rev(isotypes),bty='n',pt.cex=2,
       pch=15,col=rev(colors),text.col=rev(colors)
       )
if (save_eps) dev.off()


# Make Fig. 3A, lower panel
if (save_eps) {
  pdf(file=paste(fig_path,"Fig3A_lower.eps",sep=""), 
      width=4,height=2.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
par(mar=c(4,6,2,2))
pats = c(pats_ellison_young_cmvn[order(id_to_dummy[pats_ellison_young_cmvn])],
         pats_ellison_young_cmvp[order(id_to_dummy[pats_ellison_young_cmvp])],
         pats_ellison_old_cmvn[order(id_to_dummy[pats_ellison_old_cmvn])],
         pats_ellison_old_cmvp[order(id_to_dummy[pats_ellison_old_cmvp])])
visit = "V1"
isotypes = c("IgE","IgA","IgG","IgD","IgM")
colors = brewer.pal(n=7,name="Dark2")[3:7]
M = t(as.matrix(richnesses0[pats,visit,isotypes,"SEQCOUNT"]))
M = apply(M,2,function(s) s/sum(s))
bp = barplot(M,
             ylab=paste("Proportion of sequences (",visit_to_day(visit),")",sep=""),
             border=NA,las=2,xpd=F,
             space=c(1,rep(0.1,length(pats_ellison_young_cmvn)-1),
                     1,rep(0.1,length(pats_ellison_young_cmvp)-1),
                     1,rep(0.1,length(pats_ellison_old_cmvn)-1),
                     1,rep(0.1,length(pats_ellison_old_cmvp)-1)),
             col=colors,
             names.arg=id_to_dummy[pats],cex.names=0.8)
if (save_eps) dev.off()




# Make Fig. 3C 
if (save_eps) {
  pdf(file=paste(fig_path,"Figure3C.eps",sep=""), 
      width=3.5,height=2.55,paper="a4",
      colormodel="rgb",pointsize = 10)
}
vis = "V1"
iso1 = "IgD"
iso2 = "IgM"
y = 100*shared_richnesses[pats_ellison,vis,iso1,iso2,"SEQCOUNT"]/richnesses0[pats_ellison,vis,iso2,"SEQCOUNT"]
par(mar=c(4,6,2,2))
data_list = list(y[pats_ellison_young_cmvn],
                 y[pats_ellison_young_cmvp],
                 y[pats_ellison_old_cmvn],
                 y[pats_ellison_old_cmvp])
color=brewer.pal(n=6,name="Set1")[3]
boxplot(data_list,
        names=rep("",4),
        border=color,
        ylab="",boxwex=0.7,
        outline=F,range=Inf)
beeswarm(data_list,"swarm",add=T,pch=20,cex=1.25)
mtext(side=2,line=3, "% IgM sequences shared with IgD")
mtext(side=1,line=1.5,at=1:4,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
wilcox.test(y[c(pats_ellison_young_cmvn,pats_ellison_old_cmvn)],
            y[c(pats_ellison_young_cmvp,pats_ellison_old_cmvp)])
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure3C_legend.eps",sep=""), 
      width=3.5,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("top","day 0",text.col=color,bty='n')
if (save_eps) dev.off()


### Make Fig. 3B, middle panel
if (save_eps) {
  pdf(file=paste(fig_path,"Figure3B_middle.eps",sep=""), 
      width=3.75,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
iso = "IgD"
mut_status = "unmutated"
var = "CHAO1"
data_list = list(richnesses[pats_ellison_young_cmvn,"V1",iso,mut_status,var],
                 richnesses[pats_ellison_young_cmvn,"V2",iso,mut_status,var],
                 richnesses[pats_ellison_young_cmvn,"V3",iso,mut_status,var],
                 richnesses[pats_ellison_young_cmvp,"V1",iso,mut_status,var],
                 richnesses[pats_ellison_young_cmvp,"V2",iso,mut_status,var],
                 richnesses[pats_ellison_young_cmvp,"V3",iso,mut_status,var],
                 richnesses[pats_ellison_old_cmvn,"V1",iso,mut_status,var],
                 richnesses[pats_ellison_old_cmvn,"V2",iso,mut_status,var],
                 richnesses[pats_ellison_old_cmvn,"V3",iso,mut_status,var],
                 richnesses[pats_ellison_old_cmvp,"V1",iso,mut_status,var],
                 richnesses[pats_ellison_old_cmvp,"V2",iso,mut_status,var],
                 richnesses[pats_ellison_old_cmvp,"V3",iso,mut_status,var])
wilcox.test(unlist(data_list[c(1,4)]),unlist(data_list[c(7,10)]))
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
        ylab = "",
        axes=F,
        ylim=range(data_list))
boxplot(data_list,add=T,
        border = rep(colors,4),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",8),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=4, 'Naïve diversity')
axis(side=1,at=name_pos,label=F)
axis(side=2,at=axTicks(2),las=2,label=as.expression(sciNotation(axTicks(2), 1)))
box()
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure3B_legend.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("top",c("day 0","day 7","day 28"),text.col=colors,bty='n')
if (save_eps) dev.off()


### Make Fig. 3B, upper panel
if (save_eps) {
  pdf(file=paste(fig_path,"Figure3B_upper.eps",sep=""), 
      width=3.75,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
iso = "IgD"
mut_status = "unmutated"
var = "CHAO1"
y = 100*richnesses[,,iso,mut_status,var]
total = apply(richnesses[,,,,var],c("patients","visits"),sum)
data_list = list(y[pats_ellison_young_cmvn,"V1"]/total[pats_ellison_young_cmvn,"V1"],
                 y[pats_ellison_young_cmvn,"V2"]/total[pats_ellison_young_cmvn,"V2"],
                 y[pats_ellison_young_cmvn,"V3"]/total[pats_ellison_young_cmvn,"V3"],
                 y[pats_ellison_young_cmvp,"V1"]/total[pats_ellison_young_cmvp,"V1"],
                 y[pats_ellison_young_cmvp,"V2"]/total[pats_ellison_young_cmvp,"V2"],
                 y[pats_ellison_young_cmvp,"V3"]/total[pats_ellison_young_cmvp,"V3"],
                 y[pats_ellison_old_cmvn,"V1"]/total[pats_ellison_old_cmvn,"V1"],
                 y[pats_ellison_old_cmvn,"V2"]/total[pats_ellison_old_cmvn,"V2"],
                 y[pats_ellison_old_cmvn,"V3"]/total[pats_ellison_old_cmvn,"V3"],
                 y[pats_ellison_old_cmvp,"V1"]/total[pats_ellison_old_cmvp,"V1"],
                 y[pats_ellison_old_cmvp,"V2"]/total[pats_ellison_old_cmvp,"V2"],
                 y[pats_ellison_old_cmvp,"V3"]/total[pats_ellison_old_cmvp,"V3"])
wilcox.test(unlist(data_list[c(1,4)]),unlist(data_list[c(7,10)]))
wilcox.test(unlist(data_list[c(2,5)]),unlist(data_list[c(8,11)]))
wilcox.test(unlist(data_list[c(3,6)]),unlist(data_list[c(9,12)]))
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
mtext(side=2,line=3, '% naïve sequences')
axis(side=1,at=name_pos,label=F)
axis(side=2,at=axTicks(2),las=2,label=axTicks(2))
box()
if (save_eps) dev.off()


### Make Fig. 3B, lower panel
if (save_eps) {
  pdf(file=paste(fig_path,"Figure3B_lower.eps",sep=""), 
      width=3.75,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
isos = c("IgA","IgG")
var = "CHAO1"
data_list = list(apply(richnesses0[pats_ellison_young_cmvn,"V1",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_young_cmvn,"V2",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_young_cmvn,"V3",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_young_cmvp,"V1",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_young_cmvp,"V2",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_young_cmvp,"V3",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_old_cmvn,"V1",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_old_cmvn,"V2",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_old_cmvn,"V3",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_old_cmvp,"V1",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_old_cmvp,"V2",isos,var],"patients",sum),
                 apply(richnesses0[pats_ellison_old_cmvp,"V3",isos,var],"patients",sum))
wilcox.test(unlist(data_list[c(1,4)]),unlist(data_list[c(7,10)]))
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
        ylab = "",
        axes=F,
        ylim=range(data_list))
boxplot(data_list,add=T,
        border = rep(colors,4),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",8),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=4, 'Antigen-experienced diversity')
mtext(side=1,line=1.5,at=name_pos,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
axis(side=1,at=name_pos,label=F)
axis(side=2,at=axTicks(2),las=2,label=as.expression(sciNotation(axTicks(2), 1)))
box()
if (save_eps) dev.off()



