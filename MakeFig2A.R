library(alakazam)
library(plyr)
library(reshape2) 
library(data.table)
library(RColorBrewer)

calculate = T
save_eps = T
fig_path = "/home/charles/Ellison_paper_figures/"

### Functions
make_color_palette = function(palette,n_selected_lins) {
  max_n_colors = brewer.pal.info[palette,"maxcolors"]
  n_repeats = n_selected_lins%/%max_n_colors
  n_remainder = n_selected_lins%%max_n_colors
  all_colors = brewer.pal(n=max_n_colors,name=palette)
  result = c("grey",
             rep(all_colors,n_repeats),
             all_colors[1:n_remainder])
  return(result)
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


### Calculate values for the figure
path="/media/charles/New_Volume/Ellison_analysis_output/"
N_lin_select = 20
lineage_field = "LINEAGE_dissim0.1"
all_pats_dupcounts_plot = list()
oligo_character = list()

if (calculate) {
for (type in unique(pat_data$group)) {
  patients_dummy = subset(pat_data,group==type,DummyID,drop=T)
  all_pats_dupcounts_plot[[type]] = matrix(nrow=N_lin_select+1,ncol=0)
  for (pat_dummy in patients_dummy) {  
    pat = dummy_to_id[pat_dummy]
    file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass.tab",sep="")
    data = readChangeoDb(file)
    setnames(data,old=lineage_field,new="LINEAGE")
    
    dupcounts = dcast(data, LINEAGE ~ VISIT, value.var="DUPCOUNT",fun.aggregate=sum)
    dupcounts$allVisits = dupcounts$V1+dupcounts$V2+dupcounts$V3
    dupcounts$V1_norm = dupcounts$V1/sum(dupcounts$V1)
    dupcounts$V2_norm = dupcounts$V2/sum(dupcounts$V2)
    dupcounts$V3_norm = dupcounts$V3/sum(dupcounts$V3)
    dupcounts$allVisits_norm = dupcounts$allVisits/sum(dupcounts$allVisits)
    rownames(dupcounts) = dupcounts$LINEAGE
    exp_lins = dupcounts$LINEAGE[order(dupcounts$allVisits,decreasing=T)[1:N_lin_select]]
    dupcounts = dupcounts[order(dupcounts$allVisits,decreasing=T),]
    dupcounts_plot = subset(dupcounts,LINEAGE %in% exp_lins,
                            c("V1_norm","V2_norm","V3_norm"))
    dupcounts_plot["rest",] = 1-colSums(dupcounts_plot)
    dupcounts_plot = as.matrix(dupcounts_plot)
    dupcounts_plot = dupcounts_plot[nrow(dupcounts_plot):1,]
    colnames(dupcounts_plot) = paste(colnames(dupcounts_plot),pat_dummy,sep="_")
    all_pats_dupcounts_plot[[type]] = cbind(all_pats_dupcounts_plot[[type]],dupcounts_plot)
  }
  bp = barplot(all_pats_dupcounts_plot[[type]],
               ylab="proportion of molecules",
               border=NA,las=2,
               ylim=c(0,1),xpd=F,
               space=rep(c(1,rep(0.1,2)),length(patients_dummy)),
               col=make_color_palette("Paired",N_lin_select),
               names.arg=rep(c("day 0","day 7","day 28"),length(patients_dummy)))
  par(xpd=T)
  oligo_character[[type]] = sapply(patients_dummy,function(s) {
    M = all_pats_dupcounts_plot[[type]]
    columns = paste("V",1:3,"_norm_",s,sep="")
    exp_M = M[rownames(M)!="rest",columns]
    persistent_expanded = apply(exp_M,1,function(s) all(s>0.05))
    result = ifelse(sum(persistent_expanded)>=1,"oligocl.","")
    return(result)
  })
  text(x=bp[3*(1:length(patients_dummy))-1],y=1,pos=3,xpd=T,
       labels=paste(patients_dummy,oligo_character[[type]][patients_dummy],sep="\n"),
       col=ifelse(oligo_character[[type]][patients_dummy]=="oligocl.","red","black"))
  text(x=bp[length(bp)]+(bp[2]-bp[1])/2,y=0.5,pos=4,xpd=T,
       labels=type) 
}  
saveRDS(all_pats_dupcounts_plot,"/home/charles/R_objects/Fig1CD_all_pats_dupcounts_corrected_plot.RDS")
}

### Make Fig. 2A
all_pats_dupcounts_plot = readRDS("/home/charles/R_objects/Fig1CD_all_pats_dupcounts_corrected_plot.RDS")
colors = make_color_palette("Paired",N_lin_select)
for (type in unique(pat_data$group)) {
  patients_dummy = subset(pat_data,group==type,DummyID,drop=T)
  
  if (save_eps == TRUE) {
    pdf(file=paste(fig_path,"Figure1CD_",substr(patients_dummy[1],1,2),".eps",sep=""), 
        width=3.5,height=2.75,paper="a4",
        colormodel="rgb",pointsize = 10)
  }
  bp = barplot(all_pats_dupcounts_plot[[type]],
               ylab="",
               border=NA,las=2,
               ylim=c(0,1),xpd=F,
               space=rep(c(1,rep(0.1,2)),length(patients_dummy)),
               names.arg=rep("",3*length(patients_dummy)),
               col=colors)
  par(xpd=T)
  oligo_character[[type]] = sapply(patients_dummy,function(s) {
    M = all_pats_dupcounts_plot[[type]]
    columns = paste("V",1:3,"_norm_",s,sep="")
    exp_M = M[rownames(M)!="rest",columns]
    persistent_expanded = apply(exp_M,1,function(s) all(s>0.05))
    result = ifelse(sum(persistent_expanded)>=1,"oligocl.","")
    return(result)
  })
  text(x=bp[3*(1:length(patients_dummy))-1],y=1,pos=3,xpd=T,
       labels=paste(oligo_character[[type]][patients_dummy],patients_dummy,sep="\n")
       )
  mtext(side=1,line=0.5,at=c(-0.5,bp),
        c('Day: ',rep(c("0,","7,","28"),length(patients_dummy))))
  mtext(side=2,line=2.5,"Proportion of molecules")
  if (save_eps == TRUE) dev.off()
}

if (save_eps) {
  pdf(file=paste(fig_path,"Figure1CD_legend.eps",sep=""), 
      width=3,height=10,paper="a4",
      colormodel="rgb",pointsize = 8)
}
plot.new()
legend(0.5,0.5,legend=c(1:N_lin_select,"lower"),pch=15,col=rev(colors),ncol=1,
       xpd=T,bty='n',pt.cex=2)
text(0.55,0.52,"Lineages,\nranked by\nabundance:")
if (save_eps) dev.off()




