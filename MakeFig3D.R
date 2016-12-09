library(beeswarm)
library(data.table)
library(alakazam)
library(RColorBrewer)
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

sequence_entropy = function(seqs) {
  seq_list = lapply(unique(seqs), function(s) strsplit(s,"")[[1]])
  M = matrix(unlist(seq_list), ncol = length(seq_list[[1]]), byrow = TRUE)
  entropy_by_nucleotide = apply(M,2,function(s) {
    t = table(s)
    p = t/sum(t)
    return(-sum(p*log(p)))
  })
  mean_entropy_per_nucleotide = mean(entropy_by_nucleotide)
  return(mean_entropy_per_nucleotide)
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


### Compute within-lineage entropies:
path = "/media/charles/New_Volume/Ellison_analysis_output/"
lineage_field = "LINEAGE_dissim0.1"
visits = c("V1","V2","V3","all")
if (calculate) {
lineage_entropies = lineage_sizes = list()
for (pat in pats_ellison) {  
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab",sep="")
  data0 = readChangeoDb(file)
  setnames(data0,old=lineage_field,new="LINEAGE")
  data0$INCIDENCE = 1
  count_data = aggregate(cbind(DUPCOUNT,INCIDENCE) ~ LINEAGE+VISIT,data0,sum)
  count_data_allvisits = aggregate(cbind(DUPCOUNT,INCIDENCE) ~ LINEAGE,data0,sum)
  for (vis in visits) {
    if (vis=="all") {
      pat_vis = paste(pat,vis,sep=",")
      lineages = subset(count_data_allvisits,INCIDENCE>=10,LINEAGE,drop=T)
      data = subset(data0,LINEAGE %in% lineages)
      lineage_entropies[[pat_vis]] = aggregate(JUNCTION ~ LINEAGE,data,sequence_entropy)$JUNCTION
      lineage_sizes[[pat_vis]] = aggregate(JUNCTION ~ LINEAGE,data,function(seqs) length(unique(seqs)))$JUNCTION
    } else {
      pat_vis = paste(pat,vis,sep=",")
      lineages = subset(count_data,VISIT==vis & INCIDENCE>=10,LINEAGE,drop=T)
      data = subset(data0,VISIT==vis & LINEAGE %in% lineages)
      lineage_entropies[[pat_vis]] = aggregate(JUNCTION ~ LINEAGE,data,sequence_entropy)$JUNCTION
      lineage_sizes[[pat_vis]] = aggregate(JUNCTION ~ LINEAGE,data,function(seqs) length(unique(seqs)))$JUNCTION
    }
  }
}
saveRDS(lineage_entropies,paste(fig_path,"lineage_entropies_corrected",sep=""))
saveRDS(lineage_sizes,paste(fig_path,"lineage_sizes_corrected",sep=""))
}

### Plotting
lineage_entropies = readRDS(paste(fig_path,"lineage_entropies_corrected",sep=""))
lineage_sizes = readRDS(paste(fig_path,"lineage_sizes_corrected",sep=""))

# Make Fig. 3D
if (save_eps) {
  pdf(file=paste(fig_path,"Figure3D.eps",sep=""), 
      width=3.5,height=3.1,paper="a4",
      colormodel="rgb",pointsize = 10)
}
ma = function(x,n=5){stats::filter(x,rep(1/n,n), sides=1)}
yrange=range(unlist(lineage_entropies))
xrange=range(unlist(lineage_sizes))
visit = "V1"
for (i in 1:length(pats_ellison)) {
  pat = pats_ellison[i]
  pat_vis = paste(pat,visit,sep=",")
  x = lineage_sizes[[pat_vis]]
  y = lineage_entropies[[pat_vis]]
  df = data.frame(x,y)
  df = aggregate(y ~ x,df,mean)
  if (i==1) {
    plot(df$x,ma(df$y),
         ann=F,axes=F,xaxt="n",yaxt="n",
         xlim=c(10,2000),
         ylim=c(0,0.8),
         type='l',col=group_color(pat),log='x',
         xlab="Number of unique CDR3 sequences",
         ylab="Average entropy per nucleotide in CDR3")
  } else {
    lines(df$x,ma(df$y),
          type='l',col=group_color(pat))
  }
}
mtext(side=1,line=1.5,"Unique CDR3 sequences in lineage")
mtext(side=2,line=1.75,"Entropy per nucleotide in CDR3")
axis(side=1,at=axTicks(1),mgp=c(3, .5, 0))
axis(side=2,at=axTicks(2),mgp=c(3, .5, 0))
box()

if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure3D_legend.eps",sep=""), 
      width=3.5,height=3.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("topleft",c("young, CMV-","young, CMV+","elderly, CMV-","elderly, CMV+"),
       text.col=sapply(c(pats_ellison_young_cmvn[1],
                         pats_ellison_young_cmvp[1],
                         pats_ellison_old_cmvn[1],
                         pats_ellison_old_cmvp[1]),
                       group_color),
       bty='n')
if (save_eps) dev.off()


