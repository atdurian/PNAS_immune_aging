library(ape)
library(pheatmap)
library(RColorBrewer)
library(beeswarm)
save_eps = T


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

visit_to_day = paste("day",c(0,7,28),sep=" ")
names(visit_to_day) = c("V1","V2","V3")

### Read UniFrac data
path = "/media/charles/New_Volume/Ellison_analysis_output/"
filename = "RUBL.seed-1.Nseq-10000.GERMLINE_IMGT.tab"
D = read.table(paste(path,filename,sep=""),
               header=T,check.names=F)
D = as.matrix(D)
D0 = D

D = D[grepl("V1",rownames(D)),grepl("V1",colnames(D))]
rownames(D) = id_to_dummy[as.character(sapply(rownames(D), function(s) strsplit(s,"_")[[1]][1]))] # (rownames/colnames were in the form "participantID_visitID")
colnames(D) = id_to_dummy[as.character(sapply(colnames(D), function(s) strsplit(s,"_")[[1]][1]))]
D = D[id_to_dummy,id_to_dummy]


### Make Fig. 1B
annot = data.frame(group=pat_data$group,row.names=pat_data$DummyID)
min_nonzero = min(D[row(D)!=col(D)])
delta = 0.06
diag_value = min_nonzero-delta
D_hm = D; diag(D_hm) = diag_value

if (save_eps) {
  pdf(file=paste(fig_path,"Figure1B.eps",sep=""), 
      width=6,height=6,
      paper="a4",onefile=F,
      colormodel="rgb",pointsize = 12)
}
pheatmap(D_hm, 
         cluster_rows=F, cluster_cols=F, 
         color = c('grey',colorRampPalette(rev(brewer.pal(n=11,name ="RdBu")))(50)),
         breaks = c(diag_value,seq(min_nonzero-1e-10,0.96,length.out=51)),
         legend_breaks = c(diag_value+delta/2,seq(ceiling(min_nonzero/0.02)*0.02,0.96,0.04)),
         legend_labels = c(0,seq(ceiling(min_nonzero/0.02)*0.02,0.96,0.04)),
         gaps_row=5*1:3,gaps_col=5*1:3,
         border_color=NA,
         cellwidth=12,cellheight=12,
         fontsize=12,fontsize_number=12)
if (save_eps) dev.off()

### Compute age-sum matrix
age_sum_matrix = matrix(nrow=nrow(D),ncol=ncol(D),dimnames=list(rownames(D),colnames(D)))
for (s1 in rownames(age_sum_matrix)) {
  for (s2 in colnames(age_sum_matrix)) {
    if (s1==s2) {
      age_sum_matrix[s1,s2] = 0
    } else {
      age_sum_matrix[s1,s2] = ages_ellison[s1]+ages_ellison[s2]
    }
  }
}

### Test for correlation between RUBL and age-sum matrix
set.seed(10)
mantel.test(D,age_sum_matrix,nperm=5000,graph=F,alternative="two.sided")
cor.test(D[upper.tri(D)],age_sum_matrix[upper.tri(age_sum_matrix)])

cmvn_cond = grepl("N",rownames(D))
D_cmvn = D[cmvn_cond,cmvn_cond]
age_sum_cmvn = age_sum_matrix[cmvn_cond,cmvn_cond]
cor.test(D_cmvn[upper.tri(D_cmvn)],age_sum_cmvn[upper.tri(age_sum_cmvn)])

cmvp_cond = grepl("P",rownames(D))
D_cmvp = D[cmvp_cond,cmvp_cond]
age_sum_cmvp = age_sum_matrix[cmvp_cond,cmvp_cond]
cor.test(D_cmvp[upper.tri(D_cmvp)],age_sum_cmvp[upper.tri(age_sum_cmvp)])


### Make Fig. 1C
if (save_eps) {
  pdf(file=paste(fig_path,"Figure1C.eps",sep=""), 
      width=3.5,height=4,paper="a4",
      colormodel="rgb",pointsize = 10)
}
DYY_N = D[grepl("YN",rownames(D)),grepl("YN",colnames(D))]; DYY_N=DYY_N[upper.tri(DYY_N)]
DEE_N = D[grepl("EN",rownames(D)),grepl("EN",colnames(D))]; DEE_N=DEE_N[upper.tri(DEE_N)]
DYE_N = c(D[grepl("YN",rownames(D)),grepl("EN",colnames(D))])

DYY_P = D[grepl("YP",rownames(D)),grepl("YP",colnames(D))]; DYY_P=DYY_P[upper.tri(DYY_P)]
DEE_P = D[grepl("EP",rownames(D)),grepl("EP",colnames(D))]; DEE_P=DEE_P[upper.tri(DEE_P)]
DYE_P = c(D[grepl("YP",rownames(D)),grepl("EP",colnames(D))])

data_list = list(DYY_N,DYY_P,
                 DYE_N,DYE_P,
                 DEE_N,DEE_P)
par(mar=c(8,5,2,2))
delta1 = 0.5
delta2 = 1
start = 2
box_pos = start+delta2*rep(0:2,each=2)+delta1*c(0:1,1:2,2:3)
name_pos = 0.2+c(mean(box_pos[1],box_pos[2]),
                 mean(box_pos[3],box_pos[4]),
                 mean(box_pos[5],box_pos[6]))
colors = brewer.pal(n=6,name="Set1")[3:4]
boxplot(list(c(),c(),c()),
        at=name_pos,
        las=1,
        names = rep("",3),
        ylab = "",
        ylim=range(data_list))
boxplot(data_list,add=T,
        border = rep(colors,3),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",6),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=3, 'UniFrac(participant A, particpant B) at day 0')

wilcox.test(DYY_N,DYY_P)

YY_age_sum = age_sum_matrix[grepl("Y",rownames(age_sum_matrix)),grepl("Y",colnames(age_sum_matrix))]
YY_age_sum_range = range(YY_age_sum[upper.tri(YY_age_sum)])
YE_age_sum = age_sum_matrix[grepl("Y",rownames(age_sum_matrix)),grepl("E",colnames(age_sum_matrix))]
YE_age_sum_range = range(YE_age_sum[upper.tri(YE_age_sum)])
EE_age_sum = age_sum_matrix[grepl("E",rownames(age_sum_matrix)),grepl("E",colnames(age_sum_matrix))]
EE_age_sum_range = range(EE_age_sum[upper.tri(EE_age_sum)])
mtext(side=1,line=2,at=name_pos,
      c("A: young\nB: young",
        "A: young\nB: elderly",
        "A: elderly\nB: elderly"))
mtext(side=1,line=4,at=c(0.85,name_pos),xpd=T,
      c("age(A)\n+age(B)\nin years:", paste(YY_age_sum_range[1]," to ",YY_age_sum_range[2],sep=""),
        paste(YE_age_sum_range[1]," to ",YE_age_sum_range[2],sep=""),
        paste(EE_age_sum_range[1]," to ",EE_age_sum_range[2],sep="")))
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"Figure1C_legend.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0.5,0.5,bty='n',c('CMV-','CMV+'),
       text.col=colors)
if (save_eps) dev.off()



