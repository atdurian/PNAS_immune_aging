library(ape)
library(beeswarm)
save_eps = T

set.seed(20)

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

### Read UniFrac data
path = "/media/charles/New_Volume/Ellison_analysis_output/"
filename = "RUBL.seed-1.Nseq-10000.GERMLINE_IMGT.tab"
D = read.table(paste(path,filename,sep=""),
               header=T,check.names=F)
D = as.matrix(D)

D_V1V3 = D[grepl("V1",rownames(D)),grepl("V3",colnames(D))]
rownames(D_V1V3) = id_to_dummy[as.character(sapply(rownames(D_V1V3), function(s) strsplit(s,"_")[[1]][1]))] # (rownames/colnames were in the form "participantID_visitID")
colnames(D_V1V3) = id_to_dummy[as.character(sapply(colnames(D_V1V3), function(s) strsplit(s,"_")[[1]][1]))]
D_V1V3 = diag(D_V1V3[id_to_dummy,id_to_dummy])

D_V1V2 = D[grepl("V1",rownames(D)),grepl("V2",colnames(D))]
rownames(D_V1V2) = id_to_dummy[as.character(sapply(rownames(D_V1V2), function(s) strsplit(s,"_")[[1]][1]))]
colnames(D_V1V2) = id_to_dummy[as.character(sapply(colnames(D_V1V2), function(s) strsplit(s,"_")[[1]][1]))]
D_V1V2 = diag(D_V1V2[id_to_dummy,id_to_dummy])


### Make Fig. 1D
if (save_eps) {
  pdf(file=paste(fig_path,"Figure1D.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
D_V1V2_YN = D_V1V2[grepl("YN",names(D_V1V2))]
D_V1V2_YP = D_V1V2[grepl("YP",names(D_V1V2))]
D_V1V2_EN = D_V1V2[grepl("EN",names(D_V1V2))]
D_V1V2_EP = D_V1V2[grepl("EP",names(D_V1V2))]

D_V1V3_YN = D_V1V3[grepl("YN",names(D_V1V3))]
D_V1V3_YP = D_V1V3[grepl("YP",names(D_V1V3))]
D_V1V3_EN = D_V1V3[grepl("EN",names(D_V1V3))]
D_V1V3_EP = D_V1V3[grepl("EP",names(D_V1V3))]

data_list = list(D_V1V2_YN,D_V1V3_YN,
                 D_V1V2_YP,D_V1V3_YP,
                 D_V1V2_EN,D_V1V3_EN,
                 D_V1V2_EP,D_V1V3_EP)
par(mar=c(8,5,2,2))
box_pos = c(1,2,4,5,7,8,10,11)/2.5
name_pos = 0.2+c(mean(box_pos[1],box_pos[2]),
             mean(box_pos[3],box_pos[4]),
             mean(box_pos[5],box_pos[6]),
             mean(box_pos[7],box_pos[8]))
colors = brewer.pal(n=6,name="Set1")[3:4]
boxplot(list(c(),c(),c(),c()),
        at=name_pos,
        las=1,
        names = rep("",4),
        ylab = "",
        ylim=range(data_list))
boxplot(data_list,add=T,
        border = rep(colors,4),
        at=box_pos,
        range = Inf, outline=F,
        names = rep("",8),ann=F,axes=F,
        boxwex=0.3)
beeswarm(data_list,"swarm",add=T,pch=20,cex=0.7,
         at=box_pos)
mtext(side=2,line=4, 'UniFrac(day 0, day 7)',col=colors[1])
mtext(side=2,line=3, 'UniFrac(day 0, day 28)',col=colors[2])
mtext(side=1,line=2,at=name_pos,c("young,\nCMV-","young,\nCMV+","elderly,\nCMV-","elderly,\nCMV+"))
if (save_eps) dev.off()

wilcox.test(c(D_V1V2_YN,D_V1V2_YP),
            c(D_V1V2_EN,D_V1V2_EP))
wilcox.test(c(D_V1V3_YN,D_V1V3_YP),
            c(D_V1V3_EN,D_V1V3_EP))


