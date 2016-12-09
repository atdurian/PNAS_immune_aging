args = commandArgs(trailingOnly=TRUE)

Nseq = as.numeric(args[1])
germline_field = args[2]
files = args[3:(2+(length(args)-2)/2)]

### Packages:
package_path = '/local10G/debourcy/tools/Rpackages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
library('data.table',lib.loc=package_path)
library('stringdist', lib.loc=package_path)
library('snow',lib.loc=package_path)
library('iterators',lib.loc=package_path)
library('foreach',lib.loc=package_path)
library('doSNOW',lib.loc=package_path)
library('SDMTools',lib.loc=package_path)
library('ade4', lib.loc=package_path)
library('seqinr',lib.loc=package_path)
library('shm',lib.loc=package_path) # depreacted, need to install shazam instead
library('Rcpp',lib.loc=package_path)
library('dplyr',lib.loc=package_path)
library('ape', lib.loc=package_path)
library('phangorn', lib.loc=package_path)
library('permute', lib.loc=package_path)
library('vegan', lib.loc=package_path)
library('GUniFrac', lib.loc=package_path)
library('maps', lib.loc=package_path)
library('phytools', lib.loc=package_path)
library('plotrix', lib.loc=package_path) 
muscle_path = '/local10G/debourcy/tools/muscle'

 
### Functions:
subsampleChangeoDbSeqNumberProbDUPCOUNT = function(df,N_final) {
  if (N_final < nrow(df)) {
    indices = sample(1:nrow(df),size=N_final,prob=df$DUPCOUNT)
    return(df[indices,])  
  } else {
    return(df)  
  }
}

extract_Vgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  s=sapply(w,function(x) {paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
  return(as.vector(s))
}
extract_Jgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  return(as.vector(w))
}

rbind.all.columns = function(x,y) {
  if (length(x)==0) {return(y)} else {
  x.diff = setdiff(colnames(x),colnames(y))
  y.diff = setdiff(colnames(y),colnames(x))
  x[,c(as.character(y.diff))] = 0
  y[,c(as.character(x.diff))] = 0
  return(rbind(x,y))
  }
}

seqDist = function(seq1, seq2, dist_model) {
   # [adapted from a function in shazam]
   seq1 = unlist(strsplit(seq1, ""))
   seq2 = unlist(strsplit(seq2, ""))
   valid.idx = !(seq1 %in% c("-", ".")) | !(seq2 %in% c("-","."))
   seq1 = seq1[valid.idx]
   seq2 = seq2[valid.idx]
   d = sapply(1:length(seq1), function(x) {
     dist_model[seq1[x], seq2[x]]
   })
   return(sum(d[d >= 0]))
}

UniFrac_unweighted = function(otu.tab, tree) {
    # [this function is based on the GUniFrac function from the GUniFrac R-package, modified for our purpose]
    if (!is.rooted(tree)) 
        stop("Rooted phylogenetic tree required!")
    otu.tab <- as.matrix(otu.tab)
    row.sum <- rowSums(otu.tab)
    otu.tab <- otu.tab/row.sum
    n <- nrow(otu.tab)
    if (is.null(rownames(otu.tab))) {
        rownames(otu.tab) <- paste("comm", 1:n, sep = "_")
    }
    unshared_branchlength <- matrix(NA, n, n, dimnames = list(rownames(otu.tab),rownames(otu.tab)))
    total_branchlength <- matrix(NA, n, n, dimnames = list(rownames(otu.tab),rownames(otu.tab)))
    if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
        stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\t\tin the OTU table and the tree should match!")
    }
    absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
    if (length(absent) != 0) {
        tree <- drop.tip(tree, absent)
        warning("The tree has more OTU than the OTU table!")
    }
    tip.label <- tree$tip.label
    otu.tab <- otu.tab[, tip.label]
    ntip <- length(tip.label)
    nbr <- nrow(tree$edge)
    edge <- tree$edge
    edge2 <- edge[, 2]
    br.len <- tree$edge.length
    cum <- matrix(0, nbr, n)
    for (i in 1:ntip) {
        tip.loc <- which(edge2 == i)
        cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
        node <- edge[tip.loc, 1]
        node.loc <- which(edge2 == node)
        while (length(node.loc)) {
            cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
            node <- edge[node.loc, 1]
            node.loc <- which(edge2 == node)
        }
    }
    for (i in 2:n) {
        for (j in 1:(i - 1)) {
            cum1 <- cum[, i]
            cum2 <- cum[, j]
            ind <- (cum1 + cum2) != 0
            cum1 <- cum1[ind]
            cum2 <- cum2[ind]
            br.len2 <- br.len[ind]
            I_cum1 <- (cum1 != 0)
            I_cum2 <- (cum2 != 0)
            
            unshared_branchlength[i,j] = unshared_branchlength[j,i] = sum(abs(I_cum1-I_cum2)*abs(cum1-cum2)*br.len2)
            total_branchlength[i,j] = total_branchlength[j,i] = sum((cum1 + cum2)*br.len2)

        }
    }
    return(list(unshared_branchlength = unshared_branchlength, total_branchlength = total_branchlength))
}


### Read data in:
data = data.frame()
for (file in files) {
  data0 = readChangeoDb(file)
  pat = strsplit(file,"_")[[1]][1]
  for (vis in unique(data0$VISIT)) {
    data1 = data0[data0$VISIT==vis,]
    data1$SAMPLE = paste(pat,vis,sep="_")
	consolidation_formula = as.formula(paste("DUPCOUNT ~ SAMPLE + V_CALL_GENOTYPED + J_CALL + JUNCTION_LENGTH + SEQUENCE_IMGT + ",germline_field,sep=""))
    data1 = aggregate(consolidation_formula, data=data1, FUN=sum) 
    data1 = subsampleChangeoDbSeqNumberProbDUPCOUNT(data1,Nseq)
    data = rbind(data,data1)
  }
}


### Make sure gene/allele labeling is consistent between patients if sequence is the same:
reference_seqs = character()
for (file in files) {
pat = strsplit(file,"_")[[1]][1]
reference_seqs0 = read.FASTA(paste(pat,"personal_hIGHV.fasta",sep="_"))
reference_seqs0 = as.character(reference_seqs0)
reference_seqs0 = sapply(reference_seqs0,function(s) paste(s,collapse=""))
reference_seqs = c(reference_seqs,reference_seqs0)
}
reference_seqs = gsub("-",".",reference_seqs)
truncated_reference_seqs = sapply(1:length(reference_seqs),function(s) {
  allele = names(reference_seqs)[s]
  sequence = reference_seqs[s]
  if (grepl("IGHV1",allele)) {start=276}       # start position of sequence given the primers that were used
  if (grepl("IGHV2",allele)) {start=280}
  if (grepl("IGHV3",allele)) {start=284}
  if (grepl("IGHV4",allele)) {start=281}
  if (grepl("IGHV5",allele)) {start=279}
  if (grepl("IGHV6",allele)) {start=281}
  if (grepl("IGHV7",allele)) {start=288}
  return(substr(sequence,start,nchar(sequence)))  
})
map = character()
for (i in 1:length(truncated_reference_seqs)) {
  allele = names(truncated_reference_seqs)[i]
  seq = truncated_reference_seqs[allele]
  map[allele] = sort(names(truncated_reference_seqs)[truncated_reference_seqs == seq])[1]
}
data$V_CALL_GENOTYPED = sapply(strsplit(as.character(data$V_CALL_GENOTYPED),","),
                               function(s) paste(unique(map[s]),collapse=","))

### group by V-J-length combination:
data$V_SEGMENT = extract_Vgene(data$V_CALL_GENOTYPED)
data$J_SEGMENT = sapply(data$J_CALL,function(s) strsplit(s," ")[[1]][2])
data$J_SEGMENT = extract_Jgene(data$J_SEGMENT)
data$VJL_GROUP = paste(data$V_SEGMENT,data$J_SEGMENT,data$JUNCTION_LENGTH,sep="_")
        

### Calculate branch lengths, group by group:
samples = unique(data$SAMPLE)
sum_unshared_branchlength = sum_total_branchlength = matrix(0,nrow=length(samples),ncol=length(samples)) 
vjl_groups = unique(data$VJL_GROUP)
for (i in 1:length(vjl_groups)) {
  vjl_group = vjl_groups[i]
  subdata = data[data$VJL_GROUP==vjl_group,]
  
  ## Germlines based on IMGT are aligned from the same start position, but lengths may differ.
  ## Pad sequences to the same length: 
  germline_seqs = subdata[,germline_field]
  germ_seq_list = strsplit(tolower(germline_seqs),"")
  names(germ_seq_list) = 1:length(germline_seqs)
  seqLength = max(sapply(germ_seq_list,length))
  if (!all(lapply(germ_seq_list,length)==seqLength)) {
    germ_seq_list = lapply(germ_seq_list,function(s) c(s,rep("n",seqLength - length(s))))
  }

  ## Mask positions where germlines differ and record resulting 'consensus' germline:
  germ_seq_list = lapply(germ_seq_list,function(s) gsub("\\.","-",s))
  germ_seq_matrix = matrix(unlist(germ_seq_list), nrow = length(germ_seq_list), byrow = TRUE) 
  germ_seq_matrix[,!apply(germ_seq_matrix,2,function(s) all(s==s[1]))]='n'
  consensusGermline = germ_seq_matrix[1,]
  consensusGermline = toupper(paste(consensusGermline,collapse=""))
  data$CONSENSUS_GERMLINE[data$VJL_GROUP==vjl_group] = consensusGermline
  
  ## Apply mask (and length padding) from the consensus germline to the observed sequences as well:
  seqs = c(consensusGermline,subdata$SEQUENCE_IMGT)
  names(seqs) = c("germline",paste("seq",1:length(seqs[-1]),sep=""))
  seq_list = strsplit(tolower(seqs),"")
  names(seq_list) = names(seqs)
  seqLength = max(sapply(seq_list,length))
  if (!all(lapply(seq_list,length)==seqLength)) {
    seq_list = lapply(seq_list,function(s) c(s,rep("n",seqLength - length(s))))
  }
  seq_list = lapply(seq_list,function(s) gsub("\\.","-",s))
  seq_matrix = matrix(unlist(seq_list), nrow = length(seq_list), byrow = TRUE) 
  rownames(seq_matrix) = names(seqs)
  Nindices = which(seq_matrix["germline",]=='n')
  seq_matrix[,Nindices] = 'n'
  
  ## Compute matrix of distances between sequences in the V-J-length group:
  seq_matrix[,apply(seq_matrix,2,function(s) any(s %in% c("n","-")))]='n'
  alnd_seqs = toupper(apply(seq_matrix, 1, paste, collapse = ""))
  names(alnd_seqs) = names(seqs)
  subdata$SEQUENCE_MOD = alnd_seqs[-1]
  alnd_seqs = alnd_seqs[sort(names(alnd_seqs))]
  alnd_seqs = alnd_seqs[!duplicated(alnd_seqs)]
  dm = sapply(1:length(alnd_seqs),function(i) c(rep.int(0,i - 1),sapply(i:length(alnd_seqs),function(j) {
                                               seqDist(alnd_seqs[i],alnd_seqs[j],dist_model=HS1FDistance)
                                             })))
  dm = dm + t(dm)
  rownames(dm) = colnames(dm) = names(alnd_seqs)
  DM = as.dist(dm)
  
  if (length(DM)!=0) {
  
    ## Make phylogenetic tree from distance matrix:
    if (length(DM)==1) seq_tree=upgma(DM) else seq_tree=NJ(DM)       # If there is only 1 sequence besides the germline, the tree is trivial and it does not matter which tree-building method we use, so we use the upgma function because NJ would give an error
    nn = which(seq_tree$tip.label=="germline")
    seq_tree$edge.length[seq_tree$edge.length < 0] = 0
    seq_tree = reroot(seq_tree,nn,seq_tree$edge.length[which(seq_tree$edge[,2]==nn)])

    ## Make incidence table of sequences in the different samples:
    subdata$SEQ_NAME = sapply(subdata$SEQUENCE_MOD,function(s) names(alnd_seqs)[alnd_seqs==s])                              
    incidence_table = data.frame()
    for (sample in samples) {
      incidence_line = as.data.frame(t(c(table(subdata$SEQ_NAME[subdata$SAMPLE==sample]))))
      rownames(incidence_line) = sample
      if (!("germline" %in% names(incidence_line))) incidence_line$germline = 1
      incidence_table = rbind.all.columns(incidence_table,incidence_line)
    }
    
    ## Calculate unshared and total branch lengths for the V-J-length group:
    uf = UniFrac_unweighted(incidence_table, seq_tree)       
    sum_unshared_branchlength = sum_unshared_branchlength + uf$unshared_branchlength
    sum_total_branchlength = sum_total_branchlength + uf$total_branchlength
  }
}

### Having calculated contributions from all V-J-length groups, calculate overall RUBL:
RUBL = sum_unshared_branchlength/sum_total_branchlength
diag(RUBL) = 0


### Write results:
outfile = paste("RUBL.Nseq-",Nseq,".",germline_field,".tab",sep="")
write.table(as.data.frame(RUBL),file=outfile)




