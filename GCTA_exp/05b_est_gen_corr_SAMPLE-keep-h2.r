####by Heather E. Wheeler 20170710####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('AFA','CAU') #uncomment for testing
"%&%" = function(a,b) paste(a,b,sep="")
library(dtplyr)
library(data.table)
library(dplyr)

### shuffle the gene exp values without replacement before estimating genetic correlation
### include a seed variable to run different shufflings

###############################################
### Directories & Variables
pre <- "/home" #~/mount or /home
rna.dir <- pre %&% "/lauren/mesa_expression_files_ens/"
gt.dir <- pre %&% "/lauren/mesa_predixcan_dosages/"
out.dir <- pre %&% "/wheelerlab3/mesa_analyses/GCTA_exp/"
grm.dir <- out.dir %&% "GRMs/"

pop1 <- args[1]
pop2 <- args[2]
seednum <- as.double(args[3])

### For testing, comment out when running
#pop1 <- 'AFA'
#pop2 <- 'CAU'
#seednum <- '42'
###############

set.seed(seednum)
tis <- pop1 %&% "-" %&% pop2 %&% "_MESA_Nk-20" 

################################################
#read in expression data to get appropriate sample list, combine pops
for(pop in c(pop1, pop2)){
  expfile <- rna.dir %&% pop %&% "_MESA_Epi_GEX_data_sidno_Nk-20.txt"
  expdata <- fread(expfile,header=TRUE)
  expmat <- t(as.matrix(expdata[,-1])) #want people in rows, genes in cols
  rownames(expmat) <- colnames(expdata)[-1]
  colnames(expmat) <- expdata$PROBE_ID
  if(exists('allexpmat')){
    allexpmat <- rbind(allexpmat, expmat)
  }else{
    allexpmat <- expmat
  }
}

gencodefile <- rna.dir %&% 'MESA_GEX_ilmn_gencode.v18.annotation.txt'
gencode <- data.frame(fread(gencodefile))
#remove ensid duplicates, keep first entry
gencode <- gencode[!duplicated(gencode$gene_id),]
rownames(gencode) <- gencode[,5]
allexpmat <- allexpmat[,intersect(colnames(allexpmat),rownames(gencode))] ###pull gene expression data w/gene info

expsamplelist <- rownames(allexpmat) ###samples with exp data###

for(pop in c(pop1, pop2)){
  famfile <- gt.dir %&% "samples_gen3_" %&% pop %&% ".txt" ###samples with gt data###
  fam <- read.table(famfile)
  samplelist <- intersect(fam$V1,expsamplelist)
  samplelist <- cbind(samplelist,pop)
  if(exists('allsamplelist')){
    allsamplelist <- rbind(allsamplelist, data.frame(samplelist))
  }else{
    allsamplelist <- data.frame(samplelist)
  }
}

pop1list <- as.character(dplyr::filter(allsamplelist,pop==pop1)$samplelist)
pop2list <- as.character(dplyr::filter(allsamplelist,pop==pop2)$samplelist)
                        
####permutation scheme to keep same h2 distribution as in Brown et al. (popcorn paper)####

#get pop1 expression of samples with genotypes
pop1.exp.w.geno <- allexpmat[pop1list,]
#split into 2 sets, each with 1/2 the genes, shuffle ID labels on first set
nga <- floor(dim(pop1.exp.w.geno)[2]/2) #num genes in a
p1a <- pop1.exp.w.geno[,1:nga] #first set
p1b <- pop1.exp.w.geno[,(nga+1):dim(pop1.exp.w.geno)[2]] #second set
#shuffle p1a gene expression (shuffle sample IDs, keep same distribution of exp per gene, 
#keeps exp correlation together, will separate from genotypes)
shuffle <- sample(rownames(p1a)) #shuffle rownames (IDs)
rownames(p1a) <- shuffle #assign shuffled IDs to exp.w.geno
p1a <- p1a[ order(as.numeric(row.names(p1a))),] #sort by rownames
p1b <- p1b[ order(as.numeric(row.names(p1b))),] #sort by rownames
pop1.exp.w.geno <- cbind(p1a, p1b)

#get pop2 expression of samples with genotypes
pop2.exp.w.geno <- allexpmat[pop2list,]
#split into 2 sets, each with 1/2 the genes, shuffle ID labels on second set
nga <- floor(dim(pop2.exp.w.geno)[2]/2) #num genes in a
p2a <- pop2.exp.w.geno[,1:nga] #first set
p2b <- pop2.exp.w.geno[,(nga+1):dim(pop2.exp.w.geno)[2]] #second set
#shuffle p2b gene expression (shuffle sample IDs, keep same distribution of exp per gene, 
#keeps exp correlation together, will separate from genotypes)
shuffle <- sample(rownames(p2b)) #shuffle rownames (IDs)
rownames(p2b) <- shuffle #assign shuffled IDs to exp.w.geno
p2a <- p2a[ order(as.numeric(row.names(p2a))),] #sort by rownames
p2b <- p2b[ order(as.numeric(row.names(p2b))),] #sort by rownames
pop2.exp.w.geno <- cbind(p2a, p2b)

exp.w.geno <- rbind(pop1.exp.w.geno,pop2.exp.w.geno) ###combine pop1 and pop2###
explist <- colnames(exp.w.geno)
nsubj <- dim(exp.w.geno)[1]


### Get gene set to analyze, i.e. those with GRMs
system("ls " %&% grm.dir %&% pop1 %&% "-" %&% pop2 %&% "-local*id > tmp." %&% pop1 %&% pop2)
system("awk -F \"" %&% pop1 %&% "-" %&% pop2 %&% "-local-\" \'{print $2}\' < tmp." %&% pop1 %&% pop2 %&%
         " | awk -F \".grm\" \'{print $1}\' >localGRM.list." %&% pop1 %&% "-" %&% pop2)
localfile <- "localGRM.list."	%&% pop1 %&% "-" %&% pop2
locallist <- as.data.frame(scan(localfile,"character")) ##list of local GRMs
colnames(locallist)<-'gene'
stopifnot(locallist$gene %in% explist)

### Output matrix
loc.mat <- matrix(0,nrow=length(locallist$gene),ncol=13)
colnames(loc.mat) <- c("pop-data","N","ensid","gene","chr","start","end","pop1.h2","pop1.se","pop2.h2","pop2.se","rG","rG.se")

for(i in 1:length(locallist$gene)){
  cat(i,"/",dim(locallist)[1],"\n")
  gene <- as.character(locallist$gene[i])
  geneinfo <- gencode[gene,]

  #output expression pheno for gcta
  geneexp <- data.frame(cbind(IID=rownames(exp.w.geno),allexp=exp.w.geno[,gene]), stringsAsFactors = FALSE)
  geneexp <- mutate(geneexp, pop1list=ifelse(IID %in% pop1list, allexp, NA))
  geneexp <- mutate(geneexp, pop2list=ifelse(IID %in% pop2list, allexp, NA)) 
  geneexp <- dplyr::select(geneexp, -allexp)
  rownames(geneexp) <- geneexp$IID
  write.table(geneexp, file="tmp." %&% pop1 %&% pop2 %&% ".pheno", col.names=F, quote=F) #output pheno for gcta
  
  ## Y ~ localGRM
  runLOC <- "/usr/local/bin/gcta64 --grm " %&% grm.dir %&% pop1 %&% "-" %&% pop2 %&% "-local-" %&% gene %&% 
    " --reml-bivar --pheno tmp." %&% pop1 %&% pop2 %&% ".pheno --out tmp." %&% pop1 %&% pop2 %&% " --thread-num 10"
  system(runLOC)
  if(file.exists("tmp." %&% pop1 %&% pop2 %&% ".hsq")==TRUE){
    hsq <- scan("tmp." %&% pop1 %&% pop2 %&% ".hsq","character")
    res <- c(tis, nsubj, gene, geneinfo$gene_name, geneinfo$chr, geneinfo$start, geneinfo$end, hsq[26], hsq[27], hsq[29], hsq[30], hsq[32], hsq[33])
    system("rm tmp." %&% pop1 %&% pop2 %&% ".hsq")
  }else{
    res <- c(tis, nsubj, gene, geneinfo$gene_name, geneinfo$chr, geneinfo$start, geneinfo$end, NA, NA, NA, NA, NA, NA) ##gcta did not converge
  }
  loc.mat[i,] <- res 
}

output <- data.frame(loc.mat) %>% arrange(desc(pop1.h2))
write.table(output,file=tis %&% ".local-h2_gen-corr.SAMPLE-keep-h2." %&% seednum %&% "." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
