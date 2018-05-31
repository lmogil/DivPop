####by Heather E. Wheeler 20170710####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('AFA','CAU') #uncomment for testing
"%&%" = function(a,b) paste(a,b,sep="")
library(dtplyr)
library(data.table)
library(dplyr)

###############################################
### Directories & Variables
pre <- "/home" #~/mount or /home
rna.dir <- pre %&% "/lauren/mesa_expression_files_ens/"
gt.dir <- pre %&% "/lauren/mesa_predixcan_dosages/"
out.dir <- pre %&% "/wheelerlab3/mesa_analyses/GCTA_exp/"
grm.dir <- out.dir %&% "GRMs/"

pop1 <- args[1]
pop2 <- args[2]
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
                        
exp.w.geno <- allexpmat[as.character(allsamplelist$samplelist),] ###get expression of samples with genotypes###
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
  cat(i,"/",length(locallist),"\n")
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
write.table(output,file=tis %&% ".local-h2_gen-corr." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
