####by Heather E. Wheeler 20170710####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('AFA') #uncomment for testing
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

pop <- args[1]
tis <- pop %&% "_MESA_Nk-20" 

################################################
#read in expression data to get appropriate sample list
expfile <- rna.dir %&% pop %&% "_MESA_Epi_GEX_data_sidno_Nk-20.txt"
expdata <- fread(expfile,header=TRUE)
expmat <- t(as.matrix(expdata[,-1])) #want people in rows, genes in cols
rownames(expmat) <- colnames(expdata)[-1]
colnames(expmat) <- expdata$PROBE_ID

gencodefile <- rna.dir %&% 'MESA_GEX_ilmn_gencode.v18.annotation.txt'
gencode <- data.frame(fread(gencodefile))
#remove ensid duplicates, keep first entry
gencode <- gencode[!duplicated(gencode$gene_id),]
rownames(gencode) <- gencode[,5]
expmat <- expmat[,intersect(colnames(expmat),rownames(gencode))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(expmat) ###samples with exp data###

famfile <- gt.dir %&% "samples_gen3_" %&% pop %&% ".txt" ###samples with gt data###
fam <- read.table(famfile)
samplelist <- intersect(fam$V1,expsamplelist)
                        
exp.w.geno <- expmat[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)
nsubj <- dim(exp.w.geno)[1]

### Get gene set to analyze, i.e. those with GRMs
system("ls " %&% grm.dir %&% pop %&% "-local*id > tmp." %&% pop)
system("awk -F \"" %&% pop %&% "-local-\" \'{print $2}\' < tmp." %&% pop %&% " | awk -F \".grm\" \'{print $1}\' >localGRM.list." %&% pop)
localfile <- "localGRM.list."	%&% pop
locallist <- as.data.frame(scan(localfile,"character")) ##list of local GRMs
colnames(locallist)<-'gene'
stopifnot(locallist$gene %in% explist)

### Output matrix
loc.mat <- matrix(0,nrow=length(locallist$gene),ncol=10)
colnames(loc.mat) <- c("exp-data","N","ensid","gene","chr","start","end","local.h2","local.se","local.p")

for(i in 1:length(locallist$gene)){
  cat(i,"/",length(locallist),"\n")
  gene <- as.character(locallist$gene[i])
  geneinfo <- gencode[gene,]

  #output expression pheno for gcta
  geneexp <- cbind(rownames(exp.w.geno),exp.w.geno[,gene])
  write.table(geneexp, file="tmp.pheno." %&% pop, col.names=F, quote=F) #output pheno for gcta
  
  ## Y ~ localGRM
  runLOC <- "/usr/local/bin/gcta64 --grm " %&% grm.dir %&% pop %&% "-local-" %&% gene %&% " --reml --pheno tmp.pheno." %&% 
    pop %&% " --out tmp." %&% pop %&% " --thread-num 10"
  system(runLOC)
  if(file.exists("tmp." %&% pop %&% ".hsq")==TRUE){
    hsq <- scan("tmp." %&% pop %&% ".hsq","character")
    res <- c(tis, nsubj, gene, geneinfo$gene_name, geneinfo$chr, geneinfo$start, geneinfo$end, hsq[14], hsq[15], hsq[25])
    system("rm tmp." %&% pop %&% ".hsq")
  }else{
    res <- c(tis, nsubj, gene, geneinfo$gene_name, geneinfo$chr, geneinfo$start, geneinfo$end, NA, NA, NA) ##gcta did not converge
  }
  loc.mat[i,] <- res 
}

output <- data.frame(loc.mat) %>% arrange(desc(local.h2))
write.table(output,file=tis %&% ".local-h2." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
