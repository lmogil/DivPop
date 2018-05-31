####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','AFA')
"%&%" = function(a,b) paste(a,b,sep="")
library(dtplyr)
library(data.table)
library(dplyr)

###############################################
### Directories & Variables
pre <- "/home" #~/mount or /home
rna.dir <- pre %&% "/lauren/mesa_expression_files_ens/"
gt.dir <- pre %&% "/lauren/mesa_predixcan_dosages/"
out.dir <- pre %&% "/wheelerlab3/mesa_analyses/BSLMM_exp/"

pop <- args[2]
tis <- pop %&% "_MESA_Nk-20"  
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom

getquant <- function(x) quantile(x,c(0.5,0.025,0.975)) ##pulls the median and the 95% credible sets

################################################
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
gencode <- gencode[gencode[,2]==chrname,] ##pull genes on chr of interest
expmat <- expmat[,intersect(colnames(expmat),rownames(gencode))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(expmat) ###samples with exp data###

famfile <- gt.dir %&% "samples_gen3_" %&% pop %&% ".txt" ###samples with gt data###
fam <- read.table(famfile)
samplelist <- intersect(fam$V1,expsamplelist)
                        
exp.w.geno <- expmat[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

gtfile <- gt.dir %&% 'imp2_' %&% pop %&% '_chr' %&% chrom %&% '_ref1kg_dosage.txt.gz'
gtX <- fread('zcat ' %&% gtfile)

snpinfo <- data.frame(gtX[,1:5])
colnames(snpinfo) <- c("chr","snp","bp","ref","alt")
rownames(snpinfo) <- snpinfo$snp

gtX <- as.matrix(gtX[,-1:-6])
rownames(gtX) <- snpinfo$snp
colnames(gtX) <- fam$V1
X <- gtX[,samplelist]

resultsarray <- array(0,c(length(explist),19))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
dimnames(resultsarray)[[2]] <- resultscol

working100K <- out.dir %&% "working_" %&% tis %&% "_exp_BSLMM-s100K_iterations_chr" %&% chrom %&% "b" %&% "_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=19,sep="\t")

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  geneinfo <- gencode[gene,]
  chr <- geneinfo$chr
  c <- substr(chr,4,5)
  start <- geneinfo$start - 1e6 ### 1Mb lower bound for cis-eQTLS
  end <- geneinfo$end + 1e6 ### 1Mb upper bound for cis-eQTLs
  chrsnps <- subset(snpinfo,snpinfo[,1]==chrname) ### pull snps on same chr
  cissnps <- subset(chrsnps,chrsnps[,3]>=start & chrsnps[,3]<=end) ### pull cis-SNP info
  cisgenos <- X[intersect(rownames(X),cissnps[,2]),,drop=FALSE] ### pull cis-SNP genotypes
  if(dim(cisgenos)[1] > 0){
    cisgenos <- mutate(data.frame(cisgenos),snp=rownames(cisgenos))
    cisbim <- inner_join(snpinfo,cisgenos,by='snp')
    annotfile <- cbind(cisbim[,2],cisbim[,3],c)
    genofile <- cbind(cisbim[,2],cisbim[,4:dim(cisbim)[2]])
    phenofile <- data.frame(exp.w.geno[,gene])

    write.table(annotfile, file=out.dir %&% "tmp2.annot." %&% pop %&% chrom %&% "b" , quote=F, row.names=F, col.names=F, sep=",")
    write.table(genofile, file=out.dir %&% "tmp2.geno." %&% pop %&% chrom %&% "b" , quote=F, row.names=F, col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmp2.pheno." %&% pop %&% chrom %&% "b" , quote=F, row.names=F, col.names=F, sep=",")

    runBSLMM <- "/usr/local/bin/gemma -g " %&% out.dir %&% "tmp2.geno." %&% pop %&% chrom %&% "b" %&%  " -p " %&% out.dir %&% "tmp2.pheno." %&% 
      pop %&% chrom %&% "b" %&% " -a " %&% out.dir %&% "tmp2.annot." %&% pop %&% chrom %&% "b" %&% " -bslmm 1 -seed 12345 -s 100000 -o tmp2." %&% pop %&% chrom %&% "b"  
    system(runBSLMM)

    hyp <- read.table(out.dir %&% "output/tmp2." %&% pop %&% chrom %&% "b" %&% ".hyp.txt",header=T)
    hyp50 <- hyp[(dim(hyp)[1]/2+1):dim(hyp)[1],] #take second half of sampling iterations
    quantres <- apply(hyp50,2,getquant)
    res <- c(gene,quantres[1,],quantres[2,],quantres[3,])

  }else{
    res <- c(gene,rep(NA,18))
  }
  names(res) <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975") 
  resultsarray[gene,] <- res
  write(res,file=working100K,ncolumns=19,append=T,sep="\t")
}

write.table(resultsarray,file=out.dir %&% tis %&% "_exp_BSLMM-s100K_iterations_chr" %&% chrom %&% "b" %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
