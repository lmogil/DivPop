####by Heather E. Wheeler 20180330####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('9','CAU')
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)
library(data.table)

###############################################
### Directories & Variables
pre <- "/home" #~/mount or /home
rna.dir <- pre %&% "/lauren/mesa_expression_files_ens/"
gt.dir <- pre %&% "/lauren/mesa_predixcan_dosages/"
out.dir <- pre %&% "/wheelerlab3/mesa_analyses/BVSR_exp/"

pop <- args[2]
tis <- pop %&% "_MESA_Nk-20"  
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom

getquant <- function(x) quantile(x,c(0.5,0.025,0.975),na.rm = TRUE) ##pulls the median and the 95% credible sets

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

resultsarray <- array(0,c(length(explist),28))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","h50","hh50","p50","snp50","ar50","war50","sigma50","bf50","like50","h025","hh025","p025","snp025","ar025","war025","sigma025","bf025","like025","h975","hh975","p975","snp975","ar975","war975","sigma975","bf975","like975")
dimnames(resultsarray)[[2]] <- resultscol

working100K <- out.dir %&% "working_" %&% tis %&% "_exp_BVSR-s100K_iterations_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=28,sep="\t")

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
    genofile <- cbind(cisbim[,2],cisbim[,4:dim(cisbim)[2]])
    phenofile <- data.frame(exp.w.geno[,gene])

    write.table(genofile, file=out.dir %&% "tmp2.geno." %&% pop %&% chrom , quote=F, row.names=F, col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmp2.pheno." %&% pop %&% chrom , quote=F, row.names=F, col.names=F, sep=",")

    runBVSR <- "/usr/local/bin/pimass-lin -g " %&% out.dir %&% "tmp2.geno." %&% pop %&% chrom %&%  " -p " %&% out.dir %&% "tmp2.pheno." %&% 
      pop %&% chrom %&% " -r 12345 -w 10000 -s 100000 -num 10 -o tmp2." %&% pop %&% chrom  
    system(runBVSR)

    hyp <- read.table(out.dir %&% "output/tmp2." %&% pop %&% chrom %&% ".path.txt",header=T)
    #hyp50 <- hyp[(dim(hyp)[1]/2+1):dim(hyp)[1],] #take second half of sampling iterations, like BSLMM
    quantres <- apply(hyp,2,getquant)
    res <- c(gene,quantres[1,],quantres[2,],quantres[3,])

  }else{
    res <- c(gene,rep(NA,27))
  }
  names(res) <- c("gene","h50","hh50","p50","snp50","ar50","war50","sigma50","bf50","like50","h025","hh025","p025","snp025","ar025","war025","sigma025","bf025","like025","h975","hh975","p975","snp975","ar975","war975","sigma975","bf975","like975") 
  resultsarray[gene,] <- res
  write(res,file=working100K,ncolumns=28,append=T,sep="\t")
}

write.table(resultsarray,file=out.dir %&% tis %&% "_exp_BVSR-s100K_iterations_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
