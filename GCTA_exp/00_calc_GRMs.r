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
out.dir <- pre %&% "/wheelerlab3/mesa_analyses/GCTA_exp/"
grm.dir <- out.dir %&% "GRMs/"

pop <- args[2]
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom

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
gencode <- gencode[gencode[,2]==chrname,] ##pull genes on chr of interest
expmat <- expmat[,intersect(colnames(expmat),rownames(gencode))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(expmat) ###samples with exp data###

famfile <- gt.dir %&% "samples_gen3_" %&% pop %&% ".txt" ###samples with gt data###
fam <- read.table(famfile)
samplelist <- intersect(fam$V1,expsamplelist)
                        
exp.w.geno <- expmat[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

gtfile <- gt.dir %&% 'imp2_' %&% pop %&% '_chr' %&% chrom %&% '_ref1kg_dosage.txt.gz'
#read gtfile depending if running remotely or locally
ifelse(pre == "/home",gtX <- fread('zcat ' %&% gtfile),gtX <- fread('gzcat ' %&% gtfile))

snpinfo <- data.frame(gtX[,1:6])
colnames(snpinfo) <- c("chr","snp","bp","ref","alt","altaf")
rownames(snpinfo) <- snpinfo$snp

gtX <- as.matrix(gtX[,-1:-6])
rownames(gtX) <- snpinfo$snp
colnames(gtX) <- fam$V1
X <- gtX[,samplelist]

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
    cisgenos <- t(cisgenos)
    mldosecols <- dplyr::mutate(data.frame(cisgenos),id=rownames(cisgenos) %&% "->" %&% rownames(cisgenos),ml="MLDOSE") %>%
      dplyr::select(id,ml)
    mldose <- cbind(mldosecols, cisgenos)
    
    mlinfosnps <- data.frame(snp=colnames(cisgenos))
    mlinfocols <- left_join(mlinfosnps,snpinfo,by='snp') %>% mutate(Freq1=(1-altaf),MAF=ifelse(altaf<0.5,altaf,1-altaf))
    mlinfo <- dplyr::select(mlinfocols,snp,ref,alt,Freq1,MAF) %>% mutate(quality=1,Rsq=1) #dummy variables
    colnames(mlinfo) <- c("SNP","Al1","Al2","Freq1","MAF","Quality","Rsq")
    
    write.table(mlinfo, file=out.dir %&% "tmp2." %&% pop %&% chrom %&% ".mlinfo", quote=F, row.names=F, sep="\t")
    write.table(mldose, file=out.dir %&% "tmp2." %&% pop %&% chrom %&% ".mldose" , quote=F, row.names=F, col.names=F, sep=" ")
    
    runGCTA <- "/usr/local/bin/gcta64 --dosage-mach " %&% out.dir %&% "tmp2." %&% pop %&% chrom %&% ".mldose " %&%
      out.dir %&% "tmp2." %&% pop %&% chrom %&% ".mlinfo --make-grm-bin --thread-num 10 --out " %&% 
      grm.dir %&% pop %&% "-local-" %&% gene 
    system(runGCTA)

  }else{
    cat("no SNPs for", pop, gene)
  }
}
