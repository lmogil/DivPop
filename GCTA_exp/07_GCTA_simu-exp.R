####Goal
# simulate gene expression phenotypes with the same h2 distribution as the real data
# pull h2 estimates from origial rG estimation and use all SNPs in the gene region
# to simulate exp phenotypes with the same h2
# save results as an exp matrix to use with 03_est_gen_corr.r, to match: want people in cols, genes in rows

####by Heather E. Wheeler 20170710####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('AFA','CAU') #uncomment for testing
"%&%" = function(a,b) paste(a,b,sep="")
library(dtplyr)
library(data.table)
library(dplyr)
library(tibble)

###############################################
### Directories & Variables
pre <- "/home" #~/mount or /home
out.dir <- pre %&% "/wheelerlab3/mesa_analyses/GCTA_exp/"
rna.dir <- pre %&% "/lauren/mesa_expression_files_ens/"
gt.dir <- out.dir %&% "MESA_plink_no_dupvar/"


pop1 <- args[1]
pop2 <- args[2]
simnum <- args[3]
#pop1 <- 'AFA'
#pop2 <- 'CAU'
#simnum <- 1
resfile <- out.dir %&% pop1 %&% "-" %&% pop2 %&% "_MESA_Nk-20.local-h2_gen-corr.2017-07-21.txt"

res <- fread(resfile) %>% mutate(chrnum=as.integer(substr(chr,4,5)),pop1.snpcount=0,pop2.snpcount=0)

### Output simulated expression matrices
pop1.fam <- fread(gt.dir %&% "MESA_" %&% pop1 %&% "_w_expression.fam")
pop1.mat <- matrix(0,nrow=dim(pop1.fam)[1],ncol=dim(res)[1])
rownames(pop1.mat) <- pop1.fam$V2
colnames(pop1.mat) <- res$ensid

pop2.fam <- fread(gt.dir %&% "MESA_" %&% pop2 %&% "_w_expression.fam")
pop2.mat <- matrix(0,nrow=dim(pop2.fam)[1],ncol=dim(res)[1])
rownames(pop2.mat) <- pop2.fam$V2
colnames(pop2.mat) <- res$ensid

for(i in 1:dim(res)[1]){
    info <- res[i,]
  ###run pop1 simulation
  runPLINK <- "plink --bfile " %&% gt.dir %&% "MESA_" %&% pop1 %&% "_w_expression_maf0.01 --chr " %&% info$chrnum %&% 
    " --from-bp " %&% (max(1,(info$start - 1e6))) %&% " --to-bp " %&% (info$end + 1e6) %&% " --make-bed --out tmp" %&% simnum 
  system(runPLINK)
  #get SNP count and make tmp.causal.snplist of all snps in gene region
  if(file.exists("tmp" %&% simnum %&% ".bim")==TRUE){
    bim <- fread("tmp" %&% simnum %&% ".bim")
    res$pop1.snpcount[i] <- dim(bim)[1]
    write.table(bim$V2,"tmp" %&% simnum %&% ".causal.snplist",quote=F,row.names = F,col.names = F)
    #simulate exp phenotype with same h2 as real data
    runGCTA <- "gcta64 --bfile tmp" %&% simnum %&% " --simu-qt --simu-causal-loci tmp" %&% simnum %&% ".causal.snplist --simu-hsq " %&% info$pop1.h2 %&% 
      " --simu-rep 1 --out tmp" %&% simnum %&% "." %&% pop1
    system(runGCTA)
    system("rm tmp" %&% simnum %&% ".bim")
    #add to output matrix
    if(file.exists("tmp" %&% simnum %&% "." %&% pop1 %&% ".phen")){
      phen <- fread("tmp" %&% simnum %&% "." %&% pop1 %&% ".phen")
      pop1.mat[,i] <- phen$V3
      system("rm tmp" %&% simnum %&% "." %&% pop1 %&% ".phen")
    }
  }
      
  ###run pop2 simulation
  runPLINK <- "plink --bfile " %&% gt.dir %&% "MESA_" %&% pop2 %&% "_w_expression_maf0.01 --chr " %&% info$chrnum %&% 
    " --from-bp " %&% (max(1,(info$start - 1e6))) %&% " --to-bp " %&% (info$end + 1e6) %&% " --make-bed --out tmp" %&% simnum
  system(runPLINK)
  #get SNP count and make tmp.causal.snplist of all snps in gene region
  if(file.exists("tmp" %&% simnum %&% ".bim")==TRUE){
    bim <- fread("tmp" %&% simnum %&% ".bim")
    res$pop2.snpcount[i] <- dim(bim)[1]
    write.table(bim$V2,"tmp" %&% simnum %&% ".causal.snplist",quote=F,row.names = F,col.names = F)
    #simulate exp phenotype with same h2 as real data
    runGCTA <- "gcta64 --bfile tmp" %&% simnum %&% " --simu-qt --simu-causal-loci tmp" %&% simnum %&% ".causal.snplist --simu-hsq " %&% info$pop2.h2 %&% 
      " --simu-rep 1 --out tmp" %&% simnum %&% "." %&% pop2
    system(runGCTA)
    system("rm tmp" %&% simnum %&% ".bim")
    if(file.exists("tmp" %&% simnum %&% "." %&% pop2 %&% ".phen")){
      #add to output matrix
      phen <- fread("tmp" %&% simnum %&% "." %&% pop2 %&% ".phen")
      pop2.mat[,i] <- phen$V3
      system("rm tmp" %&% simnum %&% "." %&% pop2 %&% ".phen")
      }
  }
}

tpop1 <- t(pop1.mat)
tpop1df <- as.data.frame(tpop1) %>% rownames_to_column(var="PROBE_ID")
fwrite(tpop1df, file=out.dir %&% "MESA_simulated_exp/" %&% pop1 %&% "-" %&% pop2 %&% 
         "_MESA_Nk-20.local-h2_gen-corr.2017-07-21_" %&% pop1 %&% "-EXP-sim" %&% simnum %&% ".txt", sep="\t")

tpop2 <- t(pop2.mat)
tpop2df <- as.data.frame(tpop2) %>% rownames_to_column(var="PROBE_ID")
fwrite(tpop2df, file=out.dir %&% "MESA_simulated_exp/" %&% pop1 %&% "-" %&% pop2 %&% 
         "_MESA_Nk-20.local-h2_gen-corr.2017-07-21_" %&% pop2 %&% "-EXP-sim" %&% simnum %&% ".txt", sep="\t")

#fwrite(res, file=out.dir %&% pop1 %&% "-" %&% pop2 %&% "_MESA_Nk-20.local-h2_gen-corr.2017-07-21_with_SNP_count.txt",na=NA,sep="\t")
