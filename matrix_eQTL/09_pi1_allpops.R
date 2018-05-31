library(qvalue)
library(data.table)
library(dplyr)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")

#test qvalue package
#data(hedenfalk)
#pvalues <- hedenfalk$p
#qobj <- qvalue(p = pvalues)
#hist(qobj)
#print(qobj$pi0)


poplist2 <- c('AFA','HIS','CAU','MEX','YRI','GEU','FHS')
poplist <- c('AFA','HIS','CAU','AFHI','ALL')
#poplist <- c('AFA','YRI')
Nk<-c(0,10,20,30)
Nk2<-c(0,10,20,30)

npops <- length(poplist)
npops2 <- length(poplist2)
nks <- length(Nk)
nks2 <- length(Nk2)
pi1_matrix <- matrix(NA,nrow=npops*nks,ncol=npops2*nks2)

rownames(pi1_matrix) <- c('AFA0','AFA10','AFA20','AFA30','HIS0','HIS10','HIS20','HIS30','CAU0','CAU10','CAU20','CAU30','AFHI0','AFHI10','AFHI20','AFHI30','ALL0','ALL10','ALL20','ALL30')

colnames(pi1_matrix) <- c('AFA0','AFA10','AFA20','AFA30','HIS0','HIS10','HIS20','HIS30','CAU0','CAU10','CAU20','CAU30','MEX0','MEX10','MEX20','MEX30','YRI0','YRI10','YRI20','YRI30','GEU0','GEU10','GEU20','GEU30','FHS0','FHS10','FHS20','FHS30')


for(i in 1:length(poplist)){
        for(n in 1:length(Nk)){
                refinfile<-fread("/home/lauren/files_for_revisions_plosgen/meqtl_results/GEU/"%&% poplist[i]%&%"_meQTL_results_"%&% Nk[n] %&%"_Nk_3_PCs.txt")
		refinfile$gene<- sub("\\.[0-9]+$", "", refinfile$gene)
                name <- paste(poplist[i] %&% Nk[n],sep='')
                row<-poplist[i]%&% Nk[n]
                for(j in 1:length(poplist2)){
                        for(k in 1:length(Nk2)){
                                testinfile<-fread("/home/lauren/files_for_revisions_plosgen/meqtl_results/GEU/"%&% poplist2[j]%&%"_meQTL_results_"%&% Nk2[k] %&%"_Nk_3_PCs.txt")
                                testinfile$gene<- sub("\\.[0-9]+$", "", testinfile$gene)
                                col<-poplist2[j]%&% Nk2[k]
                                if(poplist[i] == poplist2[j]){
                                pi1 <- 1
                                 }else{
                                pop1 <- refinfile
                                pop2 <- testinfile
                                pop1fdr05 <- dplyr::filter(pop1, FDR < 0.05)
                                fdr<-dim(pop1fdr05)
                                pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
                                over<-dim(pop2tested)
                                pop2pvals <- pop2tested$pvalue.y
                                qobjpop2 <- qvalue(p = pop2pvals)
                                pi1 <- 1 - qobjpop2$pi0
                                }
                                pi12<-signif(pi1,4)
                                pi1_matrix[row,col] <- pi12
                                print(pi1_matrix)
                                write.table(pi1_matrix, file="/home/lauren/files_for_revisions_plosgen/pi1values_mesa_afhi_all_rest.txt",quote=F, sep = '\t')
                        }
                }
        }
}
write.table(pi1_matrix, file="/home/lauren/files_for_revisions_plosgen/pi1values_mesa_afhi_all_rest.txt",quote=F, sep = '\t')
