library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

for(i in 1:22){
	list<-read.table("/home/lauren/files_for_revisions_plosgen/samples_HIS_233.txt", header = F)
	chr22<-fread("/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/genotypes/HIS_" %&% i %&% "_snp.txt", header = TRUE)
	rownames(chr22)<-chr22$id
	tchr22<-t(chr22)
	colnames(tchr22)<-tchr22[1,]
	dim(tchr22)
	chrfilt<-subset(tchr22, rownames(tchr22) %in% list$V1)
	dim(chrfilt)
	chrfin<-t(chrfilt)
	write.table(chrfin, file="/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/genotype_dw/HIS_" %&% i %&% "_snp.txt",quote=F, sep = '\t')
}

for(i in 1:22){
	list<-read.table("/home/lauren/files_for_revisions_plosgen/samples_CAU_233.txt", header = F)
	chr22<-fread("/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/genotypes/CAU_" %&% i %&% "_snp.txt", header = TRUE)
	rownames(chr22)<-chr22$id
	tchr22<-t(chr22)
	colnames(tchr22)<-tchr22[1,]
	dim(tchr22)
	chrfilt<-subset(tchr22, rownames(tchr22) %in% list$V1)
	dim(chrfilt)
	chrfin<-t(chrfilt)
	write.table(chrfin, file="/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/genotype_dw/CAU_" %&% i %&% "_snp.txt",quote=F, sep = '\t')
}
