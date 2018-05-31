##### merge files for cis-eqtl summary stats
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

meqtl.dir="/home/lauren/files_for_revisions_plosgen/meqtl_results/GEU/"
genes.dir="/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/expression/"
dosage.dir="/home/lauren/files_for_revisions_plosgen/px_dosages_expression_ppl/"

afa_eqtl<- data.frame(fread(meqtl.dir %&% "AFA_meQTL_results_10_Nk_3_PCs.txt"))

gene_info<-data.frame(fread(genes.dir %&% "gencode.v18.annotation.parsed.txt"))
gene_info$gene_id<- sub("\\.[0-9]+$", "", gene_info$gene_id)
colnames(gene_info)[2]<-"gene"

afa_merge<-left_join(afa_eqtl,gene_info, by="gene")

afa_snp_info<-data.frame(fread(dosage.dir %&% "AFA_first6.txt"))
colnames(afa_snp_info)<-c("chr","snps","pos_snps","ref","alt","maf")


afa_merge_g<-left_join(afa_merge,afa_snp_info, by="snps")

afa_merge_g$chr.y<-NULL

colnames(afa_merge_g)[7]<-"chr"

afa_merge_g$maf<-NULL

write.table(afa_merge_g,"/home/lauren/files_for_revisions_plosgen/meqtl_results/AFA_cis_eqtl_summary_statistics.txt",quote=F, row.names=F, sep='\t')



his_eqtl<- data.frame(fread(meqtl.dir %&% "HIS_meQTL_results_10_Nk_3_PCs.txt"))

gene_info<-data.frame(fread(genes.dir %&% "gencode.v18.annotation.parsed.txt"))
gene_info$gene_id<- sub("\\.[0-9]+$", "", gene_info$gene_id)
colnames(gene_info)[2]<-"gene"

his_merge<-left_join(his_eqtl,gene_info, by="gene")

his_snp_info<-data.frame(fread(dosage.dir %&% "HIS_first6.txt"))
colnames(his_snp_info)<-c("chr","snps","pos_snps","ref","alt","maf")


his_merge_g<-left_join(his_merge,his_snp_info, by="snps")

his_merge_g$chr.y<-NULL

colnames(his_merge_g)[7]<-"chr"

his_merge_g$maf<-NULL

write.table(his_merge_g,"/home/lauren/files_for_revisions_plosgen/meqtl_results/HIS_cis_eqtl_summary_statistics.txt",quote=F, row.names=F, sep='\t')



cau_eqtl<- data.frame(fread(meqtl.dir %&% "CAU_meQTL_results_10_Nk_3_PCs.txt"))

gene_info<-data.frame(fread(genes.dir %&% "gencode.v18.annotation.parsed.txt"))
gene_info$gene_id<- sub("\\.[0-9]+$", "", gene_info$gene_id)
colnames(gene_info)[2]<-"gene"

cau_merge<-left_join(cau_eqtl,gene_info, by="gene")

cau_snp_info<-data.frame(fread(dosage.dir %&% "CAU_first6.txt"))
colnames(cau_snp_info)<-c("chr","snps","pos_snps","ref","alt","maf")


cau_merge_g<-left_join(cau_merge,cau_snp_info, by="snps")

cau_merge_g$chr.y<-NULL

colnames(cau_merge_g)[7]<-"chr"

cau_merge_g$maf<-NULL

write.table(cau_merge_g,"/home/lauren/files_for_revisions_plosgen/meqtl_results/CAU_cis_eqtl_summary_statistics.txt",quote=F, row.names=F, sep='\t')
