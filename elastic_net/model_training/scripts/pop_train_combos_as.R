setwd("/home/lauren/files_for_revisions_plosgen/en_v7/model_training/scripts/")
source("gtex_v7_nested_cv_elnet_combo_a0.05.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
tiss <- argv[1]
chrom <- argv[2]
peer<-argv[3]
covar<-argv[4]

snp_annot_file <- "../../prepare_data/genotypes/" %&% tiss %&% "_" %&% chrom %&% "_annot.txt"
gene_annot_file <- "../../prepare_data/expression/gencode.v18.annotation.parsed.txt"
genotype_file <- "../../prepare_data/genotypes/" %&% tiss %&% "_" %&% chrom %&% "_snp.txt"
expression_file <- "../../prepare_data/expression/meqtl_sorted_" %&% tiss %&% "_MESA_Epi_GEX_data_sidno_Nk-" %&% peer %&% ".txt"
covariates_file <- "../../prepare_data/covariates/" %&% tiss %&% "_" %&% covar %&% "_PCs.txt"
prefix <- tiss %&% "_nested_cv"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix)


