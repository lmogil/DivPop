setwd("/home/lauren/files_for_revisions_plosgen/en_v7/model_training/scripts/")
source("gtex_v7_nested_cv_elnet_n.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
tiss <- argv[1]
chrom <- argv[2]

snp_annot_file <- "../../prepare_data/genotypes/" %&% tiss %&% "_" %&% chrom %&% "_annot.txt"
gene_annot_file <- "../../prepare_data/expression/gencode.v18.annotation.parsed.txt"
genotype_file <- "../../prepare_data/genotypes/" %&% tiss %&% "_" %&% chrom %&% "_snp.txt"
expression_file <- "../../prepare_data/expression/" %&% tiss %&% "_MESA_Epi_GEX_data_sidno_Nk-20.txt"
covariates_file <- "/home/lauren/files_for_revisions_plosgen/covariates/" %&% tiss %&% "pcs10cov.txt"
prefix <- tiss %&% "_nested_cv"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)


