#!/usr/bin/env python

import subprocess
import time

CMD = 'qsub -v snp_annot_file={0},gene_annot_file={1},genotype_file={2},expression_file={3},covariates_file={4},chrom={5}, prefix2={6},' + \
    '-N {6}_model_chr{5}  train_afa_by_chr.sh'


gene_annot_file= '../../prepare_data/expression/gencode.v18.annotation.parsed.txt'
expression_file= '../../prepare_data/expression/AFA_MESA_Epi_GEX_data_sidno_Nk-20.txt'
covariates_file = '../../prepare_data/covariates/afa_PCs_sorted.txt'
prefix1= 'AFA_nested_cv'
prefix2 ='AFA_nested_cv_permuted'
for chr in range(1,23):
	genotype_file = '../../prepare_data/genotype/AFA_snp_chr1-22.chr'+str(chr)+'.txt'
	snp_annot_file = '../../prepare_data/genotype/AFA_annot_chr1-22.chr'+str(chr)+'.txt'
	cmd = CMD.format(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, str(chr), prefix2)
	print(cmd)
	subprocess.call(cmd, shell=True)
	time.sleep(2)

