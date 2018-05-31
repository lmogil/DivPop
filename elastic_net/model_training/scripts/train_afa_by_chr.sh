#!/bin/bash
###############################
# Resource Manager Directives #
###############################
#!/bin/bash
#PBS -N testjob2
#PBS -S /bin/bash
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -o ../joblogs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e ../joblogs/${PBS_JOBNAME}.e${PBS_JOBID}.err
 
#################
# Job Execution #
#################

cd  /home/lauren/PredictDB_Pipeline_GTEx_v7/model_training/scripts/

Rscript ./gtex_v7_nested_cv_elnet_n.R     

$snp_annot_file $gene_annot_file $genotype_file $expression_file $covariates_file $chrom $prefix2 
