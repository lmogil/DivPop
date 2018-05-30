# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

library(MatrixEQTL)
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
pop <- args[1] #ex. MXL, YRI
Nk <- args[2] #number of PEER factors
covar <- args[3]
my.dir <- "/home/lauren/files_for_revisions_plosgen/"
cov.dir <- "/home/lauren/files_for_revisions_plosgen/covariates/"
out.dir <- "/home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/"

for(i in c(1:22)){

  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

  # Genotype file name
  SNP_file_name = (my.dir %&% "meqtl_files/" %&% pop %&% "chr" %&% i %&% "meqtl_w_expression.txt")
  #format:
  #rsid	26189	25542	24894
  #rs4468273	0.019	0.019	0.019
  snps_location_file_name = (my.dir %&% "meqtl_files/" %&% pop %&% "_pop_dosage_chr" %&% i %&% "_snpsloc.txt")


  # Gene expression file name
  expression_file_name = (my.dir %&% "expresion/old_exp/meqtl_sorted_" %&% pop %&% "_Epi_GEX_data_sidno_Nk-" %&% Nk %&% ".txt");
  gene_location_file_name = (my.dir %&% "expresion/old_exp/MESA_ens_loc_gencode18_exp.txt");


  # Covariates file name
  # Set to character() for no covariates
  if(covar == 0){
    covariates_file_name = character()
  }else{
    covariates_file_name = (cov.dir %&% pop %&% "pcs" %&% covar %&% "cov.txt")
  }

  # Output file name
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();

  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 1;
  pvOutputThreshold_tra = 0;

  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();

  # Distance for local gene-SNP pairs
  cisDist = 1e6;

  ## Load genotype data

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name); #does not take .gz

  ## Load gene expression data

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);

  ## Load covariates

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  ## Run the analysis (can I change this to fread?)
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

  unlink(output_file_name_tra);
  unlink(output_file_name_cis);

  ## Results:

  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  #cat('Detected local eQTLs:', '\n');
  #show(me$cis$eqtls)
  #cat('Detected distant eQTLs:', '\n');
  #show(me$trans$eqtls)

  ## Make the histogram plot of local and distant p-values
  png(filename = out.dir  %&% pop %&% "_chr" %&% i %&% "_Nk_" %&% Nk %&%"_PCs_" %&% covar %&% ".png")
  plot(me)
  dev.off()
  write.table(me$cis$eqtls, out.dir %&% pop %&% "_Nk_" %&% Nk %&% "_PFs_chr" %&% i %&% "pcs_" %&% covar %&% ".meqtl.cis." %&% date %&% ".txt", quote = FALSE, row.names = FALSE, sep = "\t")
}
