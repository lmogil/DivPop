library(dplyr)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b, sep='')

tissues <- c('ALL')
driver <- dbDriver("SQLite")
gene_annot <- read.table("../../prepare_data/expression/gencode.v18.annotation.parsed.txt", header = T, stringsAsFactors = F)

for (tiss in tissues) {
  print(tiss)
  # Extra table ----
  model_summaries <- read.table('/home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/' %&% tiss %&% '_nested_cv_chr1_model_summaries_prop.txt', header = T, stringsAsFactors = F)
  tiss_summary <- read.table('/home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/' %&% tiss %&% '_nested_cv_chr1_tiss_chr_summary_prop.txt', header = T, stringsAsFactors = F)
  
  n_samples <- tiss_summary$n_samples
  
  for (i in 1:22) {
    model_summaries <- rbind(model_summaries,
                             read.table('/home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/' %&% tiss %&% '_nested_cv_chr' %&% as.character(i) %&% '_model_summaries_prop.txt', header = T, stringsAsFactors = F))
    tiss_summary <- rbind(tiss_summary,
                             read.table('/home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/' %&% tiss %&% '_nested_cv_chr' %&% as.character(i) %&% '_tiss_chr_summary_prop.txt', header = T, stringsAsFactors = F))
  }
  
  model_summaries <- rename(model_summaries, gene = gene_id)

  conn <- dbConnect(drv = driver, dbname='../../dbs/' %&% tiss %&% '_imputed_prop_tw0.5.db')
  dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")
  
  # Weights Table -----
  weights <- read.table('/home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/' %&% tiss %&% '_nested_cv_chr1_weights_prop.txt', header = T, stringsAsFactors = F)
  for (i in 1:22) {
    weights <- rbind(weights,
                       read.table('/home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/' %&% tiss %&% '_nested_cv_chr' %&% as.character(i) %&% '_weights_prop.txt', header = T, stringsAsFactors = F))
  }
  weights <- rename(weights, gene = gene_id, weight = beta, eff_allele=alt)
  dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  
  # Sample_info Table ----
  sample_info <- data.frame(n_samples = n_samples, population = 'MESA', tissue = tiss)
  dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  
  # Construction Table ----
  construction <- tiss_summary %>%
                    select(chrom, cv_seed) %>%
                    rename(chromosome = chrom)
  dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
}

