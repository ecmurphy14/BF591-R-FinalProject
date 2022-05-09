library(tidyverse)
library(GEOquery)
library(Biobase)
library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(GSEABase)
library(data.table)

#' Processing of GEO SOFT files to retrieve sample metadata
#' 
#' 
#' 

gse <- GEOquery::getGEO('GSE64810', GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))

filtered_metadata <- dplyr::select(metadata, 1:2, 6, 8:9, 54:65)

colnames(filtered_metadata) <- c('sample_name', 'geo_accession', 'type', 'source_name', 'organism',
                                 'age_at_death', 'age_at_onset', 'cag', 'diagnosis', 'duration', 'hv_cortical_score',
                                 'hv_striatal_score', 'mrna_seq_reads', 'pmi', 'rin', 'tissue', 'vonsattel_grade')

write.csv(filtered_metadata, 'sample_info.csv')



#' Running fgsea on DESeq2 results
#'

run_gsea <- function(deseq_res, min_size, max_size) {
  deseq_res <- dplyr::mutate(deseq_res, 
                                   genes= str_extract(genes, '.*\\.') 
                                   %>% str_replace('\\.', ''))

  get_hgnc <- AnnotationDbi::select(org.Hs.eg.db,
                                     key=deseq_res$genes, 
                                     columns="SYMBOL",
                                     keytype="ENSEMBL")
  
  get_hgnc <- as_tibble(get_hgnc)
  
  hgnc_res <- inner_join(deseq_res, get_hgnc, by=c("genes"="ENSEMBL"))
  
  rank_list <- hgnc_res %>% 
    dplyr::select(SYMBOL, log2FoldChange) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(log2FoldChange)) %>% 
    deframe()
  
  pathways <- gmtPathways("data/c2.cp.v7.5.1.symbols.gmt")

  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = rank_list,
                    minSize  = min_size,
                    maxSize  = max_size)
  
  return(fgseaRes)
}

deseqRes <- read_delim('data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.tsv',
                       delim = '\t', 
                       col_names = TRUE) %>%
  dplyr::rename('genes'='...1')

fgseaRes <- run_gsea(deseqRes, 
                     min_size = 15, 
                     max_size = 500)

fgseaRes <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

write_csv(fgseaRes, 'fgsea_results.csv')

