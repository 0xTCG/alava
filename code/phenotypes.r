#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
input_kegg <- args[2]
out_prefix <- args[3]
cat("Arguments:", input_data, input_kegg, out_prefix, "\n")

#%%
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
    "GSVA", "GSEABase", "GSVAdata", "zFPKM", "biomaRt", "org.Hs.eg.db",
    "TCGAbiolinks", "TCGAretriever", "GWASinspector"), update=F, ask=F)

suppressPackageStartupMessages({
  library(dplyr)
  library(GWASinspector)
  library(ggplot2)
  library(GSVA)
  library(GSEABase)
  library(GSVAdata)
  library(zFPKM)
  library(stringr)
  library(biomaRt)
  library(data.table)
  library(org.Hs.eg.db) # genome annotations for H.sapiens
})

#%% Part 1: enriching
data(c2BroadSets)

cat("Starting enrichment...\n")

# Load RNA-seq data
rnaseq_data <- read.table(input_data, header = TRUE, row.names = 1, sep = "\t")
# Only keep Entrez gene IDs
rownames(rnaseq_data) <- sapply(
  strsplit(rownames(rnaseq_data), split = "|", fixed = TRUE),
  function(x) (x[2])
)
# Optional: remove genes that are not expressed in any cell
# rnaseq_data <- rnaseq_data[rowSums(rnaseq_data > 0) > 0, ]
gene_sets <- GeneSetCollection(c2BroadSets)
enrichment_ssgsea <- gsva(
  ssgseaParam(data.matrix(rnaseq_data), gene_sets, normalize = TRUE),
  verbose = FALSE
)
cat("Enrichment (GSVA) done!\n")


#%% Part 2: filtering
data_enriched <- t(enrichment_ssgsea)
rownames(data_enriched) <- str_replace_all(rownames(data_enriched), "[.]", "-")
rownames(data_enriched) <- substr(rownames(data_enriched), 1, 12)
# This file contains list of KEGG HSAs and IDs:
# check https://rest.kegg.jp/list/pathway/hsa for details
filtered <- read.table(input_kegg, row.names = 1, sep = " ")
df <- as.data.frame(data_enriched[, filtered$V2])
# Use hsa IDs
rename_map <- setNames(rownames(filtered), filtered$V2)
names(df) <- ifelse(
  names(df) %in% names(rename_map),
  rename_map[names(df)],
  names(df)
)
rownames(df) <- str_replace_all(rownames(df), "[.]", "-")
dx <- aggregate(df, by = list(rownames(df)), FUN = mean)
dx <- cbind(FID = dx$Group.1, IID = dx$Group.1, dx)
dx$index <- NULL
dx$Group.1 <- NULL
out_path <- paste(out_prefix, "kegg-phenotypes.txt", sep = "/")
write.table(dx, out_path, sep = " ", row.names = FALSE, quote = FALSE, na = "0.0")
cat("Filtering done; wrote", out_path, "\n")


#%% Get [pathway: gene] IDs (one pathway per line)
pathway_genes <- lapply(
  rownames(filtered),
  (function(x) {
    c(x,
      filtered[x,],
      sort(unlist(as.list(org.Hs.egPATH2EG)[sub("hsa", "", x)], use.names = FALSE))
    )
  })
)
out_path <- paste(out_prefix, "kegg-pathways.txt", sep = "/")
con <- file(out_path, "w")
for (row_data in pathway_genes) {
  writeLines(paste(row_data, collapse = " "), con)
}
close(con)
cat("Pathway grouping done; wrote", out_path, "\n")

#%% Get Entrez IDs
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="hsapiens_gene_ensembl",
                        host="https://apr2020.archive.ensembl.org")
genes <- getBM(attributes=c("hgnc_symbol", "entrezgene_id"), mart=mart)
out_path <- paste(out_prefix, "entrez-ids.txt", sep = "/")
write.table(genes, out_path, row.names=FALSE, col.names=TRUE)



# Additional parts:
#
# Filter new RNA-seq file to include the genes for KEGG pathways only
# This part does filtering before enrichment; by default filtering is done after it.
#
# kegg_paths <- unlist(as.list(org.Hs.egPATH2EG), use.names = FALSE)
# rnaseq_kegg_genes <- rnaseq_data[rownames(rnaseq_data) %in% kegg_paths, ]
# write.csv(
#   cbind(gene_id = row.names(rnaseq_kegg_genes), rnaseq_kegg_genes),
#   row.names = FALSE, quote = FALSE,
#   file = "out/rnaseq-kegg-genes-only.csv"
# )
# cat("RNA-seq reading done; wrote out/rnaseq-kegg-genes-only.csv", "\n")
# enrichment_ssgsea <- gsva(
#   ssgseaParam(data.matrix(rnaseq_kegg_genes), gene_sets, normalize = TRUE),
#   verbose = FALSE
# )
# write.table(enrichment_ssgsea, "out/enrichment-ssgsea-norm.txt")
