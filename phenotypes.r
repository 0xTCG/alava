#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
print(args)

#%%
# if (!require("BiocManager", quietly=TRUE))
# install.packages("BiocManager")
# BiocManager::install(c(
#     "GSVA", "GSEABase", "GSVAdata", "zFPKM", "biomaRt", "org.Hs.eg.db",
#     "TCGAbiolinks", "TCGAretriever", "GWASinspector"), update=F, ask=F)

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

#%%
data(c2BroadSets)
gene_sets <- GeneSetCollection(c2BroadSets)

# Pathway data
kegg_paths <- unlist(as.list(org.Hs.egPATH2EG), use.names = FALSE)
# Load RNA-seq data
rnaseq_data <- read.table(args[1], header = TRUE, row.names = 1, sep = "\t")
# Only keep Entrez gene IDs
rownames(rnaseq_data) <- sapply(
  strsplit(rownames(rnaseq_data), split = "|", fixed = TRUE),
  function(x) (x[2])
)
# Optional: remove genes that are not expressed in any cell
# rnaseq_data <- rnaseq_data[rowSums(rnaseq_data > 0) > 0, ]

# Filter new RNA-seq file to include the genes for KEGG pathways only
rnaseq_kegg_genes <- rnaseq_data[rownames(rnaseq_data) %in% kegg_paths, ]
write.csv(
  cbind(gene_id = row.names(rnaseq_kegg_genes), rnaseq_kegg_genes),
  row.names = FALSE, quote = FALSE,
  file = "out/rnaseq-kegg-genes-only.csv"
)
cat("RNA-seq reading done; wrote out/rnaseq-kegg-genes-only.csv", "\n")

#%%
enrichment_gaussian <- gsva(
  ssgseaParam(data.matrix(rnaseq_kegg_genes), gene_sets), verbose = FALSE
)
write.table(enrichment_gaussian, "out/enrichment-ssgsea.txt")
enrichment_ssgsea <- gsva(
  ssgseaParam(data.matrix(rnaseq_kegg_genes), gene_sets, normalize = TRUE),
  verbose = FALSE
)
write.table(enrichment_ssgsea, "out/enrichment-ssgsea-norm.txt")
# Optional: remove genes that are not expressed in any cell
# rnaseq_data <- rnaseq_data[rowSums(data > 0) > 550, ]
enrichment_ssgsea <- gsva(
  ssgseaParam(data.matrix(rnaseq_data), gene_sets, normalize = TRUE),
  verbose = FALSE
)
write.table(enrichment_ssgsea, "out/enrichment-all-ssgsea-norm-full.txt")
cat("GSVA done; wrote out/enrichment-ssgsea.txt ",
    "out/enrichment-ssgsea-norm.txt and ",
    "out/enrichment-all-ssgsea-norm-full.txt", "\n")

#%%
data_enriched <- t(enrichment_ssgsea)
rownames(data_enriched) <- str_replace_all(rownames(data_enriched), "[.]", "-")
rownames(data_enriched) <- substr(rownames(data_enriched), 1, 12)
# This file contains list of KEGG HSAs and IDs:
# check https://rest.kegg.jp/list/pathway/hsa for details
filtered <- read.table(args[2], row.names = 1, sep = " ")
target_pathways <- data_enriched[, filtered$V2]
# write.table(target_pathways, "out/kegg-all-metabolic-phenotypes.txt")

#%%
# Original Python code:
# df = pd.read_csv("out/kegg-all-metabolic-phenotypes.txt", sep=" ")
# df = df.groupby([df.index]).mean().reset_index()
# df = df.rename(columns={'index': 'IID'})
# df.insert(loc=0, column='FID', value=df['IID'])
# df.to_csv('out/kegg-all-metabolic-phenotypes-corrected.txt', header=1, index=None, na_rep=0, sep=' ')

# Group by row index and take the mean
df <- as.data.frame(target_pathways)
df$index <- str_replace_all(rownames(df), "[.]", "-")
df <- aggregate(. ~ index, data = df, FUN = mean)
# Rename index to IID and insert FID=IID
names(df)[names(df) == "index"] <- "IID"
df <- cbind(FID = df$IID, df)

# Use hsa IDs
rename_map <- setNames(rownames(filtered), filtered$V2)
names(df) <- ifelse(
  names(df) %in% names(rename_map),
  rename_map[names(df)],
  names(df)
)

# Write output (space-separated, NA as 0, no row names)
write.table(
  df, file = "out/kegg-all-metabolic-phenotypes.txt",
  sep = " ", row.names = FALSE, quote = FALSE, na = "0"
)
cat("Grouping done; wrote out/kegg-all-metabolic-phenotypes.txt", "\n")

#%%
pathway_genes <- lapply(
  rownames(filtered),
  (function(x) {
    c(x,
      sort(unlist(
        as.list(org.Hs.egPATH2EG)[sub("hsa", "", x)],
        use.names = FALSE
      ))
    )
  })
)
con <- file("out/kegg-pathways-genes.txt", "w")
for (row_data in pathway_genes) {
  writeLines(paste(row_data, collapse = " "), con)
}
close(con)
cat("Pathway grouping done; wrote out/kegg-pathways-genes.txt", "\n")
