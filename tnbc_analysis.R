if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2"))

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

results <- getResults(query)
table(results$sample_type)

subtypes <- TCGAquery_subtype("BRCA")
table(subtypes$BRCA_Subtype_PAM50)

basal_patients <- subtypes$patient[subtypes$BRCA_Subtype_PAM50 %in% "Basal"]

basal_barcodes <- results$cases[substr(results$cases, 1, 12) %in% basal_patients &
                                  results$sample_type == "Primary Tumor"]

normal_barcodes <- results$cases[results$sample_type == "Solid Tissue Normal"]

length(basal_barcodes)
length(normal_barcodes)

basal_barcodes <- basal_barcodes[!duplicated(substr(basal_barcodes, 1, 12))]
length(basal_barcodes)

barcodes <- c(basal_barcodes, normal_barcodes)

query2 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = barcodes
)

GDCdownload(query2)

options(timeout = 1000)

GDCdownload(query2, method = "api", files.per.chunk = 10)

data <- GDCprepare(query2)
dim(data)

library(SummarizedExperiment)

data$group <- ifelse(data$sample_type == "Solid Tissue Normal", "Normal", "Basal")
table(data$group)

assayNames(data)

library(DESeq2)

data$group <- factor(data$group, levels = c("Normal", "Basal"))

dds <- DESeqDataSet(data, design = ~ group)

keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
nrow(dds)

dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "Basal", "Normal"))

res_df <- as.data.frame(res)
res_df$gene <- rowData(dds)$gene_name
res_df <- res_df[order(res_df$padj), ]

head(res_df[, c("gene", "log2FoldChange", "padj")], 10)

checks <- c("KRT5", "KRT17", "ESR1", "PIK3CA", "PTEN", "AKT1", "PIK3CB")
res_df[res_df$gene %in% checks, c("gene", "log2FoldChange", "padj")]

write.csv(res_df, "TNBC_vs_normal_DESeq2_results.csv")

up <- res_df[!is.na(res_df$padj) & res_df$padj < 0.01 & res_df$log2FoldChange > 2, ]
up$type <- rowData(dds)[rownames(up), "gene_type"]
up <- up[order(-up$log2FoldChange), ]

head(up[, c("gene", "type", "log2FoldChange", "padj")], 20)
table(up$type)

write.csv(res_df, "TNBC_vs_normal_DESeq2_results.csv")

#keep only well-expressed genes and save
up_strong <- up[up$baseMean > 100, ]
nrow(up_strong)
head(up_strong[, c("gene", "type", "baseMean", "log2FoldChange", "padj")], 20)

write.csv(res_df, "TNBC_vs_normal_DESeq2_results.csv")
write.csv(up_strong, "TNBC_candidate_genes.csv")

list.files(pattern = "csv")

getwd()
