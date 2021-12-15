# Directory
setwd("~/UCSF/Sirota Lab Rotation/SLE_RNA-seq_Counts")

# Libraries
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)


# Load in the gene matrix (merged in Python)
gene_matrix <- read.csv("gene_matrix_count.csv", row.names = 1)

# Get the original column names so that I can get the Lanes for each sample
og_col_names <- colnames(gene_matrix)

# Get the lanes from the original sample names
og_col_names <- sapply(str_split(og_col_names, "Aligned", n = 2), `[`, 1)
og_col_names <- sapply(str_split(og_col_names, "_", n = 2), `[`, 2)
sample_lanes <- sapply(str_split(og_col_names, "_", n = 2), `[`, 2) # This is identical on purpose. Now we can get the lanes.

# Get a vector with modified column names, so that we can match them to the clinical data. 
# Then place those new column names in the gene matrix
mod_col_names <- sapply(str_split(colnames(gene_matrix), "_",  n = 2), `[`, 1)
colnames(gene_matrix) <- mod_col_names # They're in the same order, so all is well


# Get only the sample IDs from the vector we made earlier
tmp <- str_split(mod_col_names, "[.]", n = 2)
sample_IDs_sorted <- sapply(tmp, tail, 1)
cell_type <- sapply(tmp, head, 1)
rm(tmp)

# Create a tmp df to join the clinical data to the modified sample names (this is so that they match the gene matrix columns)
tmp_df <- tibble(mod_col_names, sample_IDs_sorted, cell_type, sample_lanes)

# Load in the clinical data
clues_time1_rnaseq_physact <- read.csv("../clues_time_rnaseq_physact.csv", row.names = 1)
clues_time1_rnaseq_physact$subjectid <- as.character(clues_time1_rnaseq_physact$subjectid)

# Join the tmp df and the clinical data
clues_time1_rnaseq_physact_correct_sample_names <-left_join(tmp_df, clues_time1_rnaseq_physact, 
                                                            by = c("sample_IDs_sorted" = "subjectid")) # Manually checked and the order of the rows lines up with the columns in the gene counts matrix

# Find individuals for which we lack physical activity info (NAs in physact_classification column)
ppl_without_phys_act <- clues_time1_rnaseq_physact_correct_sample_names[is.na(clues_time1_rnaseq_physact_correct_sample_names$physact_classification),] # 80 because there are 20 individuals and 4 datasets for each individual
ppl_without_phys_act_vec <- ppl_without_phys_act$mod_col_names

# Drop the respective columns from the counts matrix and from the metadata/coldata
gene_matrix <- gene_matrix[,!(names(gene_matrix) %in% ppl_without_phys_act_vec)]
clues_time1_rnaseq_physact_correct_sample_names <- clues_time1_rnaseq_physact_correct_sample_names[!clues_time1_rnaseq_physact_correct_sample_names$mod_col_names %in% ppl_without_phys_act_vec,]

# Cast some of the clinical variables to strings to avoid downstream trouble
clues_time1_rnaseq_physact_correct_sample_names$physact_classification <- as.character(clues_time1_rnaseq_physact_correct_sample_names$physact_classification)
clues_time1_rnaseq_physact_correct_sample_names$female <- as.character(clues_time1_rnaseq_physact_correct_sample_names$female)

# Create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = gene_matrix,
                              colData = clues_time1_rnaseq_physact_correct_sample_names,
                              design = ~sample_lanes  + female + age + cell_type + physact_classification)

# Double check to make sure we have our labels all lined up correctly
gene_matrix_cols <- colnames(gene_matrix)
identical(gene_matrix_cols, clues_time1_rnaseq_physact_correct_sample_names$mod_col_names) # Returns True, so we've got this right



####################### Start of the DESeq2 analysis/QC steps

# Filter out genes with less than 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Calculate new size factors because many of our genes have 0s in them
dds <- estimateSizeFactors(dds, type = 'poscounts')


# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("physact_classification", "0", "1"))


# Transform the dds object for different analyses (clustering, heat maps, etc)
vsd <- vst(dds)

# Plot the PCA, coloring by various features of value
pca_cell_type <- plotPCA(vsd, intgroup = "cell_type") + ggtitle("Sample-wide PCA, after including cell type as covariate")
pca_cell_type_df <- plotPCA(vsd, intgroup = "cell_type", returnData = T)
pca_physact_classification <- plotPCA(vsd, intgroup = "physact_classification") + ggtitle("Sample-wide PCA, after including cell type as covariate")
pca_sledai_score <- plotPCA(vsd, intgroup = "sledaiscore")
pca_sample_lanes <- plotPCA(vsd, intgroup = "sample_lanes")

# Plot the PCA, but removing the batch effect
assay(vsd) <- removeBatchEffect(assay(vsd), vsd$sample_lanes)

plotPCA(vsd, intgroup = "sample_lanes") + ggtitle("Sample-wide PCA After Removal of Lane Batch Effect")


# Get some summarized statistics for the results object
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

# MA plot 
plotMA(res, ylim=c(-2,2))

# Volcano plot
par(mfrow=c(1,1))# Reset par

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Sample-wide volcano plot (Covariates: Age, Sex, Sample Lanes, Cell Type)", xlim=c(-3,3))) # Make a basic volcano plot


with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue")) # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Extract the genes that have significant differential expression
resSig <- subset(res, padj < 0.1)
resSig_df <- as.data.frame(resSig)

# Get the gene IDs
resSig_df <- rownames_to_column(resSig_df, var = "GeneID")

# Perform GO enrichment
ego_bp <- enrichGO(gene         = resSig_df$GeneID, # 16 enriched terms
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)
head(ego_bp) # Look at the terms

ego_mf <- enrichGO(gene         = resSig_df$GeneID, # 2 enriched terms
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
head(ego_mf)

ego_cc <- enrichGO(gene         = resSig_df$GeneID, # 3 enriched terms
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
head(ego_cc)

# Try out gene set enrichment and see if the results change
samplewide_gene_list <- res$log2FoldChange[res$padj < 0.1]
samplewide_gene_list <- na.omit(samplewide_gene_list)
res_a <- which(res$padj < 0.1,arr.ind=T)
names(samplewide_gene_list) <- rownames(res[res_a,])
samplewide_gene_list <- sort(samplewide_gene_list, decreasing = TRUE) # Genes must be sorted for gseGO to work

samplewide_gse <- gseGO(geneList=samplewide_gene_list, 
                 ont ="ALL", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "fdr")

# No enriched terms

