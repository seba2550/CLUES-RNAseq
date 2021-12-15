# Directory
setwd("~/UCSF/Sirota Lab Rotation/SLE_RNA-seq_Counts")

# Libraries
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(AnnotationDbi)
library(ReactomePA)


# Load in the gene matrix (merged in Python)
gene_matrix <- read.csv("gene_matrix_count.csv", row.names = 1)

# Get the original column names so that I can get the Lanes for each sample
og_col_names <- colnames(gene_matrix)

# Get the lanes from the original sample names
og_col_names <- sapply(str_split(og_col_names, "Aligned", n = 2), `[`, 1)
og_col_names <- sapply(str_split(og_col_names, "_", n = 2), `[`, 2)
sample_lanes <- sapply(str_split(og_col_names, "_", n = 2), `[`, 2)

mod_col_names <- sapply(str_split(colnames(gene_matrix), "_",  n = 2), `[`, 1)
colnames(gene_matrix) <- mod_col_names

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

# Create cell-type specific subsets of the counts data
CD4_gene_matrix <- gene_matrix[,grep("CD4", colnames(gene_matrix))]
CD14_gene_matrix <- gene_matrix[,grep("CD14", colnames(gene_matrix))]
CD19_gene_matrix <- gene_matrix[,grep("CD19", colnames(gene_matrix))]
NK_gene_matrix <- gene_matrix[,grep("NK", colnames(gene_matrix))]

# Match the metadata for each cell-type specific counts matrix
CD4_metadata <- clues_time1_rnaseq_physact_correct_sample_names[clues_time1_rnaseq_physact_correct_sample_names$cell_type == "CD4",]
CD14_metadata <- clues_time1_rnaseq_physact_correct_sample_names[clues_time1_rnaseq_physact_correct_sample_names$cell_type == "CD14",]
CD19_metadata <- clues_time1_rnaseq_physact_correct_sample_names[clues_time1_rnaseq_physact_correct_sample_names$cell_type == "CD19",]
NK_metadata <- clues_time1_rnaseq_physact_correct_sample_names[clues_time1_rnaseq_physact_correct_sample_names$cell_type == "NK",]

# Cast clinical variables to strings
CD4_metadata$physact_classification <- as.character(CD4_metadata$physact_classification)
CD14_metadata$physact_classification <- as.character(CD14_metadata$physact_classification)
CD19_metadata$physact_classification <- as.character(CD19_metadata$physact_classification)
NK_metadata$physact_classification <- as.character(NK_metadata$physact_classification)

CD4_metadata$female <- as.character(CD4_metadata$female)
CD14_metadata$female <- as.character(CD14_metadata$female)
CD19_metadata$female <- as.character(CD19_metadata$female)
NK_metadata$female <- as.character(NK_metadata$female)

# Create four DESeq objects, one for each cell type
CD4_dds <- DESeqDataSetFromMatrix(countData = CD4_gene_matrix,
                                 colData = CD4_metadata,
                                 design = ~ sample_lanes + female + age + physact_classification)

CD14_dds <- DESeqDataSetFromMatrix(countData = CD14_gene_matrix,
                                  colData = CD14_metadata,
                                  design = ~sample_lanes + female + age + physact_classification)

CD19_dds <- DESeqDataSetFromMatrix(countData = CD19_gene_matrix,
                                  colData = CD19_metadata,
                                  design = ~sample_lanes + female + age + physact_classification)

NK_dds <- DESeqDataSetFromMatrix(countData = NK_gene_matrix,
                                  colData = NK_metadata,
                                  design = ~sample_lanes + female + age + physact_classification)

#################################################################################### Start of the DE analysis
#### CD4 Cells
# Filter out genes with less than 10 reads total
CD4_keep <- rowSums(counts(CD4_dds)) >= 10
CD4_dds <- CD4_dds[CD4_keep,]

# Calculate new size factors because many of our genes have 0s in them
CD4_dds <- estimateSizeFactors(CD4_dds, type = 'poscounts')


# Differential expression analysis
CD4_dds <- DESeq(CD4_dds)
CD4_res <- results(CD4_dds, contrast = c("physact_classification", "0", "1"))

# Transform the dds object 
CD4_vsd <- vst(CD4_dds)

# CD4-specific PCA plots
CD4_pca_physact_classification <- plotPCA(CD4_vsd, intgroup = "physact_classification")
CD4_pca_sample_lanes <- plotPCA(CD4_vsd, intgroup = "sample_lanes", returnData = T)

ggplot(CD4_pca_sample_lanes, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  scale_color_brewer(palette = "Accent")

# Summary statistics
summary(CD4_res)
sum(CD4_res$padj < 0.1, na.rm=TRUE)

# CD4 MA Plot
plotMA(CD4_res, ylim=c(-2,2))

# GO analysis
CD4_gene_list <- CD4_res$log2FoldChange[CD4_res$padj < 0.1]
CD4_gene_list <- na.omit(CD4_gene_list)
a <- which(CD4_res$padj < 0.1,arr.ind=T)
names(CD4_gene_list) <- rownames(CD4_res[a,])
CD4_gene_list <- sort(CD4_gene_list, decreasing = TRUE) # Genes must be sorted for gseGO to work

CD4_gse <- gseGO(geneList=CD4_gene_list, # Gene set enrichment analysis
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

ridgeplot(CD4_gse, label_format = 60) + labs(x = "enrichment distribution")
dotplot(CD4_gse, showCategory = 30, split = ".sign", label_format = 100) + facet_grid(.~.sign)  + ggtitle("Enriched GO Terms for CD4 Cells")


# PCA again, but removing batch effect
assay(CD4_vsd) <- limma::removeBatchEffect(assay(CD4_vsd), CD4_vsd$sample_lanes)

plotPCA(CD4_vsd, intgroup = "sample_lanes") + ggtitle("CD4 Samples After Removal of Lane Batch Effect")

CD4_count_table <- assay(CD4_vsd)

#### CD14 Cells
# Filter out genes with less than 10 reads total
CD14_keep <- rowSums(counts(CD14_dds)) >= 10
CD14_dds <- CD14_dds[CD14_keep,]

# Calculate new size factors because many of our genes have 0s in them
CD14_dds <- estimateSizeFactors(CD14_dds, type = 'poscounts')


# Differential expression analysis
CD14_dds <- DESeq(CD14_dds)
CD14_res <- results(CD14_dds, contrast = c("physact_classification", "0", "1"))

# Transform the dds object 
CD14_vsd <- vst(CD14_dds)

# CD14-specific PCA plots
CD14_pca_physact_classification <- plotPCA(CD14_vsd, intgroup = "physact_classification")
CD14_pca_sample_lanes <- plotPCA(CD14_vsd, intgroup = "sample_lanes", returnData = T)

ggplot(CD14_pca_sample_lanes, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  scale_color_brewer(palette = "Accent")

# Summary statistics
summary(CD14_res)
sum(CD14_res$padj < 0.1, na.rm=TRUE)

# CD14 MA Plot
plotMA(CD14_res, ylim=c(-2,2))


# GO analysis
CD14_gene_list <- CD14_res$log2FoldChange[CD14_res$padj < 0.1]
CD14_gene_list <- na.omit(CD14_gene_list)
CD14_a <- which(CD14_res$padj < 0.1,arr.ind=T)
names(CD14_gene_list) <- rownames(CD14_res[CD14_a,])
CD14_gene_list <- sort(CD14_gene_list, decreasing = TRUE) 

CD14_gse <- gseGO(geneList=CD14_gene_list, # Gene set enrichment analysis
                 ont ="ALL", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "fdr")



dotplot(CD14_gse, showCategory = 30, split = ".sign", label_format = 100) + facet_grid(.~.sign) + ggtitle("Enriched GO Terms for CD14 Cells")

# PCA again, but removing batch effect
assay(CD14_vsd) <- limma::removeBatchEffect(assay(CD14_vsd), CD14_vsd$sample_lanes)

plotPCA(CD14_vsd, intgroup = "sample_lanes") + ggtitle("CD14 Samples After Removal of Lane Batch Effect")

#### CD19 Cells
# Filter out genes with less than 10 reads total
CD19_keep <- rowSums(counts(CD19_dds)) >= 10
CD19_dds <- CD19_dds[CD19_keep,]

# Calculate new size factors because many of our genes have 0s in them
CD19_dds <- estimateSizeFactors(CD19_dds, type = 'poscounts')


# Differential expression analysis
CD19_dds <- DESeq(CD19_dds)
CD19_res <- results(CD19_dds, contrast = c("physact_classification", "0", "1"))

# Transform the dds object 
CD19_vsd <- vst(CD19_dds)

# CD19-specific PCA plots
CD19_pca_physact_classification <- plotPCA(CD19_vsd, intgroup = "physact_classification")
CD19_pca_sample_lanes <- plotPCA(CD19_vsd, intgroup = "sample_lanes", returnData = T)

ggplot(CD19_pca_sample_lanes, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  scale_color_brewer(palette = "Accent")

# Summary statistics
summary(CD19_res)
sum(CD19_res$padj < 0.1, na.rm=TRUE)

# CD19 MA Plot
plotMA(CD19_res, ylim=c(-2,2))





# Can't do gene set enrichment analysis because there are no DE genes. Sad/


# PCA again, but removing batch effect
assay(CD19_vsd) <- limma::removeBatchEffect(assay(CD19_vsd), CD19_vsd$sample_lanes)

plotPCA(CD19_vsd, intgroup = "sample_lanes") + ggtitle("CD19 Samples After Removal of Lane Batch Effect")


#### NK Cells
# Filter out genes with less than 10 reads total
NK_keep <- rowSums(counts(NK_dds)) >= 10
NK_dds <- NK_dds[NK_keep,]

# Calculate new size factors because many of our genes have 0s in them
NK_dds <- estimateSizeFactors(NK_dds, type = 'poscounts')


# Differential expression analysis
NK_dds <- DESeq(NK_dds)
NK_res <- results(NK_dds, contrast = c("physact_classification", "0", "1"))

# Transform the dds object 
NK_vsd <- vst(NK_dds)

# NK-specific PCA plots
NK_pca_physact_classification <- plotPCA(NK_vsd, intgroup = "physact_classification")
NK_pca_sample_lanes <- plotPCA(NK_vsd, intgroup = "sample_lanes", returnData = T)

ggplot(NK_pca_sample_lanes, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  scale_color_brewer(palette = "Accent")

# Summary statistics
summary(NK_res)
sum(NK_res$padj < 0.1, na.rm=TRUE)

# NK MA Plot
plotMA(NK_res, ylim=c(-2,2))


# GO analysis
NK_gene_list <- NK_res$log2FoldChange[NK_res$padj < 0.1]
NK_gene_list <- na.omit(NK_gene_list)
NK_a <- which(NK_res$padj < 0.1,arr.ind=T)
names(NK_gene_list) <- rownames(NK_res[NK_a,])
NK_gene_list <- sort(NK_gene_list, decreasing = TRUE) 

NK_gse <- gseGO(geneList=NK_gene_list, # Gene set enrichment analysis
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "fdr")



# No enriched terms for NK cells either, despite 156 DE genes

# PCA again, but removing batch effect
assay(NK_vsd) <- limma::removeBatchEffect(assay(NK_vsd), NK_vsd$sample_lanes)

plotPCA(NK_vsd, intgroup = "sample_lanes") + ggtitle("NK Samples After Removal of Lane Batch Effect")


############################################# Reactome and KEGG
# CD4 cells
CD4_ENSEMBL_Entrez <- mapIds(org.Hs.eg.db,
       keys = names(CD4_gene_list),
       column = "ENTREZID",
       keytype = "ENSEMBL",
       multiVals = "first")

CD4_Reactome_Enrichment <- enrichPathway(gene = CD4_ENSEMBL_Entrez, pvalueCutoff = 0.05, readable = T)
CD4_KEGG_Enrichment <- enrichKEGG(gene = CD4_ENSEMBL_Entrez, organism = "hsa", pvalueCutoff = 0.05)

barplot(CD4_KEGG_Enrichment)

browseKEGG(CD4_KEGG_Enrichment, 'hsa05171') # COVID-19 term
browseKEGG(CD4_KEGG_Enrichment, 'hsa05208') # Chemical carcinogenesis-ROS
browseKEGG(CD4_KEGG_Enrichment, 'hsa00190') # Oxidative phosphorylation
browseKEGG(CD4_KEGG_Enrichment, 'hsa04932') # Non-alcoholic fatty liver disease


# CD14 cells
CD14_ENSEMBL_Entrez <- mapIds(org.Hs.eg.db,
                             keys = names(CD14_gene_list),
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

CD14_Reactome_Enrichment <- enrichPathway(gene = CD14_ENSEMBL_Entrez, pvalueCutoff = 0.05, readable = T) # 0 enriched terms
CD14_KEGG_Enrichment <- enrichKEGG(gene = CD14_ENSEMBL_Entrez, organism = "hsa", pvalueCutoff = 0.05) # 0 enriched terms


### Not doing it for CD19 cells because there were no DE genes

# NK cells
NK_ENSEMBL_Entrez <- mapIds(org.Hs.eg.db,
                             keys = names(NK_gene_list),
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

NK_Reactome_Enrichment <- enrichPathway(gene = NK_ENSEMBL_Entrez, pvalueCutoff = 0.05, readable = T) # 0 enriched terms
NK_KEGG_Enrichment <- enrichKEGG(gene = NK_ENSEMBL_Entrez, organism = "hsa", pvalueCutoff = 0.05) # 0 enriched terms




####### Heatmaps
library(RColorBrewer)
library(pheatmap)
library(genefilter)
# CD4 cells
CD4_threshold <- CD4_res$padj < 0.1 & abs(CD4_res$log2FoldChange)> 0.5
CD4_res$threshold <- CD4_threshold     
CD4_select <- rownames(subset(CD4_res, threshold == T ))
CD4_df <- as.data.frame(colData(CD4_dds)[,c("female","age", "physact_classification")])

pheatmap(assay(CD4_vsd) [ CD4_select,], show_rownames = T, annotation_col = CD4_df, fontsize = 10, fontsize_row  = 10 ,scale = "row",width = 10, height = 10)

# CD14 cells
CD14_threshold <- CD14_res$padj < 0.1 & abs(CD14_res$log2FoldChange)> 0.5
CD14_res$threshold <- CD14_threshold     
CD14_select <- rownames(subset(CD14_res, threshold == T ))
CD14_df <- as.data.frame(colData(CD14_dds)[,c("female","age", "physact_classification")])

pheatmap(assay(CD14_vsd) [ CD14_select,], show_rownames = T, annotation_col = CD14_df, fontsize = 10, fontsize_row  = 10 ,scale = "row",width = 10, height = 10)

# No DE genes in CD19 cells so skip

# NK cells
NK_threshold <- NK_res$padj < 0.1 & abs(NK_res$log2FoldChange)> 0.5
NK_res$threshold <- NK_threshold     
NK_select <- rownames(subset(NK_res, threshold == T ))
NK_df <- as.data.frame(colData(NK_dds)[,c("female","age", "physact_classification")])

pheatmap(assay(NK_vsd) [ NK_select,], show_rownames = T, annotation_col = NK_df, fontsize = 10, fontsize_row  = 10 ,scale = "row", width = 20, height = 20)


####################### Overlap of DE genes across cell types (Venn diagrams)
venn.diagram(list(CD4_select, CD14_select, NK_select), 
             category.names = c("CD4 Genes", "CD14 Genes", "NK Genes"), 
             filename = "../Cell-Type_Specific_Venn_Diagram.png")
