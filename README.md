# CLUES-RNAseq
A repo for storing the code I used to analyze clinical and RNA-seq data from the CLUES cohort as part of my rotation project in the Sirota Lab

Descriptions:
-The two csv files store the clinical data for various CLUES cohort patients at times 1 and 3, respectively.

-Correlation_Analysis.R holds the code for testing which clinical features might be correlated with physical activity outcomes. I tested this using logistic regression and also tested for enrichment using chi-square.

-Data_Wrangling_Clin_and_RNAseq.R has code for the correlation analysis I was doing but using RAPA scores. These were inconsistent at time 1 and therefore we dropped this approach. The code is still shared here for documentation purposes.

-FeatureCounts_Merge.ipynb is a Jupyter Notebook used only for merging the 400 featureCounts outputs. I tried performing this in R with an inner join but it was taking an abnormally long time. Python on the other hand had a function from the Bioinfokit library that can merge these outputs into a single table in a much more efficient and quick manner.

-RNA-seq_Initial_Analysis.R holds the sample-wide analysis I did (not stratifying by cell-type), including the data wrangling for joining the clinical data and gene counts matrices, and then performing differential expression analysis. 

-RNA-seq_Cell_Type_Analysis.R has code for doing the same things in the initial analysis script but stratifying by cell-type. Note: The enrichment stuff only works completely in the CD4 cells.

-RNAseq_samples.txt has the sample IDs to identify which patients we have both clinical data AND RNA-seq data for.
