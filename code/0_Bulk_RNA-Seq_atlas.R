# Load packages
library("ggplot2")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("DESeq2")

setwd('/hpc/group/goldsteinlab/R/Working_directory/bulk_seq_atlas')

#Code below adapted from: https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html


# First just check all available projects on TCGA
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])

#Set TCGA ID to the rownames
row.names(GDCprojects) <- GDCprojects$id

#Here now just looking at hepatocellular carcinoma
TCGAbiolinks:::getProjectSummary("TCGA-LIHC")

##########################################################################

# Can now begin pulling in 15 datasets from each cancer on TCGA

#Query TCGA specifically for bulk RNA-Seq file types
query_TCGA = GDCquery(
  project = "TCGA-DLBC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access= 'open')

#Quickly look at samples
lihc_res = getResults(query_TCGA)
summary(factor(lihc_res$sample_type))

#Randomly sample 15 sample ids from this
random_df <- lihc_res[sample(nrow(lihc_res), 15), ]
random_df$cases<-paste0("'",random_df$cases,"'")
barcodes <- gsub("\'","", random_df$cases)

#Now set query to just these cases
query_TCGA_lymphoma = GDCquery(
  project = "TCGA-DLBC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access= 'open',
  barcode = barcodes)

#Download files in dataset
GDCdownload(query = query_TCGA_lymphoma)

#Load RNA-Seq data into R
tcga_data_lymphoma = GDCprepare(query_TCGA_lymphoma)

dim(tcga_data_lymphoma)

colnames(colData(tcga_data_lymphoma))

table(tcga_data_lymphoma@colData$primary_diagnosis)
tcga_data_lymphoma@colData$primary_diagnosis <- 'Large B-cell Lymphoma'
#Standardize columns
tcga_data_lymphoma@colData <- tcga_data_lymphoma@colData[, c('primary_diagnosis', 'age_at_diagnosis', 'gender')]

tcga_data_lymphoma

##############################################

# Do again for another dataset

#Query TCGA specifically for bulk RNA-Seq file types
query_TCGA = GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access= 'open')

#Quickly look at samples
lihc_res = getResults(query_TCGA)
summary(factor(lihc_res$sample_type))

#Randomly sample 15 sample ids from this
random_df <- lihc_res[sample(nrow(lihc_res), 15), ]
random_df$cases<-paste0("'",random_df$cases,"'")
barcodes <- gsub("\'","", random_df$cases)

#Now set query to just these cases
query_TCGA = GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access= 'open',
  barcode = barcodes)

#Download files in dataset
GDCdownload(query = query_TCGA)

#Load RNA-Seq data into R
tcga_data_colon = GDCprepare(query_TCGA)

dim(tcga_data_colon)

colnames(colData(tcga_data_colon))

table(tcga_data_colon@colData$primary_diagnosis)
tcga_data_colon@colData$primary_diagnosis <- 'Colon adenocarcinoma'
#Standardize columns
tcga_data_colon@colData <- tcga_data_colon@colData[, c('primary_diagnosis', 'age_at_diagnosis', 'gender')]

tcga_data_colon

#########################################################################

# Continue doing this through all TCGA datasets
# Of course can also process via a for loop through all of them, but the method used here allows for verification at each dataset

# Once finised pulling in TCGA datasets, can concatenate

#Concatenate tcga data objects
tcga_all_combined <- cbind(tcga_data_adrenocortical, tcga_data_bladder,
                           tcga_data_breast_carcinoma, tcga_data_cervix,
                           tcga_data_cholangio, tcga_data_colon, tcga_data_endometrial_carcinoma,
                           tcga_data_GBM, tcga_data_head_and_neck, tcga_data_kidney_chromophobe,
                           tcga_data_kidney_renal_papillary, tcga_data_LGG, tcga_data_liver,
                           tcga_data_lung_adeno, tcga_data_lung_squamous, tcga_data_lymphoma,
                           tcga_data_meso, tcga_data_ovarian, tcga_data_pancreas,
                           tcga_data_Pheochrom_Paragang, tcga_data_prostate_adeno,
                           tcga_data_rectum_adeno, tcga_data_renal_clear_cell,
                           tcga_data_Sarcoma, tcga_data_skin_melanoma, tcga_data_stomach_adeno,
                           tcga_data_testicular, tcga_data_thymoma, tcga_data_thyroid_carcinoma,
                           tcga_data_uterine_sarcoma, tcga_data_uveal_melanoma)

dim(assay(tcga_combined))
head(assay(tcga_combined)[,1:10])
table(tcga_all_combined@colData$primary_diagnosis)

#Save this
saveRDS(object = tcga_all_combined,
        file = "all_tcga_data_combined.RDS",
        compress = FALSE)

#To read (if necessary)
tcga_data = readRDS(file = "all_tcga_data_combined.RDS")

########################################################################

# Now do this with TARGET, like so

#Query TARGET specifically for bulk RNA-Seq file types
TCGA_code <- "TARGET-NBL"

query_TCGA = GDCquery(
  project = TCGA_code,
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access= 'open')

#Quickly look at samples
lihc_res = getResults(query_TCGA)
summary(factor(lihc_res$sample_type))

#Randomly sample 15 sample ids from this
random_df <- lihc_res[sample(nrow(lihc_res), 15), ]
random_df$cases<-paste0("'",random_df$cases,"'")
barcodes <- gsub("\'","", random_df$cases)


#Now set query to just these cases
query_TCGA = GDCquery(
  project = TCGA_code,
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access= 'open',
  barcode = barcodes)

#Download files in dataset
GDCdownload(query = query_TCGA)

#Load RNA-Seq data into R
tcga_data = GDCprepare(query_TCGA)

dim(tcga_data)

tumor_name <- GDCprojects[TCGA_code, 'name']
table(tcga_data@colData$primary_diagnosis)
tcga_data@colData$primary_diagnosis <- tumor_name
#Standardize columns
tcga_data@colData <- tcga_data@colData[, c('primary_diagnosis', 'age_at_diagnosis', 'gender')]

Target_data_neuroblastoma <- tcga_data

# Like for TCGA, go through all desired TARGET datasets in this manner

##################################################

# Now can concatenate TCGA and TARGET datasets
tcga_target_combo <- cbind(tcga_all_combined, Target_data_neuroblastoma,
                           Target_data_osteosarcoma, Target_data_Rhabdoid,
                           Target_data_Wilms)

dim(assay(tcga_target_combo))
head(assay(tcga_target_combo)[,1:10])
table(tcga_target_combo@colData$primary_diagnosis)

#Save this
saveRDS(object = tcga_target_combo,
        file = "all_tcga_relevant_TARGET_combined_random_15.RDS",
        compress = FALSE)

#To read 
tcga_target_data = readRDS(file = "all_tcga_relevant_TARGET_combined_random_15.RDS")

############################################################

# Next, read in non TCGA and target datasets from published works

#First ONB ########################################################################################
#Read in ONB samples
ONB_df <- read.csv('Classe_et_al_19_ENB_raw_counts_output.csv')

#Remove NA gene rows
ONB_df <- ONB_df[complete.cases(ONB_df$Gene), ]

#Set rownames to genes
rownames(ONB_df) <- ONB_df$Gene

#Remove gene column
ONB_df$Gene <- NULL
dim(ONB_df)


#Now PNEC #####################################################################################
PNEC_df <- read.csv('PNEC_raw_counts.csv')

#Remove NA gene rows
PNEC_df <- PNEC_df[complete.cases(PNEC_df$Gene), ]

#Set rownames to genes
rownames(PNEC_df) <- PNEC_df$Gene

#Remove gene column
PNEC_df$Gene <- NULL
dim(PNEC_df)

#Rename columns with more clear names
num_columns <- ncol(PNEC_df)
num_rows_to_add <- 1
row_values <- paste0("Prostate_NE_", 1:num_columns)
# Create the new row with the desired values
new_row <- data.frame(t(row_values))
colnames(new_row) <- colnames(PNEC_df)
PNEC_df <- rbind(PNEC_df, new_row)
dim(PNEC_df)
colnames(PNEC_df) <- PNEC_df[25199, ]
PNEC_df <- PNEC_df[-25199, ]

#Now medulloblastoma read in #######################################################################
Medullo_df <- read.table('medulloblastoma_rnaseq.tsv', sep="\t", header=TRUE)
Medullo_df$sample <- Medullo_df$X
medullo_groups <- read.csv('medullo_subgroups.csv')

# Filter rownames in subgroups (since not all count matrices present)
medullo_groups <- medullo_groups[medullo_groups$sample %in% Medullo_df$sample, ]


medullo_full <- merge (Medullo_df, medullo_groups, by='sample')

rownames(medullo_full) <- medullo_full$sample
medullo_full$X <- NULL
medullo_full$sample <- NULL

#Move subgroup column to front
medullo_full <- medullo_full[, c("SUBGROUP", colnames(medullo_full)[colnames(medullo_full) != "SUBGROUP"])]

#Transpose
#medullo_full <- t(medullo_full)
#write.csv(medullo_full, file='medullo_raw_bulk_seq_with_subtypes.csv')


#Sort columns by subgroups (so all subgroups are together)
medullo_full <- medullo_full[apply(medullo)]

medullo_full$sample <- rownames(medullo_full)

# Now create proper sample names
unique_subtypes <- unique(medullo_full$SUBGROUP)
new_rows <- data.frame()


for (subtype in unique_subtypes) {
  # Subset dataframe for the current subtype
  subset_df <- medullo_full[medullo_full$SUBGROUP == subtype, ]
  
  # Generate sample numbers
  sample_numbers <- 1:nrow(subset_df)
  
  # Combine subtype, sample numbers, and original sample values
  combined_values <- paste("Medulloblastoma", subtype, sample_numbers, sep = "_")
  
  # Create a new row with combined values
  new_row <- data.frame(SUBTYPE = combined_values, sample = subset_df$sample)
  
  # Append the new row to the dataframe
  new_rows <- rbind(new_rows, new_row)
}

# Now combine dataframes
medullo_full_samples <- merge (medullo_full, new_rows, by='sample')
# move new rows to the front
medullo_full_samples <- medullo_full_samples[, c("SUBTYPE", colnames(medullo_full_samples)[colnames(medullo_full_samples) != "SUBTYPE"])]
rownames(medullo_full_samples) <- medullo_full_samples$sample
medullo_full_samples$sample <- NULL
colnames(medullo_full_samples)[colnames(medullo_full_samples)=='SUBTYPE'] <- 'Full_sample'

#Transpose
medullo_full_samples <- t(medullo_full_samples)

#Save this
# This is what should be read in for atlas processing
write.csv(medullo_full_samples, 'medulloblastoma_raw_data_read_for_atlas.csv')

medullo_df <- read.csv('medulloblastoma_raw_data_read_for_atlas.csv')

#Set rownames to genes
rownames(medullo_df) <- medullo_df$X

#Remove gene column
medullo_df$X <- NULL
dim(medullo_df)

#Set column names to full sample labeling
colnames(medullo_df) <- medullo_df['Full_sample', ]

#remove excess rows now
medullo_df <- medullo_df[-1, ]

#subset object into parts
medullo_SHH <- medullo_df[, medullo_df["SUBGROUP", ] == "SHH"]
medullo_WNT <- medullo_df[, medullo_df["SUBGROUP", ] == "WNT"]
medullo_Group3 <- medullo_df[, medullo_df["SUBGROUP", ] == "Group 3"]
medullo_Group4 <- medullo_df[, medullo_df["SUBGROUP", ] == "Group 4"]

#Now randomly subsetting samples per group
num_columns <- 5

selected_columns <- sample(ncol(medullo_WNT), num_columns)
medullo_WNT <- medullo_WNT[, selected_columns]
dim(medullo_WNT)

selected_columns <- sample(ncol(medullo_SHH), num_columns)
medullo_SHH <- medullo_SHH[, selected_columns]
dim(medullo_SHH)

selected_columns <- sample(ncol(medullo_Group3), num_columns)
medullo_Group3 <- medullo_Group3[, selected_columns]
dim(medullo_Group3)

selected_columns <- sample(ncol(medullo_Group4), num_columns)
medullo_Group4 <- medullo_Group4[, selected_columns]
dim(medullo_Group4)

# Can now aggregate these back into one df
subsampled_medullo_df <- cbind(medullo_SHH, medullo_WNT, medullo_Group3, medullo_Group4)
#remove subgroup row
subsampled_medullo_df <- subsampled_medullo_df[-1, ]
# Remove extra spaces at the beginning of values
row_names_medullo <- row.names(subsampled_medullo_df)
subsampled_medullo_df <- data.frame(lapply(subsampled_medullo_df, trimws, which = "l"), row.names=row_names_medullo)


# Now published SCLC data read in ###############################################################################
SCLC_df <- read.csv('sclc_bulk_published_only_subtype_labeled.csv')

# save subtype names for later
subtype_names <- SCLC_df[1,]
colnames(subtype_names)[colnames(subtype_names)=='X'] <- 'Gene'
SCLC_df <- SCLC_df[-1,]
colnames(SCLC_df)[colnames(SCLC_df)=='X'] <- 'Gene'


# Convert ensembl to gene names
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define the column name that contains Ensembl gene IDs
ensembl_col <- "Gene"

gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = SCLC_df[, ensembl_col],
                    mart = ensembl)

#Remove NA values
gene_names <- gene_names[complete.cases(gene_names),]

# Create a named vector of gene names
gene_names_vec <- setNames(gene_names$external_gene_name, gene_names$ensembl_gene_id)

# Assign the retrieved gene names to the corresponding rows
SCLC_df$Gene_name <- gene_names_vec[SCLC_df[, ensembl_col]]

# Convert blank values to NA
SCLC_df[SCLC_df == ""] <- NA

# Now remove all rows with NA values
SCLC_df <- SCLC_df[complete.cases(SCLC_df),]
rownames(SCLC_df) <- SCLC_df$Gene_name
SCLC_df$Gene_name <- NULL

# Now add back subtype info
SCLC_df_final <- rbind(SCLC_df, subtype_names)

SCLC_df_final[18013,]
SCLC_df_final[1,]
row_index <- 18013
new_rowname <- 'Full_Sample'
rownames(SCLC_df_final)[row_index] <- new_rowname
SCLC_df_final['Full_Sample',]
SCLC_df_final$Gene <- NULL
write.csv(SCLC_df_final, 'SCLC_published_raw_data_read_for_atlas.csv')

SCLC_df <- read.csv('SCLC_published_raw_data_read_for_atlas.csv')

#Set rownames to genes
rownames(SCLC_df) <- SCLC_df$X

#Remove gene column
SCLC_df$X <- NULL
dim(SCLC_df)

#Set column names to full sample labeling
colnames(SCLC_df) <- SCLC_df[18013, ]
SCLC_df <- SCLC_df[-18013, ]


# Filter all gene matrices based on genes only present in all datasets #####################################################

# Main TCGA and ONB
# Filter the gene matrix dataframe based on genes present in the genes_of_interest dataframe
filtered_counts_data <- counts_data[rownames(counts_data) %in% rownames(ONB_df), ]
dim(filtered_counts_data)
filtered_counts_data <- as.data.frame(filtered_counts_data)

# Main TCGA in and PNEC
filtered_counts_data <- filtered_counts_data[rownames(filtered_counts_data) %in% rownames(PNEC_df), ]
dim(filtered_counts_data)
filtered_counts_data <- as.data.frame(filtered_counts_data)

# Main TCGA in SCLC
filtered_counts_data <- filtered_counts_data[rownames(filtered_counts_data) %in% rownames(SCLC_df), ]
dim(filtered_counts_data)
filtered_counts_data <- as.data.frame(filtered_counts_data)

# Main TCGA in Medulloblastoma
filtered_counts_data <- filtered_counts_data[rownames(filtered_counts_data) %in% rownames(subsampled_medullo_df), ]
dim(filtered_counts_data)
filtered_counts_data <- as.data.frame(filtered_counts_data)

# ONB in main TCGA
#Now filter ONB_df based on genes present in filtereted_counts_data
filtered_ONB_df <- ONB_df[rownames(ONB_df) %in% rownames(filtered_counts_data), ]
dim(filtered_ONB_df)

# SCLC in Main TCGA
filtered_SCLC_df <- SCLC_df[rownames(SCLC_df) %in% rownames(filtered_counts_data), ]
dim(filtered_SCLC_df)

# PNEC in main TCGA
filtered_PNEC_df <- PNEC_df[rownames(PNEC_df) %in% rownames(filtered_counts_data), ]
dim(filtered_PNEC_df)

# Medulloblastoma in main TCGA
filtered_medullo_df <- subsampled_medullo_df[rownames(subsampled_medullo_df) %in% rownames(filtered_counts_data), ]
dim(filtered_medullo_df)

# Now ensuring sample names are appropriate ########################################################################

#Now set condition specific names as colnames in the tcga_target matrix
sample_order <- matrix_test$samples$primary_diagnosis
group_size <- 15
sample_order <- paste0(sample_order, "-", (seq_along(sample_order) -1) %% group_size +1)

colnames(filtered_counts_data) <- sample_order

#Need to rename rownames in PNEC df so they are PNEC-1, PNEC-2, etc
sample_numbers <- 1:21
new_column_names <- paste("Prostate_NE_", sample_numbers, sep = "")
# Rename the column names in the dataframe
colnames(PNEC_df) <- new_column_names


# Prep for concatenating all of the datasets
#First set rownames as a seperate column 
filtered_counts_data$Gene <- rownames(filtered_counts_data)
filtered_ONB_df$Gene <- rownames(filtered_ONB_df)
filtered_PNEC_df$Gene <- rownames(filtered_PNEC_df)
filtered_SCLC_df$Gene <- rownames(filtered_SCLC_df)
filtered_medullo_df$Gene <- rownames(filtered_medullo_df)
#filtered_mouseONB_df$Gene <- rownames(filtered_mouseONB_df)

library(dplyr)
#Merge with dplyr

#ONB/OE with TCGA
merged_df <- left_join(filtered_counts_data, filtered_ONB_df, by = "Gene")

#Full with PNEC
merged_df <- left_join(merged_df, filtered_PNEC_df, by = "Gene")

#Full with SCLC
merged_df <- left_join(merged_df, filtered_SCLC_df, by = "Gene")

#Full with medullo
merged_df <- left_join(merged_df, filtered_medullo_df, by = "Gene")

#Save gene names for later on
rownames(merged_df) <- merged_df$Gene

dim(merged_df)

#Remove NA gene rows (just in case there are any left)
# ideally should not be
merged_df <- merged_df[complete.cases(merged_df), ]

#Now eliminate gene column
merged_df$Gene <- NULL

#Now also need to update sample order
#Note this will not have number prefixes at end
sample_order <- matrix_test$samples$primary_diagnosis

new_samples_ONB <- rep('ONB', times=19)

new_samples_OE <- rep('Normal OE', times=3)

new_samples_PNEC <- rep('PNEC', times=21)

new_samples_SCLC <- rep('SCLC', times=16)

#THIS WILL CHANGE depending on subsampling
new_samples_medullo <- rep('Medulloblastoma', times=20)

updated_sample_order <- c(sample_order, new_samples_ONB, new_samples_OE, 
                          new_samples_PNEC, new_samples_SCLC,
                          new_samples_medullo)

ONB_OE_samples <- colnames(ONB_df)
PNEC_samples <- colnames(PNEC_df)
SCLC_samples <- colnames(SCLC_df)
Medullo_samples <- colnames(subsampled_medullo_df)
updated_sample_order_with_numbers <- c(sample_order, ONB_OE_samples, 
                                       PNEC_samples, SCLC_samples, 
                                       Medullo_samples)

sample_table <- data.frame(
  diagnosis=updated_sample_order
)

write.csv(merged_df, '604_samples_TCGA_TARGET_ONB_OE_SCLC_PNEC_Medullo_for_analysis.csv')
write.csv(updated_sample_order, 'updated_sample_order_604_samples.csv')

###########################################################################################

# Process counts dataframe

merged_df <- read.csv('604_samples_TCGA_TARGET_ONB_OE_SCLC_PNEC_Medullo_for_analysis.csv')

updated_sample_order_csv <- read.csv('updated_sample_order_604_samples.csv')
updated_sample_order <- updated_sample_order_csv$x

#Now run edgeR with updated matrix
merged_df_matrix <- as.matrix(merged_df)
class(merged_df_matrix) <- 'numeric'
dim(merged_df_matrix)

#Create edgeR object
d <- DGEList(counts=merged_df_matrix, group=factor(updated_sample_order))

# Can run edgeR pipeline as follows
#Filter data
apply(d$counts, 2, sum)

keep <- rowSums(edgeR::cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

#Reset library sizes post filtering
d$samples$lib.size <- colSums(d$counts)
d$samples

#Normalize (uses TMM normalization)
d <- calcNormFactors(d)


#Quick data exploration before DE
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)


#Estimate dispersion (helps determine variability of gene between samples)
#Here using GLM to estimate disperion
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="bin.spline")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.

d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

#Circle dendrogram
library(dendextend)
library(circlize)

counts_matrix <- edgeR::cpm(d, log = TRUE)
conditions <- factor(updated_sample_order)
dim(counts_matrix)

#save counts matrix before filtering on hvgs
write.csv(counts_matrix, '604sample_bulk_RNA_Seq_atlas.csv')

#Filter counts_matrix to HVGs
top_genes <- rownames(d2)[order(d2$tagwise.dispersion, decreasing = TRUE)][1:1000]
counts_matrix <- counts_matrix[top_genes, ]

counts_matrix <- t(counts_matrix)
dim(counts_matrix)
write.csv(counts_matrix, '1000hvgs_604sample_bulk_RNA_Seq_atlas.csv')

# Calculate the dissimilarity matrix using a suitable distance measure
dissimilarity_matrix <- dist(counts_matrix, method = "euclidean")

# Perform hierarchical clustering on the dissimilarity matrix
hclust_result <- hclust(dissimilarity_matrix, method = "complete")

# Convert the hierarchical clustering result to a dendrogram
dendro <- as.dendrogram(hclust_result)

#Create useful color palette
library(Polychrome)

# create your own color palette based on `seedcolors`

#How many unique conditions? (important for palette)
unique_conditions <- unique(conditions)
color_length <- length(unique_conditions)

P35 = createPalette(color_length,  c("#FF994E", "#0AAC00", "#438FFF"))
swatch(P35)

# Assign colors to the dendrogram branches based on conditions
dendro_colored <- color_branches(dendro, k = length(unique(conditions)), col = P35)

#Best way to save then is with export tab, somewhere around 51x48 inches for PDF will be required
dendro_circle <- circlize_dendrogram(dendro_colored, labels = TRUE, dend_track_height = 0.6, labels_track_height = .05)



