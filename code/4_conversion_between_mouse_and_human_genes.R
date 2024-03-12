# Prepare workspace
library(Seurat)
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)
library(limma)
library(EnhancedVolcano)
library(gprofiler2)

# set wd
setwd('/hpc/group/goldsteinlab/R/Working_directory')


###############################################################
# Example: human to mouse
#Load in human list of genes from human SCLC
Pou2f3_up <- read.csv('Jchan_CCell_Pou2f3_up.csv')
Pou2f3_down <- read.csv('Jchan_CCell_Pou2f3_down.csv')


#Convert to mouse
mmus_Pou2f3_up = gorth(Pou2f3_up$genes, 
                       source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_Pou2f3_down = gorth(Pou2f3_down$genes, 
                         source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

#Now can write these as .csv so that you can load into Python to calculate module score
write.csv(mmus_Pou2f3_up, 'JChan_CCell_Pou2f3_up_mouse_versions.csv')
write.csv(mmus_Pou2f3_down, 'JChan_CCell_Pou2f3_down_mouse_versions.csv')

############################################################
# Example: mouse to human

# Convert RPM and RPMA gene set lists (DE from edgeR) to human orthologs

#Load in human list of genes from human SCLC
RPM_up <- read.csv('RPM_up_genes_mouse.csv')
RPMA_up <- read.csv('RPMA_up_genes_mouse.csv')


#Convert to mouse
human_RPM_up = gorth(RPM_up$genes, 
                     source_organism = "mmusculus", target_organism = "hsapiens")$ortholog_name
human_RPMA_up = gorth(RPMA_up$genes, 
                      source_organism = "mmusculus", target_organism = "hsapiens")$ortholog_name

#Now can write these as .csv so that you can load into Python to calculate module score
write.csv(human_RPM_up, 'mouse_edgeR_RPM_up_human_versions.csv')
write.csv(human_RPMA_up, 'mouse_edgeR_RPMA_up_human_versions.csv')








