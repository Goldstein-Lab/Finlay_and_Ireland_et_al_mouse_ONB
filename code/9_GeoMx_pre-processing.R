#GeoMx data analysis
#code adapted from https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

#Library
library(ggplot2)
library(GeomxTools)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(knitr)
library(dplyr)
library(ggforce)
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)
library(umap)
library(Rtsne)
library(pheatmap)
library(ggrepel)
library(Seurat)
library(dittoSeq)

#Setwd to goldsteinlab folder 
setwd('/hpc/group/goldsteinlab/R/Working_directory/Nanostring')

#Set file path
datadir <- file.path('/hpc/group/goldsteinlab/R/Working_directory/Nanostring/outs')

#List files in each directory
DCCFiles <- dir(file.path(datadir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

PKCFiles <- dir(file.path(datadir, "pkcs"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)

SampleAnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)

#Load data
demoData <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = 'Annotation template',
                         phenoDataDccColName = "dcc_file",
                         experimentDataColNames = c("ROI_ID_2"))

#Check to make sure proper probe annotations are being used
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

#########################################################################################################
#QC

#Shift counts to 1
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

#Make flags for data
#These are the params used in this study
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       minnuclei = 20,         # Minimum # of nuclei estimated (100)
       minarea = 1000)         # Minimum segment area (5000)
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

#Visualize segment QC
col_by <- "segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)

kable(QC_Summary, caption = "QC Summary Table for each Segment")

# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

#Probe QC
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

#Create Gene-level count data
# Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))

# collapse to targets
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)
exprs(target_demoData)[1:5, 1:2]

#Calculate Limit of Quantification
# Define LOQ SD threshold and minimum value
cutoff <- 1.25
minLOQ <- 1.25

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ

##################################################################################################
#Filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

################################################

#First, find segments with very low signal
# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Tumor)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type") +
  scale_fill_brewer(palette="Set3")

#Filter out segments with low detection
#Here using 15%
target_demoData <-
  target_demoData[, pData(target_demoData)$GeneDetectionRate >= .15]

dim(target_demoData)

###################################################################################
# Perform Q3 normalization for initial unsupervised clustering
# this was just done for the first go and generation of UMAP; for all downstream analyses, quantile normalization was performed

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "Scan_name"

Stat_data <- 
  data.frame(row.names = colnames(exprs(target_demoData)),
             Segment = colnames(exprs(target_demoData)),
             Annotation = pData(target_demoData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_demoData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_demoData)[neg_probes, ])

Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")


plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_demoData <- normalize(target_demoData ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

########################################################################################################
# Initial Unsupervised clustering

#Plot UMAP
# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),  
       config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = marker, shape=Grade_hi_low)) + theme_bw() +
  geom_point(size = 3) + theme(text=element_text(family="Arial", size=20)) + scale_color_brewer(palette="Set1")


#######################################################
# for workability, convert raw Nanostring GeoMx object to Seurat
seurat_demoData <- as.Seurat(target_demoData, ident = "Tumor",
                             normData = "q_norm", forceRaw = TRUE)

# now convert Seurat object to .h5ad for usage in Python
library(SeuratDisk)

# First write h5seurat
SaveH5Seurat(seurat_demoData, filename = "GeoMx_82_ROI_3.h5Seurat")

# then convert to .h5ad
Convert("GeoMx_82_ROI_3.h5Seurat", dest = "h5ad")

# also save metadata
meta_data <- seurat_demoData@meta.data
write.csv(meta_data, 'Geomx_metadata.csv')









