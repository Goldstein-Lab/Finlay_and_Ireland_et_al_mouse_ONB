# ScRNA seq RPMA + RPM GBC allograft analysis with CellTag Data # 
# Bringing in from SCANPY analysis #

# Load necessary packages
#load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(CellTagR)
  library(viridis)
  library(ggpubr)
  library(SeuratData)
  library(SeuratDisk)
  library(zellkonverter)
})



############## Seurat conversion and analysis using xenograft+organoid FA clustering from scanpy ##########

########################################################################################################
########################################################################################################
######### Now, bring into Seurat, from anndata object processed according to onb paper methods ########################################################################################################


# Read in adata object as SCE
# note: this was created in Jupyter Notebook 7, and saved at the end of section "scvi tumors only"
setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/Finlay_Ireland_ONB_FinalProducts")
adata<-readH5AD("ONB_primary_allografts_only_scvi.h5ad")


# Double check identities
table(adata$leiden_scVI_all_allograft_1.1)
table(adata$mouse_ident)


# Generate count tables and projection coordinates to compile a Seurat object
dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
length(rownames(adata))
assay(adata,"norm")

length(rowData(adata)$gene_ids)


#Convert SCE to seurat object
rpmaIV_seurat <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
rpmaIV_seurat


#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

rpmaIV_seurat
rpmaIV_seurat[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(rpmaIV_seurat)<-'norm'
table(rpmaIV_seurat@meta.data$leiden)


# Add embeddings of umap and fa projections to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(rpmaIV_seurat)
dim(test)
rpmaIV_seurat[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
rpmaIV_seurat[['umap']]@cell.embeddings

test2<-reducedDim(adata, "X_draw_graph_fa")
head(test2)
colnames(test2)<-c("FA_1","FA_2")
rownames(test2)<-colnames(rpmaIV_seurat)
dim(test2)
rpmaIV_seurat[['fa']] <- CreateDimReducObject(test2, key="FA_", assay = "RNA")
head(rpmaIV_seurat[['fa']]@cell.embeddings)

rpmaIV_seurat@meta.data

#Define color palette
colors<-c('#F39B7F99','brown1','gold','#a1c299','#2ca02c', '#0aa6d8', '#e7298a','#984ea3', '#A0522D','gray25', 'lightblue1' ,'gray80', 'orange1','#3C548899')


#Plot by leiden cluster (Figure 5E)
DimPlot(rpmaIV_seurat, cols=colors,group.by='leiden', reduction='fa',label=TRUE,label.size=10)&NoAxes()

# Export metadata to add CellTag annotations
write.csv(rpmaIV_seurat@meta.data, "seurat_fa_metadata_RPMAandRPMGBC_021724.csv")

#Import after adding CellTag annotations and add new metadata to seurat object
md<-read.csv("seurat_fa_metadata_RPMAandRPMGBC_021724_wCTandSource.csv")
md

ct_clone<-md$CellTag_Clone
ct_bin<-md$CellTag_Binary
source<-md$Source
table(source)

rpmaIV_seurat@meta.data$Source<-source
table(rpmaIV_seurat@meta.data$Source)
# RPM_GBC_Allograft            RPM_ONB RPMA_GBC_Allograft           RPMA_ONB 
# 4161               2822              19302                527 


rpmaIV_seurat@meta.data$CellTag_Clone<-ct_clone
table(rpmaIV_seurat@meta.data$CellTag_Clone)
# No_Clone   RPM_clone_1   RPM_clone_2   RPM_clone_3  RPMA_clone_1 RPMA_clone_14 RPMA_clone_17 RPMA_clone_19  RPMA_clone_2 RPMA_clone_20 
# 26154            55             3            11           273            31            27            34            11            47 
# RPMA_clone_23 RPMA_clone_33 RPMA_clone_34 RPMA_clone_36 RPMA_clone_37 RPMA_clone_40 RPMA_clone_42  RPMA_clone_7  RPMA_clone_8  RPMA_clone_9 
# 16            10            12            22            10            11            11            17            31            26 

rpmaIV_seurat@meta.data$CellTag_Binary<-ct_bin
table(rpmaIV_seurat@meta.data$CellTag_Binary)
# CellTag_Clone      Untagged 
# 658         26154


table(rpmaIV_seurat@meta.data$leiden)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 4675 3349 3312 3303 2573 2401 2128 1981 1298  596  556  479  130   31 


table(rpmaIV_seurat@meta.data$leiden)

#Re-order source meta-data and plot as in Fig 5E
rpmaIV_seurat@meta.data$Source<-factor(rpmaIV_seurat@meta.data$Source,levels=c("RPM_ONB","RPMA_ONB","RPM_GBC_Allograft","RPMA_GBC_Allograft"))
DimPlot(rpmaIV_seurat,group.by='Source',cols=c("darkorchid4","#ff7f0e","deepskyblue4","turquoise"),reduction='fa',shuffle=TRUE)
DimPlot(rpmaIV_seurat,group.by='Source',cols=c("darkorchid4","#ff7f0e","deepskyblue4","turquoise"),reduction='fa',split.by="Source",shuffle=TRUE)


#Save seurat object with updated metadata
saveRDS(rpmaIV_seurat,"Desktop/Medgenome_RPMAGBC_082023/GBC_adata_FAmap_alltumors_seurat.rds")
rpmaIV_seurat<-readRDS("Desktop/Medgenome_RPMAGBC_082023/GBC_adata_FAmap_alltumors_seurat.rds")


########################################################################
## Visualize Cell Cycle Analysis perforemd in scanpy under "phase" (Figure S6C) ##

rpmaIV_seurat@meta.data$phase<-factor(rpmaIV_seurat@meta.data$phase,levels=c("G1","S","G2M"))
DimPlot(rpmaIV_seurat, reduction = "fa", group.by = "phase",shuffle=TRUE,label=FALSE,pt.size=.4,cols=c("indianred3","green3","royalblue4"))+NoAxes()



### Look at distribution of sample type by leiden cluster as in Figure S6A ##

Idents(rpmaIV_seurat)<-'leiden'
Idents(rpmaIV_seurat)
x<-table(Idents(rpmaIV_seurat),rpmaIV_seurat@meta.data$Source)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=c("darkorchid4","#ff7f0e","deepskyblue4","turquoise"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=12), legend.title = element_text(size=20))


### Look at distribution of cell cycle phase by leiden cluster as in Figure S6C ##

Idents(rpmaIV_seurat)<-'leiden'
Idents(rpmaIV_seurat)

x<-table(Idents(rpmaIV_seurat),rpmaIV_seurat@meta.data$phase)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

# Generate stacked barplot (Figure S6C)
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")
p + scale_fill_manual(values=c("indianred3","green3","royalblue4"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


### Look at distribution of leiden cluster by sample as in Figure 5G for RPM and RPMA allografts ##

Idents(rpmaIV_seurat)<-'Source'
Idents(rpmaIV_seurat)

x<-table(Idents(rpmaIV_seurat),rpmaIV_seurat@meta.data$leiden)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

# Generate stacked barplot (Figure S6C)
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")
p + scale_fill_manual(values=colors)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


################### Generate dot plot of expression by cell type marker for Figure S6B ##########################

OE_cell_sigs<-read.csv("OE_Cell_types_jack.csv")
oe_cells<-OE_cell_sigs$Marker_genes
oe_cells

DotPlot(rpmaIV_seurat,features = oe_cells, group.by="leiden")+ theme(axis.text.x=element_text(angle=90, hjust=1,size=8), axis.text.y=element_text(hjust=1,size=10))+scale_color_viridis(option="viridis",direction= -1)


########### Identify DEGs per leiden cluster for supplemental table #############
################################################################################

table(adata_seurat@meta.data$leiden)
Idents(adata_seurat)<-'leiden'
Cluster.Markers <- FindAllMarkers(adata_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cluster.Markers
?FindAllMarkers
df<-as.data.frame(Cluster.Markers)
df
write.csv(df, "leiden_clustermarkers_030224.csv")


####################### Assign cell states/phenotypes as in Figure 5H based on DEG's per leiden cluster and known marker genes #########################

# Phenotype assignment based on leiden cluster
rpmaIV_seurat->adata_seurat
Idents(adata_seurat)
adata_seurat@meta.data$Phenotype<- ifelse(Idents(adata_seurat) %in% c("10","7"), "GBC-like",
                                          ifelse(Idents(adata_seurat) %in% c("12","0","4","3","11"), "Neuronal/INP-like",
                                                 ifelse(Idents(adata_seurat) %in% c("8"), "Non-neuronal (Epithelial/MV/Iono-like)", 
                                                        ifelse(Idents(adata_seurat) %in% c("2","1","9","13"), "Non-neuronal (Mesenchymal)", 
                                                               ifelse(Idents(adata_seurat) %in% c("6","5"), "Stem-like", "NA")))))

adata_seurat@meta.data$Phenotype<-factor(adata_seurat@meta.data$Phenotype, c("Stem-like", "GBC-like","Neuronal/INP-like","Non-neuronal (Epithelial/MV/Iono-like)","Non-neuronal (Mesenchymal)"))

table(adata_seurat@meta.data$Phenotype)
# Stem-like                               GBC-like                      Neuronal/INP-like 
# 4529                                   2537                                  11160 
# Non-neuronal (Epithelial/MV/Iono-like)             Non-neuronal (Mesenchymal) 
# 1298                                   7288 


DimPlot(adata_seurat, reduction = "fa", group.by = "Phenotype",shuffle=TRUE,label=FALSE,pt.size=.2,cols=c("deepskyblue4","orange3","#a1c299","indianred4","#F39B7F99"))+NoAxes()
states<-c("deepskyblue4","orange3","#a1c299","indianred4","#F39B7F99")



### Look at distribution of cell state by sample as in Figure 5H for RPM and RPMA allografts ##

Idents(adata_seurat)<-'Source'
Idents(adata_seurat)

x<-table(Idents(adata_seurat), adata_seurat@meta.data$Phenotype)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

# Generate stacked barplot (Figure S6C)
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")
p + scale_fill_manual(values=states)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


########### Identify DEGs per Phenotype/Cell State cluster for supplemental table  #############
########################################################################

Idents(adata_seurat)<-"Phenotype"
Idents(adata_seurat)
Cluster.Markers <- FindAllMarkers(adata_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cluster.Markers
?FindAllMarkers
test<-as.data.frame(Cluster.Markers)
test
write.csv(test, "Pheno_clustermarkers_022324.csv")

#################################################
######## Visualizing CellTag/clone data in FA map space #########
#################################################

adata_seurat$CellTag_Binary<-factor(adata_seurat$CellTag_Binary,levels=c("Untagged","CellTag_Clone"))
table(adata_seurat$CellTag_Binary)
# Untagged CellTag_Clone 
# 26154           658 

DimPlot(adata_seurat, group.by="CellTag_Binary", cols=c(rgb(0,.9,.7,.03),"deepskyblue4"), reduction='fa', shuffle=TRUE,pt.size=1)+NoAxes()


Idents(adata_seurat)<-'CellTag_Binary'
Idents(adata_seurat)

clones<-subset(adata_seurat,idents=c("CellTag_Clone"))
Idents(adata_seurat)

table(clones@meta.data$CellTag_Clone)


########## Generate stacked bar plot of Leiden cluster per clone as in Figure 5G #########

Idents(clones)<-'leiden'

x<-table(clones@meta.data$CellTag_Clone,Idents(clones))
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

#Re-order
proportions$Var1<-factor(proportions$Var1, levels=c("RPM_clone_1","RPM_clone_2","RPM_clone_3","RPMA_clone_42","RPMA_clone_23","RPMA_clone_14","RPMA_clone_34","RPMA_clone_36","RPMA_clone_17","RPMA_clone_2","RPMA_clone_20","RPMA_clone_19","RPMA_clone_7","RPMA_clone_9","RPMA_clone_40","RPMA_clone_8","RPMA_clone_1","RPMA_clone_33","RPMA_clone_37"))
proportions$Var1

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

#Re-order based on cluster re-order
ncolors<-c('#F39B7F99','brown1','#DBB40C', '#a1c299', '#2ca02c', '#0aa6d8','#e7298a','#984ea3', '#A0522D', 'gray','orange')

p+ scale_fill_manual(values=ncolors)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20, angle=90), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


## Generate stacked bar plot of Phenotype/Cell State per clone as in Figure 5H ###

Idents(clones)<-'Phenotype'
x<-table(clones@meta.data$CellTag_Clone,Idents(clones))
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions
proportions$Var1<-factor(proportions$Var1, levels=c("RPM_clone_1","RPM_clone_2","RPM_clone_3","RPMA_clone_42","RPMA_clone_23","RPMA_clone_14","RPMA_clone_34","RPMA_clone_36","RPMA_clone_17","RPMA_clone_2","RPMA_clone_20","RPMA_clone_19","RPMA_clone_7","RPMA_clone_9","RPMA_clone_40","RPMA_clone_8","RPMA_clone_1","RPMA_clone_33","RPMA_clone_37"))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

pheno<-c("deepskyblue4","orange3","#a1c299","indianred4","#F39B7F99")
p + scale_fill_manual(values=pheno)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(angle = 45)


###########################################################################################
####### Visualize individual clones for Figure S6D ##############

## Color and visualize clones
colors2<-c(rgb(0,0,0,'0.005'),'#F39B7F99','brown1','#DBB40C', 'darkturquoise', '#2ca02c', '#0aa6d8', '#e7298a',  'navyblue',  '#A0522D', 'darkorchid4', 'gold', '#984ea3', '#ff7f0e', '#3C548899', '#377eb8',  '#d62728',"forestgreen","blue", 'black',"orange")

table(clones@meta.data$CellTag_Clone)
Idents(clones)<-'CellTag_Clone'
Idents(adata_seurat)<-'CellTag_Clone'

# Plot all clones in one FA map #
DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE,cols=colors2,pt.size=1)+ggtitle("Robust Clones") & NoAxes()


# ID cell barcodes and highlight those in FA space for Figure S6D # 

r1<-subset(clones,idents=c('RPM_clone_1'))
r1<-rownames(r1@meta.data)
r1
a<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=r1,sizes.highlight=3, cols.highlight=c('#F39B7F99'))+ggtitle("RPM Clone 1 (n=55)")+
  theme(plot.title = element_text(face = "plain")) & NoLegend() & NoAxes()
a

r2<-subset(clones,idents=c('RPM_clone_2'))
r2<-rownames(r2@meta.data)
b<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=r2,sizes.highlight=3, cols.highlight=c("brown1"))+ggtitle("RPM Clone 2 (n=3)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


r3<-subset(clones,idents=c('RPM_clone_3'))
r3<-rownames(r3@meta.data)
c<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=r3,sizes.highlight=3, cols.highlight=c("#DBB40C"))+ggtitle("RPM Clone 3 (n=11)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra1<-subset(adata_seurat,idents=c('RPMA_clone_1'))
ra1<-rownames(ra1@meta.data)
ra1
d<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra1,sizes.highlight=3, cols.highlight=c("turquoise"))+ggtitle("RPMA Clone 1 (n=273)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra14<-subset(adata_seurat,idents=c('RPMA_clone_14'))
ra14<-rownames(ra14@meta.data)
e<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra14,sizes.highlight=3, cols.highlight=c("#2ca02c"))+ggtitle("RPMA Clone 14 (n=31)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra17<-subset(adata_seurat,idents=c('RPMA_clone_17'))
ra17<-rownames(ra17@meta.data)
f<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra17,sizes.highlight=3, cols.highlight=c("#0aa6d8"))+ggtitle("RPMA Clone 17 (n=27)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()



ra19<-subset(adata_seurat,idents=c('RPMA_clone_19'))
ra19<-rownames(ra19@meta.data)
h<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra19,sizes.highlight=3, cols.highlight=c("#e7298a"))+ggtitle("RPMA Clone 19 (n=34)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra2<-subset(adata_seurat,idents=c('RPMA_clone_2'))
ra2<-rownames(ra2@meta.data)
i<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra2,sizes.highlight=3, cols.highlight=c("navyblue"))+ggtitle("RPMA Clone 2 (n=11)")+
  theme(plot.title = element_text(face = "plain")) & NoLegend() & NoAxes()


ra20<-subset(adata_seurat,idents=c('RPMA_clone_20'))
ra20<-rownames(ra20@meta.data)
j<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra20,sizes.highlight=3, cols.highlight=c("darkorchid1"))+ggtitle("RPMA Clone 20 (n=47)")+
  theme(plot.title = element_text(face = "plain")) & NoLegend() & NoAxes()



ra23<-subset(adata_seurat,idents=c('RPMA_clone_23'))
ra23<-rownames(ra23@meta.data)
k<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra23,sizes.highlight=3, cols.highlight=c("gold"))+ggtitle("RPMA Clone 23 (n=16)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra33<-subset(adata_seurat,idents=c('RPMA_clone_33'))
ra33<-rownames(ra33@meta.data)
l<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra33,sizes.highlight=3, cols.highlight=c("#984ea3"))+ggtitle("RPMA Clone 33 (n=10)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra34<-subset(adata_seurat,idents=c('RPMA_clone_34'))
ra34<-rownames(ra34@meta.data)
m<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra34,sizes.highlight=3, cols.highlight=c("#ff7f0e"))+ggtitle("RPMA Clone 34 (n=12)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra36<-subset(adata_seurat,idents=c('RPMA_clone_36'))
ra36<-rownames(ra36@meta.data)
n<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra36,sizes.highlight=3, cols.highlight=c("#3C548899"))+ggtitle("RPMA Clone 36 (n=22)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()



ra37<-subset(adata_seurat,idents=c('RPMA_clone_37'))
ra37<-rownames(ra37@meta.data)
o<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra37,sizes.highlight=3, cols.highlight=c("#377eb8"))+ggtitle("RPMA Clone 37 (n=10)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra40<-subset(adata_seurat,idents=c('RPMA_clone_40'))
ra40<-rownames(ra40@meta.data)
p<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra40,sizes.highlight=3, cols.highlight=c("#d62728"))+ggtitle("RPMA Clone 40 (n=11)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()



ra42<-subset(adata_seurat,idents=c('RPMA_clone_42'))
ra42<-rownames(ra42@meta.data)
q<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra42,sizes.highlight=3, cols.highlight=c("forestgreen"))+ggtitle("RPMA Clone 42 (n=11)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra7<-subset(adata_seurat,idents=c('RPMA_clone_7'))
ra7<-rownames(ra7@meta.data)
r<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra7,sizes.highlight=3, cols.highlight=c("blue"))+ggtitle("RPMA Clone 7 (n=17)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra8<-subset(adata_seurat,idents=c('RPMA_clone_8'))
ra8<-rownames(ra8@meta.data)
s<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra8,sizes.highlight=3, cols.highlight=c("black"))+ggtitle("RPMA Clone 8 (n=31)") +
  theme(plot.title = element_text(face = "plain"))& NoLegend() & NoAxes()


ra9<-subset(adata_seurat,idents=c('RPMA_clone_9'))
ra9<-rownames(ra9@meta.data)
t<-DimPlot(adata_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=ra9,sizes.highlight=3, cols.highlight=c("orange"))+ggtitle("RPMA Clone 9 (n=26)")+
  theme(plot.title = element_text(face = "plain")) & NoLegend() & NoAxes()


colors2<-c(rgb(0,0,0,'0.005'),'#F39B7F99','brown1','#DBB40C', 'darkturquoise', '#2ca02c', '#0aa6d8', '#e7298a',  'navyblue',  '#A0522D', 'darkorchid4', 'gold', '#984ea3', '#ff7f0e', '#3C548899', '#377eb8',  '#d62728',"forestgreen","blue", 'black',"orange")


# Plot all #
grid.arrange(a,b,c,nrow=1)
grid.arrange(d, i, r, s, t, e, f, h, j, k, l, m,n, o, p, q, nrow=4)

saveRDS(adata_seurat,"021724_final_adata2seurat_RPMRPMAGBC.rds")
adata_seurat<-readRDS("021724_final_adata2seurat_RPMRPMAGBC.rds")