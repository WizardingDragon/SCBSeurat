
data = read.table('data/MOUSE_BRAIN_DATASET_3_COUNTS.tsv',header=TRUE)
attach(data)

# # Initial -------------------------------------------------------------------
# 
# 
# # # Returns the number of genes present in Cell 1
# # sum(data$Cell.1!=0)
# 
# # # Create a new matrix that only gives us information about whether a gene is present or not
# # dataTest <- data!=0
# # geneCountPerCell <- apply(dataTest, 2, sum)
# 
# # Returns the number of gene per cell   
# geneCount <- apply(data!=0, 2, sum)
# 
# # Creates a fancy histogram of number of gene frequency
# library(ggplot2)
# geneCount_df<-data.frame(count=geneCount)
# ggplot(data=geneCount_df, aes(x=count))+geom_histogram()+
#   labs(title='Histogram of gene distribution',x='Number of genes per cell', y='Number of cells')
# # hist(geneCount, xlab='Number of genes per cell')
  
# Seurat -------------------------------------------------------------------
# Based on https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

library(Seurat)
library(dplyr)
# Create a Seurat object from the imported data
seurat_data<-CreateSeuratObject(counts=data,min.cells = 3,min.features = 200)

# It turns out that seurat_data@meta.data$nFeature_RNA already counts the data from geneCount_df
h=hist(seurat_data@meta.data$nFeature_RNA, main='Histogram of unique genes for each cell',
     xlab='Number of unique genes', ylab='Number of cells', breaks=40)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE, xlab = 'Number of unique genes', ylab = 'Percentage in the data')
# Returns the percentage of elements above 2500 unique genes
sum(h$density[15:length(h$density)], h$density[1])
# I will drop 11% of the genes by keeping a RNA count below 3000 as well as above 250
seurat_data <- subset(seurat_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
seurat_data <- subset(seurat_data, subset = percent.MT < 3)

# seurat_data@meta.data$nFeatures_RNA contains the total number of genes 
#in each cell (sum of number of copies too)

# Compute the mitochondrial gane
seurat_data[['percent.MT']]<-PercentageFeatureSet(seurat_data, pattern="^mt-")
# Creates a violin plot of the number of unique genes in each cell
VlnPlot(seurat_data, features=c("nFeature_RNA", "percent.MT"))

# Normalize the data before feature selection
seurat_data<-NormalizeData(seurat_data)

# Feature selection
seurat_data<-FindVariableFeatures(seurat_data, selection.method="vst", nFeatures=2000)

# Returns 10 most variable genes (i.e highly expressed in some cells, absent in others)
top10 <- head(VariableFeatures(seurat_data), 10)

# plot variable features with labels
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Data scaling
all.genes <- rownames(seurat_data)
seurat_data <- ScaleData(seurat_data, features = all.genes)

# Run PCA
seurat_data <- RunPCA(seurat_data, features = VariableFeatures(object = seurat_data))

# PCA visualization (?)
VizDimLoadings(seurat_data, dims = 1:2, reduction = "pca")

# PCA plot!
DimPlot(seurat_data, reduction = "pca")

# Heatmap
DimHeatmap(seurat_data, dims = 1, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat_data <- JackStraw(seurat_data, num.replicate = 100)
seurat_data <- ScoreJackStraw(seurat_data, dims = 1:15)

JackStrawPlot(seurat_data, dims = 1:15)

# Another way of plotting the same thing
# Ranks principle components based on the percentage of variance explained by each one
ElbowPlot(seurat_data)
# I choose to keep 15 PC

seurat_data <- FindNeighbors(seurat_data, dims = 1:10)
seurat_data <- FindClusters(seurat_data, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_data), 5)


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
# UMAP is neat! It's another way of doing PCA that grants more differentiated results, though
#less papers use UMAP than PCA, for the moment, since it's new
seurat_data <- RunUMAP(seurat_data, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat_data, reduction = "umap",cols=rainbow(9))
# SEE?! The clusters! They are nicely differentiated!

# Find cluster biomarkers -------------------------------------------------

# find all markers of cluster 1
cluster1.markers <- FindMarkers(seurat_data, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# Cluster makers ----------------------------------------------------------

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_data.markers <- FindAllMarkers(seurat_data, only.pos = TRUE, min.pct = 0.25, 
                                      logfc.threshold = 0.25)
# returns the gene that segregagtes the cells in clusters
seurat_data.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
# %>% required the dplyr library
gene_segregate<-seurat_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


# I chose Hexb as it was the first from cluster1.markers, same for Itm2b and Ugt8a, 
#but Itm2b returned blue for everything, thus also chose Ugt8a
FeaturePlot(seurat_data, features=c('Hexb','Itm2b','Ugt8a','C1qa'))

# I chose the genes from gene_segregate
FeaturePlot(seurat_data, features=c('Hexb','C1qa'))
FeaturePlot(seurat_data, features=c('Acta2','Myl9'))
FeaturePlot(seurat_data, features=c('Slc6a11','Fabp7'))
FeaturePlot(seurat_data, features=gene_segregate)


# Cluster identification --------------------------------------------------
# Gene origin given by mousebrain.org

# Returns the list of genes marking cluster 0 in a comma separated line easy to input on 
#mousebrain.org
# paste0(gene_segregate$gene[gene_segregate$cluster==0],collapse=",")

# Ctss,C1qa,P2ry12,C1qb,Selplg,C1qc,Tyrobp,Tmem119,Fcer1g,Hexb
# Cluster 0 : Microglia
# Acta2,Rgs5,Higd1b,Myl9,Igfbp7,Ndufa4l2,Tagln,Crip1,Cald1,Mgp
# Cluster 1: Vascular endothelial (or smooth muscle) cells
# Slc6a11,Gm3764,Bcan,Igfbp2,Ntsr2,Fabp7,Slc6a1,Slc1a3,Aldoc,Mt3
# Cluster 2: Astrocytes
# Ugt8a,Ermn,Fa2h,Etv1,Tubb4a,Mobp,S100b,Hmgcs1,Plp1,Fth1
# Cluster 3: Oligodendrocytes
# Mog,Apod,Cldn11,Stmn4,Klk6,Cnp,Fez1,Aplp1,Trf,Gsn
# Cluster 4: Mature oligodendrocytes
# Rab3c,Nptxr,Pgm2l1,Calb2,Nsg2,Celf4,Ntng1,Syp,Ndrg4,Syn2
# Cluster 5: Excitatory neurons
# Mpz,Prx,Sema3b,Ncmap,Egfl8,Fam178b,Spp1,Pmp22,Mbp,Cd9
# Cluster 6: Schwann Cells
# Eln,Crispld2,Tinagl1,Fos,Ppp1r15a,Tpm2,Egr1,Acta2,Btg2,Junb
# Cluster 7: Vascular Smooth Muscle cells (arterial)
# Ly6h,Chd3os,Sncb,Tcf7l2,Snhg11,Ahi1,Nap1l5,Pcp4,Meg3,6330403K07Rik
# Cluster 8: Excitatory Neurons
# Slc6a2,Dbh,Gal,Npy,Prph,Cyb561,S100a10,Stmn2,Tubb3,Sncg
# Cluster 9: Noradrenergic neurons

new.cluster.ids <- c("Microglia","Vascular cells, arterial","Astrocytes","Oligodendrocytes",
                     "Mature oligodendrocytes","Excitatory neurons","Schwann Cells",
                     "Vascular cells, arterial (2)","Excitatory neurons (2)", "Noradrenergic neurons")
names(new.cluster.ids) <- levels(seurat_data)
seurat_data <- RenameIdents(seurat_data, new.cluster.ids)
DimPlot(seurat_data, reduction = "umap", label = TRUE, pt.size = 0.5)
saveRDS(seurat_data, file = "output/seurat_final.rds")

# Returns information about the gene presence for markers gene of microglia
filter(seurat_data.markers, gene=="P2ry12" | gene=="Ccr6" | gene=="Tmem119" | 
         gene=="Casp8" | gene== "Ccl4" | gene=="Tnf")

library(ggthemes)
ggplot(data.frame(seurat_data@meta.data), aes(nFeature_RNA)) + 
  geom_histogram(bins=40, colour='steelblue3', fill='steelblue2')+theme_economist()+
  scale_x_continuous(name="Unique genes") + scale_y_continuous(name="Frequency")+
  geom_vline(aes(xintercept = mean(nFeature_RNA)),col='red',size=2)+
  geom_vline(aes(xintercept = median(nFeature_RNA)),col='orange',size=2)+
  ggtitle("Histogram of unique gene count per cell")
# Returns the mean number of genes for all cells
mean(seurat_data@meta.data$nFeature_RNA)
# Returns the number of cells for each cluster
summary(seurat_data$seurat_clusters)
# Returns the mean number of genes for each cell from cluster 0
mean(seurat_data@meta.data$nFeature_RNA[seurat_data$seurat_clusters==0])

VlnPlot(seurat_data, features=c("nFeature_RNA")) + geom_boxplot()+
  theme(legend.position = 'none')+ggtitle("Distribution of unique genes in each cluster")+
  ylab("Unique genes count")
