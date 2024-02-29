#### Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(clustree)

#### load and normalize the data
scRNA <- readRDS("./scRNA.rds")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#### Identification of highly variable features
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

####Scaling the data
scale.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes, vars.to.regress = c("nCount_RNA", "percent.mt"))
                   
#### Perform linear dimensional reduction
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
                   
#### harmony batch correction
scRNA = RunHarmony(scRNA, "sampleinf",plot_convergence = TRUE)
                   
#### Cluster the cells
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims =1:30)
for (res in c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1)) {
  scRNA=FindClusters(scRNA,resolution = res)
}
                   
#### Run non-linear dimensional reduction (UMAP)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)
DimPlot(scRNA, reduction = "umap", group.by = "RNA_snn_res.0.3")
                   
#### Plot cluster tree
Clustering_tree = clustree(scRNA@meta.data, prefix = "RNA_snn_res.")

#### Identify differential expression gene
diff.wilcox = FindAllMarkers(scRNA,logfc.threshold = 0.25)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)

#### Plot the top10 gene heatmap
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = scRNA,features = top10$gene)

#### Plot the violin plot of specific gene
VlnPlot(scRNA,features = 'PTPRC')

#### Plot the feature plot of specific gene
FeaturePlot(object = scRNA,features = 'PTPRC')
