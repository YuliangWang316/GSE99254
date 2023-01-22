NSCLC_normlize_centered<-read.table("G:/GSE99254/GSE99254_NSCLC.TCell.S11769.norm.centered.txt",sep = "\t",header = TRUE,row.names = 1)
NSCLC_count<-read.table("G:/GSE99254/GSE99254_NSCLC.TCell.S12346.count.txt",sep = "\t",header = TRUE,row.names = 1)
NSCLC_TPM<-read.table("G:/GSE99254/GSE99254_NSCLC.TCell.S12346.TPM.txt",header = TRUE,sep = "\t",row.names = 1)

library(dplyr)
NSCLC_count_PTR<-select(NSCLC_count,starts_with("PTR"))
NSCLC_count_TTR<-select(NSCLC_count,starts_with("TTR"))
NSCLC_TPM_PTR<-select(NSCLC_TPM,starts_with("PTR"))
NSCLC_TPM_TTR<-select(NSCLC_TPM,starts_with("TTR"))
NSCLC_normlize_centered_PTR<-select(NSCLC_normlize_centered,starts_with("PTR"))
NSCLC_normlize_centered_TTR<-select(NSCLC_normlize_centered,starts_with("TTR"))
NSCLC_count_PTR_TTR<-select(NSCLC_count,starts_with("PTR"),starts_with("TTR"))
NSCLC_TPM_PTR_TTR<-select(NSCLC_TPM,starts_with("PTR"),starts_with("TTR"))
NSCLC_normlize_centered_PTR_TTR<-select(NSCLC_normlize_centered,starts_with("PTR"),starts_with("TTR"))

library(cowplot)
library(Seurat)

NSCLC_count_PTR_seurat <- CreateSeuratObject(counts = NSCLC_count_PTR, project = "IMMUNE_CTRL", min.cells = 5)
NSCLC_count_PTR_seurat$stim <- "CTRL"
NSCLC_count_PTR_seurat <- subset(NSCLC_count_PTR_seurat, subset = nFeature_RNA > 500)
NSCLC_count_PTR_seurat <- NormalizeData(NSCLC_count_PTR_seurat, verbose = FALSE)
NSCLC_count_PTR_seurat <- FindVariableFeatures(NSCLC_count_PTR_seurat, selection.method = "vst", nfeatures = 2000)

NSCLC_count_TTR_seurat <- CreateSeuratObject(counts = NSCLC_count_TTR, project = "IMMUNE_STIM", min.cells = 5)
NSCLC_count_TTR_seurat$stim <- "STIM"
NSCLC_count_TTR_seurat <- subset(NSCLC_count_TTR_seurat, subset = nFeature_RNA > 500)
NSCLC_count_TTR_seurat <- NormalizeData(NSCLC_count_TTR_seurat, verbose = FALSE)
NSCLC_count_TTR_seurat <- FindVariableFeatures(NSCLC_count_TTR_seurat, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(NSCLC_count_PTR_seurat, NSCLC_count_TTR_seurat), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim",label = TRUE)
Cluster0.markers <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
Cluster1.markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
Cluster2.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
Cluster3.markers <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
Cluster4.markers <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
write.table(Cluster0.markers,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster0.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster1.markers,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster1.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster2.markers,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster2.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster3.markers,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster3.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster4.markers,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster4.conservemarkers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
Cluster0.response <- FindMarkers(immune.combined, ident.1 = "0_STIM", ident.2 = "0_CTRL", verbose = FALSE)
Cluster1.response <- FindMarkers(immune.combined, ident.1 = "1_STIM", ident.2 = "1_CTRL", verbose = FALSE)
Cluster2.response <- FindMarkers(immune.combined, ident.1 = "2_STIM", ident.2 = "2_CTRL", verbose = FALSE)
Cluster3.response <- FindMarkers(immune.combined, ident.1 = "3_STIM", ident.2 = "3_CTRL", verbose = FALSE)
Cluster4.response <- FindMarkers(immune.combined, ident.1 = "4_STIM", ident.2 = "4_CTRL", verbose = FALSE)

write.table(Cluster0.response,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster0.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster1.response,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster1.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster2.response,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster2.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster3.response,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster3.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(Cluster4.response,file = "c:/Users/Administrator/Desktop/GSE99254/Cluster4.response.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

FeaturePlot(immune.combined, features = c("221037"), split.by = "stim", max.cutoff = 3, 
            label = TRUE)
VlnPlot(immune.combined, features = c("221037"), split.by = "stim", 
        pt.size = 0, combine = FALSE)

NSCLC_count_PTR_TTR_Seurat<-CreateSeuratObject(counts = NSCLC_count_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200)
NSCLC_count_PTR_TTR_Seurat <- subset(NSCLC_count_PTR_TTR_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
NSCLC_count_PTR_TTR_Seurat <- NormalizeData(NSCLC_count_PTR_TTR_Seurat)
NSCLC_count_PTR_TTR_Seurat <- FindVariableFeatures(NSCLC_count_PTR_TTR_Seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NSCLC_count_PTR_TTR_Seurat)
NSCLC_count_PTR_TTR_Seurat <- ScaleData(NSCLC_count_PTR_TTR_Seurat, features = all.genes)
NSCLC_count_PTR_TTR_Seurat <- RunPCA(NSCLC_count_PTR_TTR_Seurat, features = VariableFeatures(object = NSCLC_count_PTR_TTR_Seurat))
ElbowPlot(NSCLC_count_PTR_TTR_Seurat)
NSCLC_count_PTR_TTR_Seurat <- FindNeighbors(NSCLC_count_PTR_TTR_Seurat, dims = 1:10)
NSCLC_count_PTR_TTR_Seurat <- FindClusters(NSCLC_count_PTR_TTR_Seurat, resolution = 0.5)

PTR<-colnames(NSCLC_count_PTR)
TTR<-colnames(NSCLC_count_TTR)

Idents(NSCLC_count_PTR_TTR_Seurat, cells = PTR)<-'PTR'
Idents(NSCLC_count_PTR_TTR_Seurat, cells = TTR)<-'TTR'

NSCLC_count_PTR_TTR_Seurat <- RunUMAP(NSCLC_count_PTR_TTR_Seurat, dims = 1:10)
NSCLC_count_PTR_TTR_Seurat <- RunTSNE(NSCLC_count_PTR_TTR_Seurat, dims = 1:10)

DimPlot(NSCLC_count_PTR_TTR_Seurat, reduction = "umap")
DimPlot(NSCLC_count_PTR_TTR_Seurat, reduction = "tsne")

NSCLC.markers <- FindAllMarkers(NSCLC_count_PTR_TTR_Seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(NSCLC.markers,file = "c:/Users/Administrator/Desktop/GSE99254/NSCLC.markers.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

FeaturePlot(NSCLC_count_PTR_TTR_Seurat,features = "221037",label = TRUE)
VlnPlot(NSCLC_count_PTR_TTR_Seurat,features = "221037",pt.size = 0)

write.table(NSCLC_count_PTR_TTR,file = "c:/Users/Administrator/Desktop/GSE99254/NSCLC_count_PTR_TTR_.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(NSCLC_TPM_PTR_TTR,file = "c:/Users/Administrator/Desktop/GSE99254/NSCLC_TPM_PTR_TTR_.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(NSCLC_normlize_centered_PTR_TTR,file = "c:/Users/Administrator/Desktop/GSE99254/NSCLC_normlize_centered_PTR_TTR.txt",sep = "\t",row.names = TRUE,col.names = TRUE)

FeaturePlot(NSCLC_count_PTR_TTR_Seurat,features = "221037")

NSCLC_normlize_centered_PTR_TTR_Seurat<-CreateSeuratObject(counts = NSCLC_normlize_centered_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200)
NSCLC_normlize_centered_PTR_TTR_Seurat@assays[["RNA"]]@scale.data<-as.matrix(NSCLC_normlize_centered_PTR_TTR)


PTR<-colnames(NSCLC_normlize_centered_PTR)
TTR<-colnames(NSCLC_normlize_centered_TTR)

Idents(NSCLC_normlize_centered_PTR_TTR_Seurat, cells = PTR)<-'PTR'
Idents(NSCLC_normlize_centered_PTR_TTR_Seurat, cells = TTR)<-'TTR'
NSCLC_normlize_centered_PTR_TTR_Seurat <- RunPCA(NSCLC_normlize_centered_PTR_TTR_Seurat, features = rownames(NSCLC_normlize_centered_PTR_TTR))
NSCLC_normlize_centered_PTR_TTR_Seurat <- RunUMAP(NSCLC_normlize_centered_PTR_TTR_Seurat, dims = 1:10)
NSCLC_normlize_centered_PTR_TTR_Seurat <- RunTSNE(NSCLC_normlize_centered_PTR_TTR_Seurat, dims = 1:10)

DimPlot(NSCLC_normlize_centered_PTR_TTR_Seurat, reduction = "umap")
DimPlot(NSCLC_normlize_centered_PTR_TTR_Seurat, reduction = "tsne")

FeaturePlot(NSCLC_normlize_centered_PTR_TTR_Seurat,features = "221037")
VlnPlot(NSCLC_normlize_centered_PTR_TTR_Seurat,features = "221037",pt.size = 0)

