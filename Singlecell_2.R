NSCLC_normlize_centered<-read.table("D:/Jmjd1c_Treg_Tumor/GSE99254/GSE99254_NSCLC.TCell.S11769.norm.centered.txt",sep = "\t",header = TRUE,row.names = 1)
NSCLC_count<-read.table("D:/Jmjd1c_Treg_Tumor/GSE99254/GSE99254_NSCLC.TCell.S12346.count.txt",sep = "\t",header = TRUE,row.names = 1)
NSCLC_TPM<-read.table("D:/Jmjd1c_Treg_Tumor/GSE99254/GSE99254_NSCLC.TCell.S12346.TPM.txt",header = TRUE,sep = "\t",row.names = 1)

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

IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_entrezID.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_entrezID.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

NSCLC_count_PTR_TTR_Seurat <- CellCycleScoring(NSCLC_count_PTR_TTR_Seurat, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
NSCLC_count_PTR_TTR_Seurat <- AddModuleScore(NSCLC_count_PTR_TTR_Seurat,features = IL2_stat5_list,name = "IL2_stat5")
NSCLC_count_PTR_TTR_Seurat <- AddModuleScore(NSCLC_count_PTR_TTR_Seurat,features = IL6_stat3_list,name = "IL6_stat3")

VlnPlot(NSCLC_count_PTR_TTR_Seurat,features = c("IL2_stat51","IL6_stat31","221037","8829","5133","3458"),pt.size = 0)
FeaturePlot(NSCLC_count_PTR_TTR_Seurat,features = c("IL2_stat51","IL6_stat31","221037","8829","5133","3458"))


NSCLC_count_PTR_TTR_marker<-FindMarkers(NSCLC_count_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = "DESeq2")
write.table(NSCLC_count_PTR_TTR_marker,file = "c:/Users/xjmik/Desktop/NSCLCdeseq2markers.txt",sep = "\t")
avg.NSCLC_count_PTR_TTR<-AverageExpression(NSCLC_count_PTR_TTR_Seurat)
write.table(avg.NSCLC_count_PTR_TTR,file = "G:/GSE99254/avg.NSCLC_count_PTR_TTR.txt",sep = "\t")
GSEA_TTR_PTR<-as.data.frame(NSCLC_count_PTR_TTR_Seurat@assays[["RNA"]]@scale.data)
write.table(GSEA_TTR_PTR,file = "G:/GSE99254/GSEA_TTR_PTR.txt",sep = "\t")

VlnPlot(NSCLC_count_PTR_TTR_Seurat,features = c("221037","5133","6776"),pt.size = 0)
