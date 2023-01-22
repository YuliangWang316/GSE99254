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

library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
epimarker<-read.table("D:/GSE139325_T_S_epigenetic_enzyme_for_volcano_plot_new.txt",sep = "\t",header = TRUE,row.names = 1)
for (i in c("wilcox","bimod","roc","t","negbinom","poisson","LR","DESeq2")) {
  if(i == "negbinom" | i == "poisson" | i == "DESeq2"){
    Tregmarker<-FindMarkers(NSCLC_count_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = i,slot = "counts")
  }else{
    Tregmarker<-FindMarkers(NSCLC_count_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = i)
  }
  
  genename<-rownames(Tregmarker)
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  Tregmarker<-Tregmarker[g$ENTREZID,]
  rownames(Tregmarker)<-g$SYMBOL
  Newname<-intersect(toupper(rownames(epimarker)),rownames(Tregmarker))
  Tregmarker_new<-Tregmarker[Newname,]
  write.table(Tregmarker_new,file = paste0("c:/Users/xjmik/Desktop/NSCLCmarker",i,".txt"),sep = "\t")
  remove(Tregmarker,Tregmarker_new,Newname,g,genename)
  
  
}
