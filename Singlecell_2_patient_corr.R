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

metadata<-read.table("D:/Jmjd1c_Treg_Tumor/GSE99254/GSE99254_family.txt",sep = "\t",header = TRUE)
PTR_metadata<-filter(metadata,sampleType == "PTR")
TTR_metadata<-filter(metadata,sampleType == "TTR")
PTR_TTR_metadata<-rbind(PTR_metadata,TTR_metadata)
rownames(PTR_TTR_metadata)<-PTR_TTR_metadata[,1]
NSCLC_count_PTR_TTR_Seurat<-CreateSeuratObject(counts = NSCLC_count_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200,meta.data = PTR_TTR_metadata)
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

B_TTR_seurat<-subset(NSCLC_count_PTR_TTR_Seurat,idents = "TTR")

Idents(B_TTR_seurat)<-B_TTR_seurat@meta.data$Patient

for (i in unique(B_TTR_seurat@meta.data[["Patient"]])) {
  assign(i,subset(B_TTR_seurat,idents = i))
  assign("a",slot(get(i),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@scale.data))))
  assign("c",select(b,one_of("221037","3458")))
  colnames(c)<-c("JMJD1C","IFNG")
  assign(paste("scaldata_",i,"_TTR",sep = ""),data.frame(sum(c$JMJD1C),sum(c$IFNG)))
}
scaldata<-rbind(scaldata_P0616A_TTR,scaldata_P0617_TTR,scaldata_P0706_TTR,scaldata_P0729_TTR,scaldata_P0913_TTR,scaldata_P1010_TTR,scaldata_P1011_TTR,scaldata_P1120_TTR,scaldata_P1202_TTR)
colnames(scaldata)<-c("JMJD1C","IFNG")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=IFNG)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

