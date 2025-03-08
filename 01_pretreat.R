#####################
# doublet removal
library(Seurat)
library(DoubletFinder)

obj_filterdoublet <- function(input.mat, input.level = "Weak"){
  obj_tmp <- CreateSeuratObject(input.mat, min.features = 100)
  obj_tmp <- NormalizeData(obj_tmp, verbose = F)
  obj_tmp <- FindVariableFeatures(obj_tmp, selection.method = "vst", nfeatures = 1000, verbose = F)
  obj_tmp <- ScaleData(obj_tmp, verbose = F)
  obj_tmp <- RunPCA(obj_tmp, npcs = 20, verbose = F)
  obj_tmp <- FindNeighbors(obj_tmp, dims = 1:20, verbose = F)
  
  invisible(sweep.res.list_kidney <- paramSweep_v3(obj_tmp, PCs = 1:10, sct = FALSE))
  invisible(sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE))
  invisible(bcmvn_kidney <- find.pK(sweep.stats_kidney))
  
  annotations <- Idents(obj_tmp)
  invisible(homotypic.prop <- modelHomotypic(annotations))
  nExp_poi <- round(0.075*nrow(obj_tmp@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  invisible(obj_tmp <- doubletFinder_v3(obj_tmp, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE))
  pANN.col.name <- colnames(obj_tmp@meta.data)[4]
  invisible(obj_tmp <- doubletFinder_v3(obj_tmp, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pANN.col.name, sct = FALSE))
  
  obj_tmp@meta.data[,"DF_hi.lo"] <- obj_tmp@meta.data[,5]
  obj_tmp@meta.data$DF_hi.lo[which(obj_tmp@meta.data$DF_hi.lo == "Doublet" & obj_tmp@meta.data[,6] == "Singlet")] <- "Doublet_lo"
  obj_tmp@meta.data$DF_hi.lo[which(obj_tmp@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
  
  print(table(obj_tmp@meta.data$DF_hi.lo))
  
  colnames(obj_tmp@meta.data)[4:6] <- c("pANN", "DF.classifications.1", "DF.classifications.2")
  
  meta_tmp = obj_tmp@meta.data
  
  if (input.level == "Strong") meta_tmp = subset(meta_tmp, meta_tmp$DF_hi.lo == "Singlet")
  if (input.level == "Weak") meta_tmp = subset(meta_tmp, meta_tmp$DF_hi.lo != "Doublet_hi")
  
  mat_tmp = obj_tmp@assays$RNA@counts[,row.names(meta_tmp)]
  return(mat_tmp)
}

library(S4Vectors)
merge.sparse <- function(...) {
  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()
  
  for (M in list(...)) {
    cnold <- colnames(M)
    rnold <- rownames(M)
    
    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)
    
    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)
    ind <- unname(which(M != 0,arr.ind=T))
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,M@x)
  }
  
  sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
}

# 读取数据
mat1 = Read10X("/path/to/data/con/filtered_feature_bc_matrix")
mat2 = Read10X("/path/to/data/exp/filtered_feature_bc_matrix")

mat1_dblt = obj_filterdoublet(mat1, input.level = "Weak")
mat2_dblt = obj_filterdoublet(mat2, input.level = "Weak")

colnames(mat1_dblt) = paste0("Con_", colnames(mat1_dblt))
colnames(mat2_dblt) = paste0("Exp_", colnames(mat2_dblt))

All_mat <- merge.sparse(
  mat1_dblt,
  mat2_dblt
)

save(All_mat, file = "/path/to/output/All_mat_rmdoublet.rda")

# 重启并加载Seurat 5
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
load("/path/to/output/All_mat_rmdoublet.rda")

All_obj = CreateSeuratObject(All_mat, min.features = 300, min.cells = ncol(All_mat) * 0.0001)
dim(All_mat)
dim(All_obj)

All_obj$orig.ident = "Con"
All_obj$orig.ident[colnames(All_obj) %like% "Exp_"] = "Exp"
table(All_obj$orig.ident)

All_obj@meta.data[1:3,]
All_obj[["percent.mt"]] <- PercentageFeatureSet(All_obj, pattern = "^mt-")
Idents(All_obj) = All_obj$orig.ident
VlnPlot(All_obj, features = "percent.mt", ncol = 1, pt.size = 0) & 
  geom_boxplot(outlier.shape = NA) &   
  theme(legend.position="none")
quantile(All_obj$percent.mt, 0.99)

VlnPlot(All_obj, features = c("nCount_RNA","nFeature_RNA"), ncol = 2, pt.size = 0) & 
  geom_boxplot(outlier.shape = NA) 
quantile(All_obj$nCount_RNA, 0.99)
quantile(All_obj$nFeature_RNA, 0.99)
quantile(All_obj$nCount_RNA)
quantile(All_obj$nFeature_RNA)
sum(All_obj$nFeature_RNA < 9000)
sum(All_obj$nCount_RNA < 20000)

All_obj = subset(All_obj, nCount_RNA < 20000)
All_obj = subset(All_obj, nFeature_RNA < 10000)

All_obj <- NormalizeData(All_obj)
All_obj <- FindVariableFeatures(All_obj)
All_obj <- ScaleData(All_obj)
All_obj <- RunPCA(All_obj, verbose = FALSE)

library(harmony)
table(All_obj$orig.ident)
All_obj <- All_obj %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(All_obj, 'harmony')
harmony_embeddings[1:5, 1:5]

All_obj <- All_obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

All_obj <- All_obj %>% FindClusters(resolution = 0.1) 
All_obj <- All_obj %>% FindClusters(resolution = 1) 
All_obj <- All_obj %>% FindClusters(resolution = 1.5) 
All_obj <- All_obj %>% FindClusters(resolution = 2) 
All_obj <- All_obj %>% FindClusters(resolution = 5) 

DimPlot(All_obj, reduction = "umap", label = TRUE, pt.size = .0001, ncol = 3,
        group.by = c("celltype_l1","RNA_snn_res.1.5","RNA_snn_res.5"), raster=F) & theme(aspect.ratio=1)

DimPlot(All_obj, reduction = "umap", label = TRUE, pt.size = .0001, ncol = 1,
        group.by = c("RNA_snn_res.0.1"), raster=F) & theme(aspect.ratio=1)

library(dittoSeq)
dittoBarPlot(All_obj, "celltype_l1", group.by = "orig.ident", scale = "percent") +
  dittoBarPlot(All_obj, "celltype_l1", group.by = "orig.ident", scale = "count") 

library(scProgram)
Idents(All_obj) = All_obj$RNA_snn_res.1.5
FeatureMatrix = GetFeatures(obj = All_obj, group.by = "RNA_snn_res.1.5", genenumber = 100, pct_exp = 0.05, mode = "fast")

Idents(All_obj) = All_obj$RNA_snn_res.5
FeatureMatrix.2 = GetFeatures(obj = All_obj, group.by = "RNA_snn_res.5", genenumber = 100, pct_exp = 0.05, mode = "fast")

Idents(All_obj) = All_obj$RNA_snn_res.0.1
FeatureMatrix.0 = GetFeatures(obj = All_obj, group.by = "RNA_snn_res.0.1", genenumber = 100, pct_exp = 0.05, mode = "fast")

All_obj$celltype_l1 = as.character(All_obj$RNA_snn_res.1.5)
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(20,23)] = "T"
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(30,14,24,0,4,5,27,17,13,18,19,9,15,27)] = "Macro/Mono"
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(29)] = "Endo"
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(3,25,16)] = "Cancer"
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(26)] = "NK"
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(21,22,11,12,8,7,1,2,6,10,1)] = "NEU"
All_obj$celltype_l1[All_obj$celltype_l1 %in% c(28)] = "B"
DimPlot(All_obj, reduction = "umap", label = TRUE, pt.size = .1, group.by = c("celltype_l1"), ncol = 1) & theme(aspect.ratio=1)

load(file = "/path/to/output/All_obj.rda")

# SingleR注释
library(SingleR)
load("/path/to/reference/singler_MouseRNAseqData.rda")

ncell = 5000
cellid <- sample(1:ncol(All_obj), ncell, replace=F)
length(cellid)
obj_random <- All_obj[,cellid]
dim(obj_random)

anno.cell.main = SingleR(test=GetAssayData(obj_random, slot = "counts"), ref = ref.se, labels = ref.se$label.main)  
obj_random$pred_singler <- as.character(anno.cell.main$labels)
table(as.character(anno.cell.main$labels))
pred_singler.stat <- data.frame(table(obj_random$pred_singler))
pred_singler.stat <- subset(pred_singler.stat, pred_singler.stat$Freq > 20)
obj_random$pred_singler[!(obj_random$pred_singler %in% pred_singler.stat$Var1)] <- "other"

table(obj_random$pred_singler)
tmp <- data.frame(table(obj_random$pred_singler))
DimPlot(obj_random, reduction = "umap", label = TRUE, pt.size = .5, group.by = c("pred_singler"), ncol = 1)

anno.cell.sub = SingleR(test=GetAssayData(obj_random, slot = "counts"), ref = ref.se, labels = ref.se$label.fine)  
obj_random$pred_singler_sub <- as.character(anno.cell.sub$labels)
table(as.character(anno.cell.sub$labels))
pred_singler.stat <- data.frame(table(obj_random$pred_singler_sub))
pred_singler.stat <- subset(pred_singler.stat, pred_singler.stat$Freq > 5)
obj_random$pred_singler_sub[!(obj_random$pred_singler_sub %in% pred_singler.stat$Var1)] <- "other"

table(obj_random$pred_singler_sub)
tmp <- data.frame(table(obj_random$pred_singler_sub))
DimPlot(obj_random, reduction = "umap", label = TRUE, pt.size = .5, group.by = c("pred_singler_sub"), ncol = 1)

# CD45阳性细胞
table(All_obj$celltype_l1)
obj_cd45 = subset(All_obj, celltype_l1 %in% c("B","Macro/Mono","NEU","NK","T"))

library(dittoSeq)
dittoBarPlot(obj_cd45, "celltype_l1", group.by = "orig.ident", scale = "percent") +
  dittoBarPlot(obj_cd45, "celltype_l1", group.by = "orig.ident", scale = "count")
