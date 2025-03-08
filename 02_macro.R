obj_macro = subset(All_obj, celltype_l1 == "Macro/Mono")

obj_macro <- obj_macro %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

obj_macro <- obj_macro %>% FindClusters(resolution = 0.1) 
obj_macro <- obj_macro %>% FindClusters(resolution = 0.2) 

DimPlot(obj_macro, reduction = "umap", label = TRUE, pt.size = .0001, ncol = 2,
        group.by = c("RNA_snn_res.0.1","RNA_snn_res.0.2"), raster=F) & theme(aspect.ratio=1)

library(scProgram)
Idents(obj_macro) = obj_macro$RNA_snn_res.0.2
FeatureMatrix.3 = GetFeatures(obj = obj_macro, group.by = "RNA_snn_res.0.2", genenumber = 100, pct_exp = 0.05, mode = "fast")

obj_macro$celltype_l2 = as.character(obj_macro$RNA_snn_res.0.2)
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(0)] = "Cx3cr1+Macrophages"
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(1)] = "Mrc1+Macrophages"
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(2)] = "Monocytes"
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(3)] = "DC"
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(4)] = "Mki67+Macrophages"
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(5)] = "Ifit1+Macrophages"
obj_macro$celltype_l2[obj_macro$RNA_snn_res.0.2 %in% c(6)] = "MHCII+Macrophages"

DimPlot(obj_macro, reduction = "umap", label = TRUE, pt.size = .1, group.by = c("celltype_l2"), ncol = 1) & theme(aspect.ratio=1)

load(file = "/path/to/data/obj_macro.rda")

library(dittoSeq)
dittoBarPlot(obj_macro, "celltype_l2", group.by = "orig.ident", scale = "percent") +
  dittoBarPlot(obj_macro, "celltype_l2", group.by = "orig.ident", scale = "count") 

library(UCell)

signature <- qusage::read.gmt("/path/to/signature/APC.gmt")
signature2 <- qusage::read.gmt("/path/to/signature/M1M2_signature.gmt")

signature = c(signature, signature2)

obj_DC = subset(obj_macro, celltype_l2 == "DC")

exprmat = GetAssayData(obj_DC, slot = "counts")
row.names(exprmat) = toupper(row.names(exprmat))
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(exprmat, features = signature, maxRank = 3000, name = "")))

obj_macro2 = obj_DC
obj_macro2@meta.data = cbind(obj_macro2@meta.data, signature.matrix.neu)

Idents(obj_macro2) = obj_macro2$orig.ident
VlnPlot(obj_macro2, features = colnames(signature.matrix.neu), ncol = 5, pt.size = 0)& 
  theme(aspect.ratio=1) &
  geom_boxplot(outlier.shape = NA)

wilcox.test(obj_macro2$Antigen.processing..cross.presentation ~ obj_macro2$orig.ident)

obj_macrophage = subset(obj_macro, celltype_l2 != "DC")

exprmat = GetAssayData(obj_macrophage, slot = "counts")
row.names(exprmat) = toupper(row.names(exprmat))
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(exprmat, features = signature, maxRank = 3000, name = "")))

obj_macro2 = obj_macrophage
obj_macro2@meta.data = cbind(obj_macro2@meta.data, signature.matrix.neu)

Idents(obj_macro2) = obj_macro2$orig.ident
VlnPlot(obj_macro2, features = colnames(signature.matrix.neu), ncol = 5, pt.size = 0)& 
  theme(aspect.ratio=1) &
  geom_boxplot(outlier.shape = NA)

t.test(obj_macro2$Antigen.processing..cross.presentation ~ obj_macro2$orig.ident)

library(SingleR)
load("/path/to/reference/singler_MouseRNAseqData.rda")

ncell = 5000
cellid <- sample(1:ncol(obj_macro), ncell, replace=F)
length(cellid)
obj_random <- obj_macro[,cellid]
dim(obj_random)

anno.cell.main=SingleR(test=GetAssayData(obj_random, slot = "counts"), ref = ref.se, labels = ref.se$label.main)  
obj_random$pred_singler <- as.character(anno.cell.main$labels)
table(as.character(anno.cell.main$labels))
pred_singler.stat <- data.frame(table(obj_random$pred_singler))
pred_singler.stat <- subset(pred_singler.stat, pred_singler.stat$Freq > 20)
obj_random$pred_singler[!(obj_random$pred_singler %in% pred_singler.stat$Var1)] <- "other"

table(obj_random$pred_singler)
tmp <- data.frame(table(obj_random$pred_singler))
DimPlot(obj_random, reduction = "umap", label = TRUE, pt.size = .5, group.by = c("pred_singler"), ncol = 1)
