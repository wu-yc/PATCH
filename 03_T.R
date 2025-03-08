obj_T = subset(All_obj, celltype_l1 == "T")

obj_T <- obj_T %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

obj_T <- obj_T %>% FindClusters(resolution = 0.1) 
obj_T <- obj_T %>% FindClusters(resolution = 0.2) 
obj_T <- obj_T %>% FindClusters(resolution = 0.6) 

DimPlot(obj_T, reduction = "umap", label = TRUE, pt.size = .1, ncol = 2,
        group.by = c("RNA_snn_res.0.6"), raster=F) & theme(aspect.ratio=1)

library(scProgram)
Idents(obj_T) = obj_T$RNA_snn_res.0.6
FeatureMatrix.3 = GetFeatures(obj = obj_T, group.by = "RNA_snn_res.0.6", genenumber = 100, pct_exp = 0.05, mode = "fast")

library(SingleR)
load("/path/to/reference/singler_data.rda")

obj_random<-obj_T

anno.cell.main=SingleR(test=GetAssayData(obj_random, slot = "counts") , ref = ref.se, labels = ref.se$label.main)  
obj_random$pred_singler<-as.character(anno.cell.main$labels)
table(as.character(anno.cell.main$labels))
pred_singler.stat<-data.frame(table(obj_random$pred_singler))
pred_singler.stat<-subset(pred_singler.stat, pred_singler.stat$Freq > 20)
obj_random$pred_singler[!(obj_random$pred_singler %in% pred_singler.stat$Var1)] <- "other"

table(obj_random$pred_singler)
tmp<-data.frame(table(obj_random$pred_singler))
DimPlot(obj_random, reduction = "umap", label = TRUE, pt.size = .5, group.by = c("pred_singler"), ncol = 1)

anno.cell.sub=SingleR(test=GetAssayData(obj_random, slot = "counts") , ref = ref.se, labels = ref.se$label.fine)  
obj_random$pred_singler_sub<-as.character(anno.cell.sub$labels)
table(as.character(anno.cell.sub$labels))
pred_singler.stat<-data.frame(table(obj_random$pred_singler_sub))
pred_singler.stat<-subset(pred_singler.stat, pred_singler.stat$Freq > 5)
obj_random$pred_singler_sub[!(obj_random$pred_singler_sub %in% pred_singler.stat$Var1)] <- "other"

table(obj_random$pred_singler_sub)
tmp<-data.frame(table(obj_random$pred_singler_sub))
DimPlot(obj_random, reduction = "umap", label = TRUE, pt.size = .5, group.by = c("pred_singler_sub"), ncol = 1)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')

obj_T$celltype_l2 = as.character(obj_T$RNA_snn_res.0.6)
obj_T$celltype_l2[obj_T$RNA_snn_res.0.6 %in% c(0,5)] = "Tcf7+Naive"
obj_T$celltype_l2[obj_T$RNA_snn_res.0.6 %in% c(1,3)] = "Foxp3+Treg"
obj_T$celltype_l2[obj_T$RNA_snn_res.0.6 %in% c(2)] = "Gzmb+Effector"
obj_T$celltype_l2[obj_T$RNA_snn_res.0.6 %in% c(4)] = "Pdcd1+Exhausted"
obj_T$celltype_l2[obj_T$RNA_snn_res.0.6 %in% c(6)] = "Mki67+Proliferative"

DimPlot(obj_T, reduction = "umap", label = F, pt.size = .6, group.by = c("celltype_l2"), ncol = 1, cols = my36colors) & theme(aspect.ratio=1)

load(file = "/path/to/data/obj_T.rda")
library(dittoSeq)
dittoBarPlot(obj_T, "celltype_l2", group.by = "orig.ident", scale = "percent") +scale_fill_manual(values=my36colors)
dittoBarPlot(obj_T, "celltype_l2", group.by = "orig.ident", scale = "count") +scale_fill_manual(values=my36colors) +NULL

signature <- qusage::read.gmt("/path/to/signature/T_basic_signature.gmt")
input.genes = c("Ifng","Gnly","Gzmb","Gzma","Cxcl13","Tigit","Tnfrsf9","Pdcd1","Ctla4","Lag3","Havcr2")
input.genes = c("Tigit","Tnfrsf9","Pdcd1","Ctla4","Lag3","Havcr2")
input.genes = c("Pdcd1","Lag3","Havcr2","Batf", "Ccr7","Tcf7","Lef1","Sell")
input.genes = unique(c("Tigit","Tnfrsf9","Pdcd1","Ctla4","Lag3","Havcr2","Pdcd1","Lag3","Havcr2","Icos","Batf"))
input.genes = c("Pdcd1","Lag3","Havcr2","Batf")

DotPlot(obj_T, features = input.genes, cols = c("lightgrey", "darkgreen"), group.by = "orig.ident") &RotatedAxis()

input.genes = c("Ccr7","Tcf7","Lef1","Sell","Il7r")
DotPlot(obj_T, features = input.genes, cols = c("lightgrey", "darkgreen"), group.by = "orig.ident") &RotatedAxis()

input.genes = signature$Exhaustion
DotPlot(obj_T, features = input.genes, cols = c("lightgrey", "darkgreen"), group.by = "orig.ident")

Idents(obj_T) = obj_T$orig.ident
VlnPlot(obj_T, features = input.genes, ncol = 5, pt.size = 0)& 
  theme(aspect.ratio=1) &
  geom_boxplot(outlier.shape = NA)

Idents(obj_T) = obj_T$orig.ident
diff.metab.celltype <- FindMarkers(obj_T, ident.1 ="Group1", ident.2 ="Group2", logfc.threshold = -1, return.thresh = 1)
diff.metab.celltype$gene = row.names(diff.metab.celltype)

library(UCell)
signature <- qusage::read.gmt("/path/to/signature/T_basic_signature.gmt")

exprmat = GetAssayData(obj_T, slot = "counts")
row.names(exprmat) = toupper(row.names(exprmat))
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(exprmat, features = signature, maxRank = 3000, name = "")))

obj_T2 = obj_T
obj_T2@meta.data = cbind(obj_T2@meta.data, signature.matrix.neu)

Idents(obj_T2) = obj_T2$orig.ident
VlnPlot(obj_T2, features = colnames(signature.matrix.neu), ncol = 2, pt.size = 0)& 
  theme(aspect.ratio=1) &
  geom_boxplot(outlier.shape = NA)
