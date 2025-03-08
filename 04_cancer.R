table(All_obj$celltype_l1)

obj_cancer = subset(All_obj, celltype_l1 == "Cancer")

library(UCell)
signature <- qusage::read.gmt("/path/to/signature/file.gmt")
signature = signature[names(signature)[names(signature) %like% "HALLMARK"]]

exprmat = GetAssayData(obj_cancer, slot = "counts")
row.names(exprmat) = toupper(row.names(exprmat))
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(exprmat, features = signature, maxRank = 3000, name = "")))

obj_cancer2 = CreateSeuratObject(t(signature.matrix.neu), meta.data = obj_cancer@meta.data)
obj_cancer2 = NormalizeData(obj_cancer2)

Idents(obj_cancer2) = obj_cancer2$orig.ident
diff.metab.celltype <- FindMarkers(obj_cancer2, ident.1 ="Group1", ident.2 ="Group2", logfc.threshold = -1, return.thresh = 1)
diff.metab.celltype$gene = row.names(diff.metab.celltype)

diff.metab.celltype$logp = -log(diff.metab.celltype$p_val)

diff.metab.celltype$sig = "ns"
diff.metab.celltype$sig[diff.metab.celltype$avg_log2FC > 0.1 & diff.metab.celltype$p_val < .05] = "up"
diff.metab.celltype$sig[diff.metab.celltype$avg_log2FC < -0.1 & diff.metab.celltype$p_val < .05] = "down"
table(diff.metab.celltype$sig)

diff.metab.celltype$size = abs(diff.metab.celltype$avg_log2FC)

diff.metab.celltype_label = subset(diff.metab.celltype, sig != "ns")

quantile(diff.metab.celltype$avg_log2FC)

ggplot(data = diff.metab.celltype, aes(x = avg_log2FC, y = logp, color = sig))+
  geom_point(size = 2)+
  ggrepel::geom_text_repel(data = diff.metab.celltype_label, mapping = aes(label = gene), show.legend = F, max.overlaps = 100000, size = 3) + 
  theme_bw()+
  theme(aspect.ratio=1)+
  scale_color_manual(values = c(my36colors[3], "grey50", my36colors[2])) +
  theme(legend.position="right", legend.box = "right", axis.text.x = element_text(angle = 45, hjust = 1))+
  xlim(-.31, .31)+
  NULL

Idents(obj_cancer2) = obj_cancer2$orig.ident
VlnPlot(obj_cancer2, features = diff.metab.celltype_label$gene, ncol = 4, pt.size = 0)& 
  theme(aspect.ratio=1) &
  geom_boxplot(outlier.shape = NA)