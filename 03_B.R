library(UCell)

signature <- qusage::read.gmt("/path/to/signature/APC.gmt")
signature2 <- qusage::read.gmt("/path/to/signature/M1M2_signature.gmt")

signature = c(signature, signature2)

exprmat = GetAssayData(obj_B, slot = "counts")
row.names(exprmat) = toupper(row.names(exprmat))
signature.matrix.neu <- (data.frame(ScoreSignatures_UCell(exprmat, features = signature, maxRank = 3000, name = "")))

obj_macro2 = obj_B
obj_macro2@meta.data = cbind(obj_macro2@meta.data, signature.matrix.neu)

Idents(obj_macro2) = obj_macro2$orig.ident
VlnPlot(obj_macro2, features = colnames(signature.matrix.neu), ncol = 5, pt.size = 0)& 
  theme(aspect.ratio=1) &
  geom_boxplot(outlier.shape = NA)

t.test(obj_macro2$Antigen.processing..cross.presentation ~ obj_macro2$orig.ident)
