#### Rccd1 OE -

x <- read.csv("~/Documents/BRCA1_OE/BR_RCCD1OE_1000DEG.csv")
rownames(x) <- x$unique_id
x <- x[, -c(1:2)]
y <- read.csv("~/Documents/BRCA1_OE/BR_RCCD1OE_logFC.csv")
y <- y[y$X1 %in% rownames(x), ]
x <- x[rownames(x) %in% y$X1,]
x$FDR <- y$padj_RCCD1OE.CONTROL[match(rownames(x), y$X1)]
x$FC <- y$logFC_RCCD1OE.CONTROL[match(rownames(x), y$X1)]

aux <- x[order(x$FDR), ]
aux$Genes <- rownames(aux)
aux <- aux[ ,c(7,1:6)]
write.table(aux, sep = "\t", quote = F, row.names = F, file = "~/brca1_oe/TableDEG_BR_top1000.txt")


## top 100 genes for pathway
write.table(rownames(x[order(x$FDR), ])[1:100], quote = F, row.names = F, col.names = F,
            file = "~/brca1_oe/top100_OE_BR.txt")

aux <- x[order(x$FDR)[1:100], ]
aux$AveExp <- apply(aux[, 1:4], 1, mean)
tmp <- data.frame(ID = rownames(aux), logFC = aux$FC, AveExpr = aux$AveExp, adj.P.Val = aux$FDR)

pdf("~/RCCD1_paper/HM_RCCD1OE_BR_TOP100.pdf", 15,17)
pheatmap(t(scale(t(aux[1:100, 1:4]))), color = viridis(100), main = 
           "Top 50 genes - RCCD1 Over-expression vs control - RNASeq", angle_col = 45, border_color = NA)
dev.off()         

xx <- read.csv("~/Documents/BRCA1_OE/FT_RCCD1OE_DEG1000.csv")
rownames(xx) <- xx$unique_id 
xx <- xx[, -c(1:2)]
yy <- read.csv("~/Documents/BRCA1_OE/FT-RCCD1OE_logFC.csv")
yy <- yy[yy$X1 %in% rownames(xx), ]
xx <- xx[rownames(xx) %in% yy$X1, ]
xx$FDR <- yy$padj_RCCD1OE.CONTROL[match(rownames(xx), yy$X1)]
xx$FC <- yy$logFC_RCCD1OE.CONTROL[match(rownames(xx), yy$X1)]

aux <- xx[order(xx$FDR)[1:100], ]
pdf("~/RCCD1_paper/HM_RCCD1OE_FT_TOP100.pdf", 15,17)
pheatmap(t(scale(t(aux[1:100, 1:4]))), color = viridis(100), main = 
           "Top 50 genes - FT - RCCD1 Over-expression vs control - RNASeq", angle_col = 45, border_color = NA)
dev.off()        

# aux <- xx[order(xx$FDR), ]
# aux$Genes <- rownames(aux)
# aux <- aux[ ,c(7,1:6)]
# write.table(aux, sep = "\t", quote = F, row.names = F, file = "~/brca1_oe/TableDEG_FTSEC_top1000.txt")



######## HEATMAP OF SHARED DEG GENES ###########
oe.counts <- fread("~/Documents/BRCA1_OE/RCCD1_OE_COUNTS.csv")
oe.counts$V1 <- substr(oe.counts$V1, 20, nchar(oe.counts$V1))
h.shared.deg <- oe.counts[oe.counts$V1 %in% shared.genes, ]
pdf("~/RCCD1_paper/HM_RCCD1OE_Shared_DEG.pdf", 15,17)
pheatmap(t(scale(t(h.shared.deg[, 2:dim(h.shared.deg)[2]]))), color = viridis(100), angle_col = 45,
         border_color = NA, main = "Shared DEG - RCCD1 Over-expression vs Control")
dev.off()






### pathway analysis

## first check if any of the p53 genes in the pathway are DEGs
p53 <- fread('~/RCCD1_paper/p53_GO_pathway.txt')
p53.path <- toupper(p53$`MGI Gene/Marker ID`)

table(rownames(x) %in% p53.path)
x[rownames(x) %in% p53.path,]












BiocManager::install(c("pathview", "gage", "gageData"))
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)

kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

## get entrez id from gene symbol to run pathway
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = rownames(x), mart = mart)

br_oe <- x[rownames(x) %in% shared.genes, ]
br_oe$entrez <- genes.mart$entrezgene_id[match(rownames(br_oe), genes.mart$hgnc_symbol)]
##

## prepare data for pathway 
fc.vector <- br_oe$FC
names(fc.vector) <- br_oe$entrez
fc.vector <- na.omit(fc.vector)



keggres = gage(fc.vector, gsets=kegg.sets.hs, same.dir=TRUE)


pathways = data.frame(id=rownames(keggres$greater), keggres$less)
head(pathways)

keggrespathways <- rownames(keggres$greater)[1:10]

# Extract the IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

pathview(gene.data=fc.vector, pathway.id=keggresids, species="hsa")


## go enrichment
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$MF]

gobpres = gage(fc.vector, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)
#

library(GOplot)
david <- read.delim("~/Documents/BRCA1_OE/david_shared_DEGs_OE.txt")

aux <- x[rownames(x) %in% b$b,]
aux <- aux[, c(5:6)]
aux$ID <- rownames(aux)



circ <- circle_dat(david, aux)

##
library(DOSE)
# data(geneList)
# de <- names(geneList)[abs(geneList) > 2]
genes.mart[genes.mart$hgnc_symbol %in% b$b, ]$entrezgene_id





edo <- enrichDGN(genes.mart[genes.mart$hgnc_symbol %in% b$b, ]$entrezgene_id)
library(enrichplot)
barplot(edo, showCategory=20)
dotplot(edo, showCategory=30)

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p3, ncol=2, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

# edo2 <- gseNCG(sort(fc.vector, decreasing = T), nPerm=100, pvalueCutoff = 1)
# p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
# p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
# plot_grid(p1, p2, ncol=2)



library(clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
p1 <- emapplot(xx)
p2 <- emapplot(xx,legend_n=2) 
p3 <- emapplot(xx,pie="count")
p4 <- emapplot(xx,pie="count", pie_scale=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


###############################################################
### Crispr data ###

rccd1.crispr <- read.csv("~/4C_RCCD1/RCCD1_COUNTS_edited_072120.csv", skip = 1)
rccd1.crispr2 <- read.csv("~/4C_RCCD1/BD-8579–01–06–2020_COUNTS_072120.csv", skip = 1)
aux <- read.csv("~/4C_RCCD1/BRCA_COUNTS_EDIT.csv", skip= 1)
aux <- aux[, c(1,9:12)]


rccd1.crispr <- rccd1.crispr[, c(1:5,18:23)]
rccd1.crispr <- merge(rccd1.crispr, rccd1.crispr2[, c(1,2)], by = "NEW.NAME") ### merge to get one replicate missing from the file -- jasmine/steph sent email 7/22/20
rccd1.crispr <- merge(rccd1.crispr, aux, by.x = "NEW.NAME", by.y = "X")



pca.crispr <- prcomp(t(rccd1.crispr[, 2:dim(rccd1.crispr)[2]]))
aux <- as.data.frame(pca.crispr$x)

ggplot(aux, aes(PC1, PC2, label = rownames(aux))) + geom_point() + geom_label_repel() + labs(title = "PCA - RCCD1 KO RNASeq") +
  xlab(paste0("PC1 ", prettyNum(summary(pca.crispr)$importance[2,1]*100, digits =2), "%")) +
  ylab(paste0("PC2 ", prettyNum(summary(pca)$importance[2,2]*100, digits =2), "%"))



library(DESeq2)
coldata.ko <- data.frame(row.names = colnames(rccd1.crispr)[2:16], condition = c(rep("FT282_KO",4), rep("MCF12A_KO",2), rep("MCF12A_OR",4), "MCF12A_KO", rep("FT282_OR",4)))

###                   FT282_KO vs FT282 OR                             ###
coldata.ko.sub <- subset(coldata.ko, coldata.ko$condition %in% "FT282_KO" | coldata.ko$condition %in% "FT282_OR")
tmp <- rccd1.crispr[, rownames(coldata.ko.sub)]
rownames(tmp) <- rccd1.crispr$NEW.NAME
dds <- DESeqDataSetFromMatrix(countData = round(tmp),
                              colData = coldata.ko.sub,
                              design =  ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)
res$STATUS <- "Not Sig"
res[res$log2FoldChange < -1 & res$padj < 0.05, ]$STATUS <- "Down"
res[res$log2FoldChange > 1 & res$padj < 0.05, ]$STATUS <- "Up"

ggplot(as.data.frame(res), aes(log2FoldChange, -log10(padj), col = STATUS)) + geom_point() + 
  scale_colour_manual(values = c("darkgreen", "grey", "darkred")) + labs(title = "DEGs - FT282_KO vs FT282_OR")
vst <- vst(dds)
hm <- assay(vst)
hm <- assay(dds)
hm <- hm[rownames(res)[!res$STATUS %in% "Not Sig"] , ]

pheatmap(t(scale(t(hm))), color = viridis(50), show_rownames = F, main= "FT282_KO vs FT282_OR")
genes.ko.ft282 <- rownames(hm)


###                                                             ###

###                   MCF12A_KO vs MCF12A OR                             ###
coldata.ko.sub1 <- subset(coldata.ko, coldata.ko$condition %in% "MCF12A_KO" | coldata.ko$condition %in% "MCF12A_OR")
tmp <- rccd1.crispr[, rownames(coldata.ko.sub1)]
rownames(tmp) <- rccd1.crispr$NEW.NAME
dds <- DESeqDataSetFromMatrix(countData = round(tmp),
                              colData = coldata.ko.sub1,
                              design =  ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)
res$STATUS <- "Not Sig"
res[res$log2FoldChange < -1 & res$padj < 0.05, ]$STATUS <- "Down"
res[res$log2FoldChange >1 & res$padj < 0.05, ]$STATUS <- "Up"

ggplot(as.data.frame(res), aes(log2FoldChange, -log10(padj), col = STATUS)) + geom_point() + 
  scale_colour_manual(values = c("darkgreen", "grey", "darkred")) + labs(title = "DEGs - MCF12A_KO vs MCF12A_OR")
vst <- vst(dds)
#hm <- assay(vst)
hm <- assay(dds)
hm <- hm[rownames(res)[!res$STATUS %in% "Not Sig"] , ]

pheatmap(t(scale(t(hm))), color = viridis(50), show_rownames = F, main= "MCF12A_KO vs MCF12A_OR")
genes.ko.mcf12a <- rownames(hm)


###                                                             ###

aux <- GenomicRanges::intersect(genes.ko.ft282, genes.ko.mcf12a)
aux <- gsub(".*_", "",aux)

tmp <- read.table("~/test/home/feliped/brca1_oe/TableDEG_BR_top1000.txt", stringsAsFactors = F, header = T)
tmp <- tmp[, tmp$FDR < 0.05 && (tmp$FC < 0 | tmp$FC > 0),]

tmp2 <- read.table("~/test/home/feliped/brca1_oe/TableDEG_FTSEC_top1000.txt", stringsAsFactors = F, header = T)
tmp2 <- tmp2[, tmp2$FDR < 0.05 && (tmp2$FC < 0 | tmp2$FC > 0),]

GenomicRanges::intersect(tmp$Genes, gsub(".*_", "", genes.ko.mcf12a))
GenomicRanges::intersect(tmp2$Genes,gsub(".*_", "", genes.ko.ft282))





## table for the paper with the shared deg
tt <- read_excel("~/RCCD1_paper/Shared breast ovarian cancer genes.xlsx")
ft.deg.oe <- read.table( "~/RCCD1_paper/Table_DEG_FT_RCCD1OE.txt" , stringsAsFactors = F)
br.deg.oe <- read.table("~/RCCD1_paper/Table_DEG_BR_RCCD1OE.txt", stringsAsFactors = F)

aux <- br.deg.oe[rownames(br.deg.oe) %in% tt$`GENE ID`, ]
aux <- aux[match(tt$`GENE ID`, rownames(aux)),]
tt$`Direction of differential expression` <- 
  paste0("logFC = ", round(aux$logFC_RCCD1OE.CONTROL, digits = 3), " Adj Pvalue = ", aux$padj_RCCD1OE.CONTROL)

aux <- ft.deg.oe[rownames(ft.deg.oe) %in% tt$`GENE ID`, ]
aux <- aux[match(tt$`GENE ID`, rownames(aux)),]
tt$...4 <- 
  paste0("logFC = ", round(aux$logFC_RCCD1OE.CONTROL, digits = 3), " Adj Pvalue = ", aux$padj_RCCD1OE.CONTROL)

write.table(tt, row.names = F, sep = "\t", quote = F, file = "~/RCCD1_paper/Shared_breast_ovarian_cancer_genes_FC_Pval.csv")






