# Creating results (first run Alignment script!)

## STEP 0 (Re-do for results)
{
  library(ORFik)
  library(RiboCrypt)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(org.Dr.eg.db)
  library(clusterProfiler)
  library(ggrepel)
  library(pheatmap)
  conf <- config.exper(experiment = "Hassan_Danio_rerio",
                       type = "RNA-seq",
                       assembly = "Danio_rerio_GRCz11")
}

## STEP 5 (QC)
# Run non batch corrected correlation and pca
df <- read.experiment(conf["exp RNA-seq"])
pair_data <- countTable(df, "mrna", type = "fpkm")
colnames(pair_data) <- gsub("RNA_|_|r|f", "", colnames(pair_data))
colnames(pair_data) <- gsub("^E", "Mt", colnames(pair_data))
colnames(pair_data) <- gsub("^W", "Wt", colnames(pair_data))


fig1_a <- ORFik:::correlation.plots(df, output.dir = QCfolder(df), complex.correlation.plots = FALSE,
                                    as_gg_list = TRUE, data_for_pairs = pair_data,
                                    label_size = 3.9)[[1]] + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.border = element_blank())
pcaExperiment(df, output.dir = QCfolder(df), table = countTable(df, "mrna", type = "fpkm"),
              title = "PCA analysis of mRNA fpkm")
# STEP 6 Prepare DESeq
RNA_counts <- countTable(df, type = "summarized")
colData(RNA_counts)[1:4, "replicate"] <- c(2,3,4,1)
dds_RNA <- DESeq(DESeqDataSet(RNA_counts, design = design(df, as.formula = TRUE)))
contrast <- c(design(df)[1], unique(df[,design(df)[1]]))
# dds_RNA_rep <- DESeq(DESeqDataSet(RNA_counts, design = design(df,
#                       as.formula = TRUE, batch.correction.design = TRUE)))
dds_RNA_rep <- DESeq(DESeqDataSet(RNA_counts, design = ~ condition + fraction + replicate))
## STEP 7 (DESEQ plots, PCA and MA)
# Batch correction plots
b <- vst(dds_RNA_rep)
colnames(b) <- gsub("RNA_|_|r|f", "", colnames(b)); colnames(b) <- gsub("^E", "Mt", colnames(b)); colnames(b) <- gsub("W", "Wt", colnames(b))
levels(colData(b)$condition) <- c("Mt", "Wt")
B <- b
colnames(colData(B))[4:5] <- c("Condition", "Sex")
colData(B)$Sex <- factor(rep(c("Female(F)", "Male(M)", "Female(F)", "Male(M)"),  each = 4))
heat_plot <- pheatmap(assay(B), kmeans_k = 500, cluster_rows=FALSE, show_rownames=FALSE, annotation_col=as.data.frame(colData(B)[, c("Condition", "Sex")]));heat_plot
B_pca <- B
plot <- plotPCA(B_pca, intgroup=c("Sex", "Condition")) + geom_point(aes(shape=B_pca$replicate), size = 4); plot
ggsave(filename = file.path(QCfolder(df), "PCAplot_raw_Hassan_Danio_rerio_RNA-seq.pdf"), plot)
assay(B_pca) <- limma::removeBatchEffect(assay(B_pca), batch = B_pca$replicate)
plotPCA(B_pca, intgroup=c("Condition", "Sex")) + geom_point(aes(shape=B_pca$replicate), size = 4)
assay(B_pca) <- limma::removeBatchEffect(assay(B_pca), B_pca$Sex)
plotPCA(B_pca, intgroup=c("Condition", "Sex")) + geom_point(aes(size = 4)) +
  ggtitle("PCA analysis of mRNA fpkm", paste("Numer of genes/regions:", nrow(B_pca)))
plot <- plotPCA(B_pca, intgroup=c("Condition")) + geom_point(size = 5) +
  guides(color=guide_legend(title="Condition")) + theme_classic()
ggsave(filename = file.path(QCfolder(df), "PCAplot_BatchCor_Hassan_Danio_rerio_RNA-seq.pdf"), plot)

fig1 <- cowplot::plot_grid(plotlist = list(fig1_a,
                                           cowplot::plot_grid(heat_plot$gtable, plot, ncol = 1, labels = c("B", "C"))),
                           ncol = 2, align = "h", rel_heights = c(4,1), rel_widths = c(2,1),
                           labels = c("A", NULL));fig1
ggsave(filename = file.path(QCfolder(df), "fig1_QC.pdf"), fig1, width = 12, heigh = 7)
ggsave(filename = file.path(QCfolder(df), "fig1_QC.png"), fig1, width = 12, heigh = 7)
# Without replicate in design
res_deseq <- results(dds_RNA, contrast = contrast)
res_deseq <- lfcShrink(dds_RNA, contrast=contrast,
                       res=res_deseq, type = "normal")
summary(res_deseq)
plotMA(res_deseq)
saveRDS(res_deseq, file.path(QCfolder(df), "DESeq_results.rds"))
plotMA(res_deseq)


DEG <- ORFik:::DEG.analysis(df, output.dir = NULL, batch.effect = TRUE)
DEG[, GeneID := txNamesToGeneNames(id, df)]
DEG <-DEG[!duplicated(GeneID),]
#DEG[, symbol := clusterProfiler::bitr(GeneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=org.Dr.eg.db)]
symbols <- clusterProfiler::bitr(DEG$GeneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=org.Dr.eg.db)
symbols <- symbols[!duplicated(symbols$ENSEMBL),]
DEG <- DEG[GeneID %in% symbols$ENSEMBL,]
DEG[, symbols := symbols$SYMBOL[chmatch(symbols$ENSEMBL, GeneID)]]
give_alpha <- (DEG$padj < 0.01) + 1
volcano_plot <- ggplot(DEG, aes(x = LFC, y = -log10(padj), label = symbols)) +
  ylim(c(0,30)) + xlim(c(-2.8, 2.8)) +
  theme_bw() +
  geom_point(alpha = c(0.1, 0.7)[give_alpha], color = c("gray", "green4")[give_alpha]) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(-1.7, 1.7), linetype="dotted", alpha = 0.9, col = "red") +
  geom_text_repel(data = subset(DEG, abs(LFC) > 1.7 & padj < 0.000001),
                  box.padding   = 0.35,
                  point.padding = 0.5,
                  segment.color = 'grey50'); volcano_plot
ggsave(filename = file.path(QCfolder(df), "fig3_volcano.pdf"), volcano_plot, width = 8, heigh = 6)
fig1_alt <- cowplot::plot_grid(plotlist = list(fig1, cowplot::plot_grid(plotlist = list(volcano_plot), labels = "D")),
                           ncol = 1, align = "v");fig1_alt
ggsave(filename = file.path(QCfolder(df), "fig1_QC_volc.pdf"), fig1_alt, width = 12, heigh = 12)
ggsave(filename = file.path(QCfolder(df), "fig1_QC_volc.png"), fig1_alt, width = 12, heigh = 12)

fig1_alt <- cowplot::plot_grid(plotlist = list(fig1_a, heat_plot$gtable, plot, volcano_plot),
                               ncol = 2, align = "v", );fig1_alt
# gridExtra::grid.arrange(arrangeGrob(fig1_a, left = textGrob("A", x = unit(1, "npc"),
#                                                        y = unit(.95, "npc"))),
#                    arrangeGrob(heat_plot$gtable, left = textGrob("B", x = unit(1, "npc"),
#                                                        y = unit(.95, "npc"))),
#                    arrangeGrob(plot, left = textGrob("C", x = unit(1, "npc"),
#                                                        y = unit(.95, "npc"))),
#                    arrangeGrob(volcano_plot, left = textGrob("D", x = unit(1, "npc"),
#                                                        y = unit(.95, "npc"))),
#                    layout_matrix = rbind(c(1,1,NA,2,2),
#                                          c(1,1,3,3,3),
#                                          c(4,4,4,4,4),
#                                          c(4,4,4,4,4)))

# With replicate in design
res_deseq_rep <- results(dds_RNA_rep, contrast = contrast)
res_deseq_rep <- lfcShrink(dds_RNA_rep, contrast=contrast,
                           res=res_deseq_rep, type = "normal")
summary(res_deseq_rep)
plotMA(res_deseq_rep)
saveRDS(res_deseq_rep, file.path(QCfolder(df), "DESeq_results_rep.rds"))
res_deseq_rep <- ORFik::DEG.analysis(df)
plot <- RiboCrypt::DEG_plot(res_deseq_rep)
htmlwidgets::saveWidget(plot,
                        file.path(QCfolder(df), "DESeq_interactive_sig_only.html"))
## Load res
res_deseq <- readRDS(file.path(QCfolder(df), "DESeq_results.rds"))
res_deseq_rep <- readRDS(file.path(QCfolder(df), "DESeq_results_rep.rds"))
## Check GLA and results by gene
res_deseq_rep_gene <- res_deseq_rep
rownames(res_deseq_rep_gene) <- txNamesToGeneNames(rownames(res_deseq_rep), df)
res_deseq_rep_gene <- res_deseq_rep_gene[!duplicated(rownames(res_deseq_rep_gene)),]
res_deseq_rep_gene["ENSDARG00000036155",]
res_deseq_rep_gene_high <- res_deseq_rep_gene[which(abs(res_deseq_rep_gene$log2FoldChange) > 1),]
# Save final expression data as csv
saveRDS(res_deseq_rep_gene, file.path(QCfolder(df), "DESeq_results_rep_gene_all.rds"))
fwrite(as.data.table(res_deseq_rep_gene, keep.rownames = TRUE),
       file.path(QCfolder(df), "DESeq_results_rep_gene_all.csv"))
saveRDS(res_deseq_rep_gene_high, file.path(QCfolder(df), "DESeq_results_rep_gene_high.rds"))
fwrite(as.data.table(res_deseq_rep_gene_high, keep.rownames = TRUE),
       file.path(QCfolder(df), "DESeq_results_rep_gene_high.csv"))
# GLA levels LFC ~ -1.3, a 60% reduction in transcript levels.
summary(res_deseq_rep_gene)
summary(res_deseq_rep_gene_high)
res_deseq_rep_gene_high["ENSDARG00000036155",]

## STEP 10 (GO analysis)
## For GOrilla DB online:
gorilla <- res_deseq_rep$log2FoldChange
names(gorilla) <- txNamesToGeneNames(rownames(res_deseq_rep), df)
gorilla <- gorilla[which(res_deseq_rep$padj < 0.1)]
gorilla <- sort(gorilla, decreasing = TRUE)
gorilla <- gorilla[!duplicated(names(gorilla))]
symbols <- clusterProfiler::bitr(names(gorilla), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=org.Dr.eg.db)
gorilla <- gorilla[as.integer(rownames(symbols))]
names(gorilla) <- symbols$SYMBOL
write.table(names(gorilla), file = file.path(QCfolder(df), "Gorilla.txt"),
            col.names = FALSE, sep = "", row.names = FALSE, quote = FALSE)
write.table(names(gorilla[abs(gorilla) > 1]), file = file.path(QCfolder(df), "Gorilla_high.txt"),
            col.names = FALSE, sep = "", row.names = FALSE, quote = FALSE)
## Ranked list, all:
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/49fpsq4h/GOResults.html
## Unranked lists: Target LFC > 1 vs background (all):
#http://cbl-gorilla.cs.technion.ac.il/GOrilla/4i2ucshv/GOResults.html

## Cluster profiler
gene_list <- res_deseq_rep$log2FoldChange
names(gene_list) <- txNamesToGeneNames(rownames(res_deseq_rep), df)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]
# All
gse <- gseGO(geneList=gene_list,
             ont ="ALL",
             keyType = "ENSEMBL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Dr.eg.db,
             pAdjustMethod = "BH")
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
# Compartment only (CC)
gse <- gseGO(geneList=gene_list,
             ont ="CC",
             keyType = "ENSEMBL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Dr.eg.db,
             pAdjustMethod = "BH")
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

## STEP 11 (Subcellular compartment)
# heatmap cluster of differentiated genes (control vs mutant) indicated by the affected subcellular compartment:
# lysosome, mitochondria, tight junctions, cytoskeleton, oxidative stress, peroxisome, endoplasmic reticulum, Golgi app. Cell cycle
# These next steps are done ONLINE in GO BROWSER, described in methods section
