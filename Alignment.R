## STEPS
# RNA-seq data analysis steps:
# Organism: Danio rerio, samples 16
# Author: Hassan
# 0. Set directories
# 1A. Download annotation and create STAR genome index
# 1B. Download read data
# 2. Align with STAR (paired end reads)
# 3. Create ORFik experiment, link Danio rerio annotation, design is (Condition and Sex, with 4x replicates)
# 4. Create count tables for DESeq (mRNA coverage)
# 5. Create QC plots (correlation and PCA (batch corrected and non batch corrected)
# 6. Run DESeq (default test and p-value cut off)
# 7. Plot DESeq results (MA)
# 8. Volcano plot
# 9. Run GSEA analysis (clusterProfiler)
# 10. Run GO analysis (GOrilla, online tool)

## STEP 0
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
## STEP 1A
annotation <- getGenomeAndAnnotation(organism = "Danio rerio",
                                     output.dir = conf["ref"])
index <- STAR.index(annotation)
## STEP 1B
## When data is published to read archive do this to get it:
#study <- download.SRA.metadata("PRJEB55250")
#download.SRA(study, conf["fastq RNA-seq"], rename = FALSE)

## STEP 2 (STAR alignment)
STAR.align.folder(conf["fastq RNA-seq"],
                  conf["bam RNA-seq"],
                  index, paired.end = TRUE,
                  steps = "ge",
                  include.subfolders = TRUE)
## STEP 3 (ORFik experiment)
create.experiment(dir = paste0(conf["bam RNA-seq"], "aligned"),
                  exper = conf["exp RNA-seq"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  fa = annotation["genome"],
                  organism = "Danio rerio", pairedEndBam = TRUE,
                  libtype = "RNA",
                  condition = rep(c("E", "W"), each = 8),
                  rep = rep(seq(4), 4),
                  fraction = rep(rep(c("F", "M"), each = 4), 2),
                  author = "Hassan")
df <- read.experiment(conf["exp RNA-seq"])
# STEP 4 (count tables)
dir.create(QCfolder(df), showWarnings = FALSE)
filepaths <- filepath(df, "default")
mrna <- loadRegion(df, "mrna")
Sys.time()
res <- lapply(filepaths, function(x) {
  print(x)
  return(countOverlaps(mrna, fimport(x)))
})
Sys.time()
data.table::setDT(res)
mat <- as.matrix(res); colnames(mat) <- NULL
# Make summarized experiment
varNames <- bamVarName(df)
colData <- DataFrame(SAMPLE = as.factor(bamVarName(df, TRUE)),
                     row.names=varNames)
if (!is.null(df$rep)) colData$replicate <- as.factor(df$rep)
if (!is.null(df$stage) & !all(is.na(df$stage))) colData$stage <- as.factor(df$stage)
if (!is.null(df$libtype)) colData$libtype <- as.factor(df$libtype)
if (!is.null(df$condition)) colData$condition <- as.factor(df$condition)
if (!is.null(df$fraction)) colData$fraction <- as.factor(df$fraction)
summarized <- SummarizedExperiment(assays=list(counts=mat), rowRanges=mrna,
                            colData=colData)
saveRDS(summarized, file = file.path(QCfolder(df), "countTable_mrna.rds"))
summarized_dt <- as.data.table(assay(summarized))
summarized_dt <- cbind(GeneID = txNamesToGeneNames(rownames(summarized), df), TxID = rownames(summarized), summarized_dt)
summarized_dt_unique <- summarized_dt[!duplicated(GeneID),]
fwrite(summarized_dt, file = file.path(QCfolder(df), "mrna_expression_all_isoforms.csv"))
fwrite(summarized_dt_unique, file = file.path(QCfolder(df), "mrna_expression_unique_per_gene.csv"))
summarized <- countTable(df, type = "summarized")
summarized_dt <- countTable(df, type = "fpkm")
summarized_dt <- cbind(GeneID = txNamesToGeneNames(rownames(summarized), df), TxID = rownames(summarized), summarized_dt)
summarized_dt_unique <- summarized_dt[!duplicated(GeneID),]
fwrite(summarized_dt, file = file.path(QCfolder(df), "mrna_FPKM_all_isoforms.csv"))
fwrite(summarized_dt_unique, file = file.path(QCfolder(df), "mrna_FPKM_unique_per_gene.csv"))

