#############################################################################
# This script is for doing DESeq analysis using R
# Adapted from the course on RNA Seq

#############################################################################
# Step - 1: Alignment 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")

library(Rsubread)

# Specify the correct path to the fastq.gz files
fastq.files <- list.files(path = "path/to/fastq.gz/files", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files

# Specify the correct path to the reference genome
buildindex(basename="chr1_mm10", reference="path/to/chr1.fa")

# Specify the correct path to the index file
align(index="path/to/index/file", readfile1=fastq.files)
args(align)
bam.files <- list.files(path = "path/to/BAM/files", pattern = ".BAM$", full.names = TRUE)
bam.files
props <- propmapped(files=bam.files)
props

# Extract quality scores
qs <- qualityScores(filename="path/to/fastq/files", nreads=100)

# Counting
fc <- featureCounts(bam.files, annot.inbuilt="mm10")
fc$stat
dim(fc$counts)

#############################################################################
# Step - 2: TPM preprocessing

# Specify correct paths to the expression data files
tumor_expr <- read.table("path/to/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt", header = TRUE, row.names = 1)
normal_expr <- read.table("path/to/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt", header = TRUE, row.names = 1)

expr <- cbind(tumor_expr, normal_expr)
print(dim(expr))

tumor_samples <- read.table("path/to/GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt", header = TRUE)
normal_samples <- read.table("path/to/GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt", header = TRUE)

print(length(tumor_samples$CancerType)) 
print(length(normal_samples$CancerType))

sample_labels <- c(rep("tumor", ncol(tumor_expr)), rep("normal", ncol(normal_expr)))
print(length(sample_labels))

col_data <- data.frame(row.names = colnames(expr), condition = sample_labels)
print(dim(col_data))

#############################################################################
# Step - 3: Diff Expression Analysis

# Specify correct path to the differential expression dataset
dds <- read.delim("path/to/GSE236050_BMN_RNAseq.txt")
# View(dds) # Uncomment if running interactively

# Remove the first column and set rownames
rownames(dds) <- dds[,1]
dds <- dds[,-1] 

# Optional: Convert rownames to gene symbols
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id', 'entrezgene', 'hgnc_symbol'), mart = mart)

# Create phenotype table
samples <- c("veh.1", "veh.2", "veh.3", "BMN.1", "BMN.2", "BMN.3", "BMN.4")
group <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
g <- data.frame(samples, group)

# Conduct DE analysis using DESeq2
library(DESeq2)
d <- DESeqDataSetFromMatrix(countData = dds, colData = g, design= ~ group)

d <- DESeq(d)
resultsNames(d)
res <- results(d)
res <- as.data.frame(res)
# View(res) # Uncomment if running interactively

# Filter for significant DEGs (padj < 0.05)
res_sig <- subset(res, padj < 0.05)
