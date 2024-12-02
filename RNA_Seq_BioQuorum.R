#############################################################################
# This script is for doing DESeq analysis using R
# Adapted from the course on RNA Seq

#############################################################################
# Step - 1: Alignment 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")

library(Rsubread)

fastq.files <- list.files(path = "path to the fastq.gz files", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files

buildindex(basename="chr1_mm10",reference="chr1.fa")

align(index="path to index file",readfile1=fastq.files)
args(align)
bam.files <- list.files(path = "path to BAM files", pattern = ".BAM$", full.names = TRUE)
bam.files
props <- propmapped(files=bam.files)
props

# Extract quality scores
qs <- qualityScores(filename="path to the fastq files",nreads=100)

# Counting
fc <- featureCounts(bam.files, annot.inbuilt="mm10")
fc$stat
dim(fc$counts)

#############################################################################
# Step - 2: TPM preprocessing

tumor_expr <- read.table("~/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt", header = TRUE, row.names = 1)
normal_expr <- read.table("~/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt", header = TRUE, row.names = 1)

expr <- cbind(tumor_expr, normal_expr)
print(dim(expr))

tumor_samples <- read.table("~/GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt", header = TRUE)
normal_samples <- read.table("~/GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt", header = TRUE)

print(length(tumor_samples$CancerType)) 
print(length(normal_samples$CancerType))

sample_labels <- c(rep("tumor", ncol(tumor_expr)), rep("normal", ncol(normal_expr)))
print(length(sample_labels))

col_data <- data.frame(row.names = colnames(expr), condition = sample_labels)
print(dim(col_data))

#############################################################################
# Step - 3: Diff Expression Analysis

#url: GSE236050
dds <- read.delim("~/GSE236050_BMN_RNAseq.txt")
View(dds) #to see the file - you will see that the column 1 is the ENSEMBL ID
## remove the first column
rownames(dds) <- dds[,c(1)]
dds <- dds[,-c(1)] ## now rownames are the gene id's 

## convert rownames to gene symbols (Optional)

library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id',
                              'entrezgene', 'hgnc_symbol'),mart = mart)

## create your phenotype table
##First three columns vs the last four columns

samples =c("veh.1","veh.2","veh.3","BMN.1","BMN.2","BMN.3","BMN.4")
group <- c("Control","Control","Control","Treatment","Treatment","Treatment","Treatment")
g <- data.frame(samples,group)


## conduct DEG using DESeq2 since we have the raw gene counts
library(DESeq2)
d <- DESeqDataSetFromMatrix(countData = dds,
                              colData = g,
                              design= ~ group)

d <- DESeq(d)
resultsNames(d)
res <- results(d)
res <- as.data.frame(res)
View(res)

##Filter for significant DEG (padj < 0.05)

res_sig <- res %>% filter(padj < 0.05)
