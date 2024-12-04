# =======================================
# Script Name: PICB_demo.R
# Purpose: Demonstration of the PICB workflow using a processed subset 
#          of small RNA data from Drosophila ovaries.
# Authors: PICB Team of the Haase Lab
# Date: 2024-12-04
#
# Description:
# This script demonstrates analysis on a reduced dataset, 
# Fly_Ov1_chr2L_20To21mb_filtered.bam, created from small RNA reads 
# aligned to Drosophila melanogaster genome (Dm6). The data preparation 
# includes trimming, filtering, mapping, and downsampling as detailed below.
#
# Compatible Environments:
# - ✅ R versions 4.2 - 4.4
# - ✅ R with R base
# - ✅ R with IDE RStudio
# - ✅ R with IDE Cursor AI
#
# =======================================

# Data Preparation Context:
# - Processed data: Fly_Ov1_chr2L_20To21mb_filtered.bam
# - Details:
#   - Trimmed 3' adaptors, filtered reads ≥24 nt.
#   - Removed abundant cellular RNAs (rRNAs, tRNAs, snRNAs, snoRNAs).
#   - Mapped reads to Dm6 genome (chr2L:20000000-21000000) using STAR.
#   - Retained essential tags (NH, NM) and compressed BAM file.
#   - Downsampled to 22% of reads for a smaller, manageable dataset.
#   - Indexed and renamed the file.

# Disclaimer:
# The processed data is intended for demonstration purposes only and 
# is not representative of biologically relevant samples due to 
# extensive downsampling and regional restriction.

# Additional Notes:
# The Jupyter Notebook version of this script is available in the 
# corresponding folder. Data preparation commands can be found on the 
# [PICB GitHub page](https://github.com/HaaseLab/PICB).

# =======================================

# Check R version compatibility
message("Using R version: ", R.version.string)

# Install required packages
# =======================================
# CRAN packages
cran_packages <- c("data.table", "seqinr", "openxlsx", "dplyr")
if (length(setdiff(cran_packages, installed.packages()[,"Package"])) > 0) 
    install.packages(setdiff(cran_packages, installed.packages()[,"Package"]), 
    quiet = TRUE) else 
    message("All required CRAN packages are already installed.")
# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bio_packages <- c("BiocManager", "IRanges", "GenomicRanges", "GenomicAlignments", "Rsamtools", "Biostrings", "GenomeInfoDb", "BSgenome", "rtracklayer")
if (length(setdiff(bio_packages, installed.packages()[,"Package"])) > 0) BiocManager::install(setdiff(bio_packages, installed.packages()[,"Package"]), ask = FALSE) else message("All required Bioconductor packages are already installed.")


# Install PICB
# =======================================
# Select one of the following installation methods for the PICB package:
# Uncomment the preferred method.
# devtools::install_github("HaaseLab/PICB")
# remotes::install_github("HaaseLab/PICB")
# install.packages("../PICB", repos = NULL, type = "source")
# Use terminal to clone and build if required:
# git clone https://github.com/HaaseLab/PICB.git

# Load PICB package
library("PICB") 

#Core-PICB

# Load reference genome
# =======================================
# Provide the reference genome. Select one option:
# a. Install and load a BSgenome
# BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
library("BSgenome.Dmelanogaster.UCSC.dm6")
myGenome <- "BSgenome.Dmelanogaster.UCSC.dm6"
# b. Provide a Seqinfo object containing the chromosome names and lengths according to the BAM file.
# myGenome <- GenomeInfoDb::Seqinfo(
#     seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"),
#     seqlengths = c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271, 3667352)
# )
# c. Provide supported genome as Seqinfo object.
# myGenome <- GenomeInfoDb::Seqinfo(genome = "dm6")
# d. Provide a fasta file containing the reference genome sequence.
# myGenome <- PICBloadfasta("dm6.fa")


# Load mapped piRNA data
# =======================================
# Load the demo BAM file and process it into GRanges format
bam_file <- system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", package = "PICB")
myAlignments <- PICBload(
    BAMFILE = bam_file,
    REFERENCE.GENOME = myGenome
)


# Build piRNA clusters
# =======================================
# Run PICBbuild to construct piRNA clusters.
myClusters <- PICBbuild(
    IN.ALIGNMENTS = myAlignments,
    REFERENCE.GENOME = myGenome,
    LIBRARY.SIZE = 12799826
)$clusters
# Note: Library size usually doesn't need to be set, but we are dealing with a sample that does not contain all reads, so normalization is altered.

# Optimize PICB parameters
# =======================================
parameterExploration <- PICBoptimize(
    IN.ALIGNMENTS = myAlignments,
    REFERENCE.GENOME = myGenome,
    MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2, 3, 4, 5),
    LIBRARY.SIZE = 12799826,
    VERBOSITY = 2
)
# Note: Expect a warning due to the modified LIBRARY.SIZE, which does not match the number of reads in IN.ALIGNMENTS.

parameterExploration

# A thorough parameter exploration for a subset of alignments is not particularly useful here. Nevertheless, we provide the code below to visualize the exploration of different values for MIN.UNIQUE.ALIGNMENTS.PER.WINDOW for the reads ("library") and the genome space covered by piRNA clusters.
library(ggplot2)

#Specify parameter to optimize for.
x_column <- "MIN.UNIQUE.ALIGNMENTS.PER.WINDOW" #change parameter to optimize, if applicable

scaling_factor <- max(parameterExploration$fraction.of.library.explained.by.clusters) / max(parameterExploration$clusters.fraction.of.genome.space)
ggplot(parameterExploration, aes(x = .data[[x_column]])) +
    geom_line(aes(y = fraction.of.library.explained.by.clusters * 100, color = "piRNAs Explained"), linewidth = 1) +
    geom_point(aes(y = fraction.of.library.explained.by.clusters * 100, color = "piRNAs Explained"), size = 3) +
    geom_line(aes(y = clusters.fraction.of.genome.space * scaling_factor * 100, color = "Genome Space"), linewidth = 1) +
    geom_point(aes(y = clusters.fraction.of.genome.space * scaling_factor * 100, color = "Genome Space"), size = 3) +
    scale_y_continuous(name = "piRNAs Explained by Clusters (%, piRNA sample)", sec.axis = sec_axis(~ . / scaling_factor, name = "Total piRNA cluster-length (Genome, %)")) +
    scale_x_reverse(name = paste0("Parameter chosen: ", x_column), breaks = parameterExploration[[x_column]], labels = parameterExploration[[x_column]]) +
    scale_color_manual(name = "Metrics", values = c("piRNAs Explained" = "#00a100", "Genome Space" = "black")) +
    theme_classic() +
    theme(axis.title.y.left = element_text(color = "#00a100"), axis.title.y.right = element_text(color = "black"), legend.position = "top")

# if needed modify paramater in PICBbuild (not needed for PICB-demonstration)


# Rank piRNA clusters
# =======================================
# one option is to rank by primary alignments FPM:
myClusters <- myClusters[order(-mcols(myClusters)$all_reads_primary_alignments_FPM)]
mcols(myClusters)$rank_all_reads_primary_alignments_FPM <- rank(-mcols(myClusters)$all_reads_primary_alignments_FPM)


# Export results
# =======================================
# RDS format 
saveRDS(myClusters, file="clusters_PICB_demo.RDS")
# Microsoft Excel format
PICBexporttoexcel(IN.RANGES = myClusters, EXCEL.FILE.NAME = "clusters_PICB_demo.xlsx")
# GTF format 
myClusters$transcript_id <- paste0("PICB_ rank_", myClusters$rank_all_reads_primary_alignments_FPM)
rtracklayer::export(myClusters, "clusters_PICB_demo.gtf")

