# PICB Workflow - R Script
# ✅ R versions 4.2 - 4.4
# ✅ R with IDE R base
# ✅ R with IDE RStudio
# ✅ R with IDE Cursor AI
# The other file in this folder provides the same version of this script but in Jupyter Notebook (R) format.

## Data preparation
# For the demo, we are using a subset of our mapped small RNAs from Drosophila ovaries (SRR7346241). 
# In short, we trimmed the 3’ adaptor and filtered for reads ≥24 nucleotides. Next, we removed abundant cellular RNAs (rRNAs, tRNAs, snRNAs, snoRNAs) by annotation. 
# Finally, we mapped the processed reads to the Dm6 genome using STAR (PMID: 23104886) allowing up to one mismatch and 100 alignments per read. Using samtools (PMID: 19505943), 
# we extracted alignments from chr2L:20000000-21000000 (samtools view -b original.bam chr2L:20000000-21000000 > subset.bam) and generated a subset drosophila genome containing only chromosome 2 (samtools faidx Dm6_genome.fa chr2L chr2R > subset_genome.fa).

# Check successful installation of supported R version 
R.version.string

# Install required R packages from CRAN and Bioconductor
# a. CRAN packages
cran_packages <- c("data.table", "seqinr", "openxlsx", "dplyr")
if (length(setdiff(cran_packages, installed.packages()[,"Package"])) > 0) install.packages(setdiff(cran_packages, installed.packages()[,"Package"]), quiet = TRUE) else message("All required CRAN packages are already installed.")

# b. Install BiocManager
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# c. Install Bioconductor packages
bio_packages <- c("BiocManager", "IRanges", "GenomicRanges", "GenomicAlignments", "Rsamtools", "Biostrings", "GenomeInfoDb", "BSgenome", "rtracklayer")
if (length(setdiff(bio_packages, installed.packages()[,"Package"])) > 0) BiocManager::install(setdiff(bio_packages, installed.packages()[,"Package"]), ask = FALSE) else message("All required Bioconductor packages are already installed.")


# Install PICB - Different options 
#Just run one of the following options. If you need to install the demo BAM and BAI files to follow this demonstration with the provided sample, follow the installation instructions in options c or d.
# a. Through devtools-package
devtools::install_github("HaaseLab/PICB")
# b. Through remotes-package
#remotes::install_github("HaaseLab/PICB")
# c. Download GitHub repository, unzip the file and run: 
#install.packages("../PICB", repos=NULL, type="source")
# d. Clone GitHub repository in terminal
#cd [folder name]
#git clone https://github.com/HaaseLab/PICB.git

# Load PICB
library("PICB") 

#Core-PICB

#Load the reference genome
# Provide reference genome (just run 1 option)
# a. Install and load a BSgenome
#BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
myGenome <- "BSgenome.Dmelanogaster.UCSC.dm6"
# b. Provide a Seqinfo object containing the chromosome names and lengths according to the BAM file.
#myGenome <- GenomeInfoDb::Seqinfo(genome = "dm6")
# c. Provide a fasta file containing the reference genome sequence.
#myGenome <- PICBloadfasta("dm6.fa")

# Load mapped piRNA data
# Run the function PICBload to import and organize the mapped reads in GRanges format.
myAlignments <- PICBload(BAMFILE="demo/Fly_Ov1_chr2L_20To21mb_filtered.bam", REFERENCE.GENOME=myGenome)

# Building piRNA clusters 
# Run PICBbuild to construct piRNA clusters.
myClusters <- PICBbuild(IN.ALIGNMENTS=myAlignments, REFERENCE.GENOME= myGenome, LIBRARY.SIZE = 12799826)$clusters
# Note: Library size usually doesn't need to be set, but we are dealing with a sample that does not contain all reads, so normalization is altered.

# Optimize PICB parameters
parameterExploration <- PICBoptimize(IN.ALIGNMENTS = myAlignments, REFERENCE.GENOME=myGenome, MIN.UNIQUE.ALIGNMENTS.PER.WINDOW=c(1,2,4,5,8,10,15,20), LIBRARY.SIZE = 12799826, VERBOSITY=0)
#Expect a warning due to the modified LIBRARY.SIZE, which does not match the number of reads in IN.ALIGNMENTS.

parameterExploration

# A thorough parameter exploration for a subset of alignments is not particularly useful here. Nevertheless, we provide the code below to visualize the exploration of different values for MIN.UNIQUE.ALIGNMENTS.PER.WINDOW for the reads ("library") and the genome space covered by piRNA clusters.
library(ggplot2)

#Specify parameter to optimize for.
x_column <- "MIN.UNIQUE.ALIGNMENTS.PER.WINDOW" #change parameter to optimize, if applicable

#Determine a scaling factor for the secondary axis
scaling_factor<-max(parameterExploration$fraction.of.library.explained.by.clusters)/max(parameterExploration$clusters.fraction.of.genome.space)

#plot graph
ggplot(parameterExploration, aes(x = .data[[x_column]])) +
    geom_line(aes(y = fraction.of.library.explained.by.clusters*100, color = "piRNAs Explained"), linewidth = 1) +
    geom_point(aes(y = fraction.of.library.explained.by.clusters*100, color = "piRNAs Explained"), size = 3) +
    geom_line(aes(y = clusters.fraction.of.genome.space * scaling_factor*100, color = "Genome Space"), linewidth = 1) +
    geom_point(aes(y = clusters.fraction.of.genome.space * scaling_factor*100, color = "Genome Space"), size = 3) +
    scale_y_continuous(name = "piRNAs Explained by Clusters (%, piRNA sample)", limits = c(0,100), breaks=seq(0,100,10), sec.axis = sec_axis(~./scaling_factor, name = "Total piRNA cluster-length (Genome, %)")) +
    scale_x_reverse(name = stringr::str_to_title(gsub("\\.", " ", x_column)), breaks = parameterExploration[[x_column]], labels = parameterExploration[[x_column]]) +
    scale_color_manual(name = "Metrics", values = c("piRNAs Explained" = "#00a100", "Genome Space" = "black")) +
    theme_classic() + 
    theme(axis.title.y.left = element_text(color = "#00a100"), axis.title.y.right = element_text(color = "black"),legend.position = "top")


# if needed modify paramater in PICBbuild (not needed for PICB-demonstration)


# Ranking piRNA clusters 
myClusters <- myClusters[order(-mcols(myClusters)$all_reads_primary_alignments_FPM)]
mcols(myClusters)$rank_all_reads_primary_alignments_FPM <- rank(-mcols(myClusters)$all_reads_primary_alignments_FPM)


# Exporting results
# RDS format 
saveRDS(myClusters, file="clusters_PICB_demo.RDS")
# Microsoft Excel format
PICBexporttoexcel(RANGES= myClusters, EXCEL.FILE.NAME="clusters_PICB_demo.xlsx")
# GTF format 
myClusters$transcript_id <- paste0("PICB_ rank_", myClusters$rank_all_reads_primary_alignments_FPM)
rtracklayer::export(myClusters, "clusters_PICB_demo.gtf")

