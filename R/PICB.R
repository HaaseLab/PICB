#' PICB: piRNA Cluster Builder
#'
#' piRNAs (short for PIWI-interacting RNAs) and their PIWI protein partners play a key role in fertility 
#' and maintaining genome integrity by restricting mobile genetic elements (transposons) in germ cells. 
#' piRNAs originate from genomic regions known as piRNA clusters. The piRNA Cluster Builder (PICB) is 
#' a versatile toolkit designed to identify genomic regions with a high density of piRNAs. It constructs 
#' piRNA clusters through a stepwise integration of unique and multimapping piRNAs and offers wide-ranging 
#' parameter settings, supported by an optimization function that allows users to test different parameter 
#' combinations to tailor the analysis to their specific piRNA system. The output includes extensive metadata 
#' columns, enabling researchers to rank clusters and extract cluster characteristics. 
#'
#' @section Main Functions:
#' The package provides several core functions:
#' \itemize{
#'   \item \code{PICBload}: Load and preprocess BAM files containing piRNA alignments
#'   \item \code{PICBbuild}:  Build piRNA seeds/cores/clusters from alignments
#'   \item \code{PICBoptimize}: Optimize parameters for cluster building
#'   \item \code{PICBstrandanalysis}: Add sense/antisense ratio of unique piRNAs per piRNAcluster
#'   \item \code{PICBannotate}: Annotate GRanges according to a piRNA library
#'   \item \code{PICBloadfasta}: Get SeqInfo object from a fasta file
#'   \item \code{PICBexporttoexcel}: Export cluster object into an Excel file
#'   \item \code{PICBimporttoexcel}: Import cluster object from an Excel file
#' }
#'
#' @section Workflow:
#' A typical PICB workflow consists of:
#' \enumerate{
#'   \item Loading alignments with \code{PICBload}
#'   \item Building clusters with \code{PICBbuild}
#'   \item Optional parameter optimization with \code{PICBoptimize}
#'   \item Optional strand analysis with \code{PICBstrandanalysis}
#'   \item Exporting results with \code{PICBexporttoexcel}
#' }
#' 
#' @examples
#' # 0. Load PICB
#' library(PICB)
#' 
#' # 1. Load Required Genome
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' myGenome <- "BSgenome.Dmelanogaster.UCSC.dm6"
#'
#' # Alternative: Create genome from Seqinfo
#' myGenome2 <- GenomeInfoDb::Seqinfo(
#'     seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"),
#'     seqlengths = c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271, 3667352)
#' )
#'
#' # 2. Load Example Data and Process Alignments
#' bam_file <- system.file("extdata", "Fly_Ov1_chr2L_20To21mb_filtered.bam", 
#'     package = "PICB")
#'
#' myAlignments <- PICBload(
#'     BAMFILE = bam_file,
#'     REFERENCE.GENOME = myGenome,
#'     VERBOSE = FALSE
#' )
#'
#' # 3. Build piRNA Clusters
#' myClusters <- PICBbuild(
#'     IN.ALIGNMENTS = myAlignments,
#'     REFERENCE.GENOME = myGenome,
#'     LIBRARY.SIZE = 12799826,  # Usually calculated automatically
#'     VERBOSITY = 0
#' )$clusters
#'
#' # 4. Optimize Parameters (Optional)
#' parameterExploration <- PICBoptimize(
#'     IN.ALIGNMENTS = myAlignments,
#'     REFERENCE.GENOME = myGenome,
#'     MIN.UNIQUE.ALIGNMENTS.PER.WINDOW = c(1, 2, 3, 4, 5),
#'     LIBRARY.SIZE = 12799826,  # Usually calculated automatically
#'     VERBOSITY = 1
#' )
#'
#' # 5. Perform Strand Analysis
#' myClustersWithStrand <- PICBstrandanalysis(
#'     IN.ALIGNMENTS = myAlignments,
#'     IN.RANGES = myClusters
#' )
#'
#' # 6. Export Clusters
#' PICBexporttoexcel(
#'     IN.RANGES = myClustersWithStrand,
#'     EXCEL.FILE.NAME = "myClusters_demonstration.xlsx"
#' )
#' # 7. Import Ranges
#' importedClusters <- PICBimportfromexcel(
#'     EXCEL.FILE.NAME = system.file("extdata", "myClusters_demonstration.xlsx", package = "PICB")
#' )
#' 
#' @docType package
#' @name PICB
"_PACKAGE"
