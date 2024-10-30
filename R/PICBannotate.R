#' Annotate Granges according to a piRNA library
#'
#' @param INPUT.GRANGES Granges (seeds/cores/clusters) to annotate
#' @param ALIGNMENTS list of alignments from PICBload
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param REPLICATE.NAME name of the replicate. NULL by default.
#' @param LIBRARY.SIZE number of reads in the library. By default computed as number of unique mapping alignments + number of primary multimapping alignments.
#' @param PROVIDE.NON.NORMALIZED provide annotations in non-normalized format. False by default.
#' @param COMPUTE.1U.10A.FRACTIONS for each locus and each alignments type (unique mapping, primary multimapping, secodnary multimapping) compute fraction 1U and 10A containing reads overlapping the locus. Default FALSE.
#'
#' @author Aleksandr Friman
#' @return Granges object with extra annotation columns
#' @export
#'
#' @examples outputOfPICBbuild$seeds<-PICBannotate(outputOfPICBbuild$seeds, outputOfPICBload, REPLICATE.NAME="Replicate1")
PICBannotate<-function(INPUT.GRANGES, ALIGNMENTS, REFERENCE.GENOME = NULL, REPLICATE.NAME = NULL,
                           LIBRARY.SIZE=length(ALIGNMENTS$unique)+length(ALIGNMENTS$multi.primary),
                           PROVIDE.NON.NORMALIZED = FALSE, SEQ.LEVELS.STYLE = "UCSC", COMPUTE.1U.10A.FRACTIONS = FALSE){

  if(is.null(REFERENCE.GENOME)) stop("Please provide REFERENCE.GENOME")
  if (!is.null(REPLICATE.NAME)){
    SUFFIX = paste0(".",REPLICATE.NAME)
  }else{
    SUFFIX = ""
  }

  if (typeof(REFERENCE.GENOME)=="character"){
    SI<-PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)
  }else{
    SI<-REFERENCE.GENOME
  }

  PICBfixCoverage<-function(inCov, SI){
    for (i in 1:length(SI)){
      chrname<-GenomeInfoDb::seqnames(SI)[i]
      if (chrname %in% names(inCov)){
        targetLen<-GenomeInfoDb::seqlengths(SI)[chrname]
        currentLen<-length(inCov[[chrname]])
        deltaLen<-targetLen  - currentLen
        if (deltaLen>0){
          inCov[[chrname]]<-c(inCov[[chrname]],rep(0,deltaLen))
        }
      }
    }
    return(inCov)
  }

  PICBannotateGranges<-function(INPUT.GRANGES, ALIGNMENTS,
                                SI, LIBRARY.SIZE,
                                PROVIDE.NON.NORMALIZED, SEQ.LEVELS.STYLE, COMPUTE.1U.10A.FRACTIONS){
    GenomicRanges::mcols(INPUT.GRANGES)[["width_in_nt"]]<-GenomicRanges::width(INPUT.GRANGES)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads",SUFFIX)]]<-GenomicRanges::countOverlaps(INPUT.GRANGES, ALIGNMENTS$unique)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments",SUFFIX)]]<-GenomicRanges::countOverlaps(INPUT.GRANGES, ALIGNMENTS$multi.primary)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("all_reads_primary_alignments",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads",SUFFIX)]]+GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments",SUFFIX)]]

    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads_FPM",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads",SUFFIX)]]*1e6/LIBRARY.SIZE
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments_FPM",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments",SUFFIX)]]*1e6/LIBRARY.SIZE
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("all_reads_primary_alignments_FPM",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("all_reads_primary_alignments",SUFFIX)]]*1e6/LIBRARY.SIZE

    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads_FPKM",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads_FPM",SUFFIX)]]*1e3/GenomicRanges::width(INPUT.GRANGES)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments_FPKM",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments_FPM",SUFFIX)]]*1e3/GenomicRanges::width(INPUT.GRANGES)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("all_reads_primary_alignments_FPKM",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("all_reads_primary_alignments_FPM",SUFFIX)]]*1e3/GenomicRanges::width(INPUT.GRANGES)
    #coverage
    UniqueMappingCoveragePlus<-GenomicRanges::coverage(ALIGNMENTS$unique[GenomicRanges::strand(ALIGNMENTS$unique)=="+"])
    UniqueMappingCoveragePlus<-PICBfixCoverage(UniqueMappingCoveragePlus, SI)
    UniqueMappingCoverageMinus<-GenomicRanges::coverage(ALIGNMENTS$unique[GenomicRanges::strand(ALIGNMENTS$unique)=="-"])
    UniqueMappingCoverageMinus<-PICBfixCoverage(UniqueMappingCoverageMinus, SI)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("width_covered_by_unique_alignments",SUFFIX)]][as.vector(GenomicRanges::strand(INPUT.GRANGES)=="+")]<-sum(UniqueMappingCoveragePlus[INPUT.GRANGES[GenomicRanges::strand(INPUT.GRANGES)=="+"]]>0)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("width_covered_by_unique_alignments",SUFFIX)]][as.vector(GenomicRanges::strand(INPUT.GRANGES)=="-")]<-sum(UniqueMappingCoverageMinus[INPUT.GRANGES[GenomicRanges::strand(INPUT.GRANGES)=="-"]]>0)
    GenomicRanges::mcols(INPUT.GRANGES)[[paste0("fraction_of_width_covered_by_unique_alignments",SUFFIX)]]<-GenomicRanges::mcols(INPUT.GRANGES)[[paste0("width_covered_by_unique_alignments",SUFFIX)]]/GenomicRanges::width(INPUT.GRANGES)

    if (COMPUTE.1U.10A.FRACTIONS){
      for ( alignmentType in c('unique', 'multi.primary', 'multi.secondary')){
        alignmentsAtWork = ALIGNMENTS[[alignmentType]]
        if ("seq" %in% colnames(GenomicRanges::mcols(alignmentsAtWork))){
          GenomicRanges::mcols(alignmentsAtWork)[['seq']]= as.character(GenomicRanges::mcols(alignmentsAtWork)[['seq']])
          pairsAtWork = GenomicRanges::findOverlaps(INPUT.GRANGES, alignmentsAtWork)
          pairsAtWork.DF = as.data.frame(pairsAtWork)
          pairsAtWork.DF$seq = GenomicRanges::mcols(alignmentsAtWork)[['seq']][pairsAtWork.DF$subjectHits]
          pairsAtWork.DF$oneT = (substring(pairsAtWork.DF$seq, 1, 1) == "T")
          pairsAtWork.DF$tenA = (substring(pairsAtWork.DF$seq, 10, 10) == "A")
          pairsAtWork.DF_by_cluster <- dplyr::group_by(pairsAtWork.DF, queryHits)
          summedDF = dplyr::summarise(pairsAtWork.DF_by_cluster, oneT = sum(oneT), tenA = sum(tenA), normCoef = dplyr::n())
          summedDF$oneT = summedDF$oneT/summedDF$normCoef
          summedDF$tenA = summedDF$tenA/summedDF$normCoef

          GenomicRanges::mcols(INPUT.GRANGES)[[paste0("oneU.frac.",alignmentType,SUFFIX)]] = 0
          GenomicRanges::mcols(INPUT.GRANGES)[[paste0("oneU.frac.",alignmentType,SUFFIX)]][ summedDF$queryHits ] = summedDF$oneT
          GenomicRanges::mcols(INPUT.GRANGES)[[paste0("tenA.frac.",alignmentType,SUFFIX)]] = 0
          GenomicRanges::mcols(INPUT.GRANGES)[[paste0("tenA.frac.",alignmentType,SUFFIX)]][ summedDF$queryHits ] = summedDF$tenA

        }else{
          stop("Alignments '",alignmentType,"' does not contain sequence ('seq' field). Run PICBload with parameter GET.ORIGINAL.SEQUENCE=TRUE to get 'seq' column. Cannot COMPUTE.1U.10A.FRACTIONS.")
        }
      }
    }

    if (! PROVIDE.NON.NORMALIZED){
      GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_reads",SUFFIX)]] <- NULL
      GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_sequences",SUFFIX)]]  <- NULL
      GenomicRanges::mcols(INPUT.GRANGES)[[paste0("multimapping_reads_primary_alignments",SUFFIX)]] <- NULL
      GenomicRanges::mcols(INPUT.GRANGES)[[paste0("all_reads_primary_alignments",SUFFIX)]] <- NULL
      GenomicRanges::mcols(INPUT.GRANGES)[[paste0("width_covered_by_unique_alignments",SUFFIX)]] <- NULL
    }else{
      GenomicRanges::mcols(INPUT.GRANGES)[[paste0("uniq_sequences",SUFFIX)]]<-GenomicRanges::countOverlaps(INPUT.GRANGES, unique(ALIGNMENTS$unique))
    }
    return(INPUT.GRANGES)
  }
  #identification of cluster types
  PICBannotateTypesOfClusters<-function(INPUT.CLUSTERS, INPUT.CORES, SI){
    INPUT.CLUSTERS$cores_intersected<-GenomicRanges::countOverlaps(INPUT.CLUSTERS, INPUT.CORES)
    INPUT.CORES.plus<-INPUT.CORES[GenomicRanges::strand(INPUT.CORES)=="+"]
    INPUT.CORES.minus<-INPUT.CORES[GenomicRanges::strand(INPUT.CORES)=="-"]
    CORES.COV.PLUS<-GenomicRanges::coverage(INPUT.CORES.plus)
    CORES.COV.PLUS<-PICBfixCoverage(inCov = CORES.COV.PLUS, SI = SI)
    CORES.COV.MINUS<-GenomicRanges::coverage(INPUT.CORES.minus)
    CORES.COV.MINUS<-PICBfixCoverage(inCov = CORES.COV.MINUS, SI = SI)
    INPUT.CLUSTERS.plus<-INPUT.CLUSTERS[GenomicRanges::strand(INPUT.CLUSTERS)=="+"]
    INPUT.CLUSTERS.minus<-INPUT.CLUSTERS[GenomicRanges::strand(INPUT.CLUSTERS)=="-"]
    INPUT.CLUSTERS.plus$core_coverage<-sum(CORES.COV.PLUS[INPUT.CLUSTERS.plus])
    INPUT.CLUSTERS.minus$core_coverage<-sum(CORES.COV.MINUS[INPUT.CLUSTERS.minus])
    INPUT.CLUSTERS<-c(INPUT.CLUSTERS.plus, INPUT.CLUSTERS.minus)
    INPUT.CLUSTERS$core_coverage_fraction<-INPUT.CLUSTERS$core_coverage/GenomicRanges::width(INPUT.CLUSTERS)
    GenomicRanges::mcols(INPUT.CLUSTERS)[[paste0("type",SUFFIX)]]=NA
    GenomicRanges::mcols(INPUT.CLUSTERS)[[paste0("type",SUFFIX)]][INPUT.CLUSTERS$cores_intersected>1]="MultiCore"
    GenomicRanges::mcols(INPUT.CLUSTERS)[[paste0("type",SUFFIX)]][(INPUT.CLUSTERS$cores_intersected==1) & (INPUT.CLUSTERS$core_coverage_fraction<1)]="ExtendedCore"
    GenomicRanges::mcols(INPUT.CLUSTERS)[[paste0("type",SUFFIX)]][(INPUT.CLUSTERS$cores_intersected==1) & (INPUT.CLUSTERS$core_coverage_fraction==1)]="SingleCore"
    #cleaning up
    INPUT.CLUSTERS$cores_intersected<-NULL
    INPUT.CLUSTERS$core_coverage<-NULL
    INPUT.CLUSTERS$core_coverage_fraction<-NULL
    return(INPUT.CLUSTERS)
  }
  #processing the input genomic intervals
  if (typeof(INPUT.GRANGES)=="list"){
    for (t in intersect(c(uniqueonly, uniqueandprimary, allalignments),names(INPUT.GRANGES))){
      if (length(INPUT.GRANGES[[t]])>0){
        INPUT.GRANGES[[t]]<-PICBannotateGranges(INPUT.GRANGES[[t]], ALIGNMENTS,
                                              SI, LIBRARY.SIZE,
                                              PROVIDE.NON.NORMALIZED, SEQ.LEVELS.STYLE, COMPUTE.1U.10A.FRACTIONS = COMPUTE.1U.10A.FRACTIONS)
      }
    }
    if (allalignments %in% names(INPUT.GRANGES) && length(INPUT.GRANGES[[allalignments]])>0){
      INPUT.GRANGES[[allalignments]]<-PICBannotateTypesOfClusters(INPUT.GRANGES[[allalignments]], INPUT.GRANGES[[uniqueandprimary]], SI)
    }
  }else{
    if (length(INPUT.GRANGES)>0){
      INPUT.GRANGES<-PICBannotateGranges(INPUT.GRANGES, ALIGNMENTS,
                        SI, LIBRARY.SIZE,
                        PROVIDE.NON.NORMALIZED, SEQ.LEVELS.STYLE, COMPUTE.1U.10A.FRACTIONS = COMPUTE.1U.10A.FRACTIONS)
    }
  }

  return(INPUT.GRANGES)
}
