#' Indicate Granges dual stranded clusters and bi-directional clusters
#'
#' @param INPUT.GRANGES Granges (seeds/cores/clusters) to annotate
#' @param DUAL.STRAND.MIN.OVERLAP min overlap for dual stranded clusters
#' @param PROD.FOLD max difference in productivity (all 1ry alignments FPKM) for dual-strand clusters
#' @param BIDECTIONAL.MAX.OVERLAP max overlap for bi-directional clusters
#' @param TSS.LEN length of promoter
#'
#' @return Granges object with extra annotation columns
#' @export
#'
#' @examples outputOfPICBbuild$clusters<-PICBstrandanalysis(outputOfPICBbuild$clusters)
PICBstrandanalysis<-function(INPUT.GRANGES, DUAL.STRAND.MIN.OVERLAP=1, PROD.FOLD=3,
                              BIDECTIONAL.MAX.OVERLAP=2000, TSS.LEN=1000){
  FLIP.INPUT.GRANGES<-GenomicRanges::invertStrand(INPUT.GRANGES)
  DT<-data.table::as.data.table(GenomicRanges::findOverlaps(query = INPUT.GRANGES, subject = FLIP.INPUT.GRANGES, minoverlap = DUAL.STRAND.MIN.OVERLAP, type = "any", ignore.strand=FALSE))
  DT$productivity_ratio<-INPUT.GRANGES[DT$queryHits]$all_reads_primary_alignments_FPKM/FLIP.INPUT.GRANGES[DT$subjectHits]$all_reads_primary_alignments_FPKM
  DT2<-DT[(DT$productivity_ratio <= PROD.FOLD) & (DT$productivity_ratio >= 1/PROD.FOLD),]
  INPUT.GRANGES$dual_stranded<-FALSE
  INPUT.GRANGES[unique(DT2$queryHits)]$dual_stranded<-TRUE

  get5prime<-function(INPUT.GRANGES){
    fiveprime<- ifelse(GenomicRanges::strand(INPUT.GRANGES) == "+", GenomicRanges::start(INPUT.GRANGES), GenomicRanges::end(INPUT.GRANGES))
    return(fiveprime)
  }

  get3prime<-function(INPUT.GRANGES){
    threeprime<- ifelse(GenomicRanges::strand(INPUT.GRANGES) == "+", GenomicRanges::end(INPUT.GRANGES), GenomicRanges::start(INPUT.GRANGES))
    return(threeprime)
  }


  a<-as.integer((TSS.LEN-1)/2)
  b<-as.integer((BIDECTIONAL.MAX.OVERLAP+TSS.LEN)/2)
  bidirGRstep1<-INPUT.GRANGES
  bidirGRstep1[(GenomicRanges::start(bidirGRstep1)>=a) ]<-GenomicRanges::resize(bidirGRstep1[GenomicRanges::start(bidirGRstep1)>=a],
                                                                                GenomicRanges::width(bidirGRstep1[GenomicRanges::start(bidirGRstep1)>=a])+a,
                                                                                fix = "end")
  bidirGRstep2<-bidirGRstep1
  bidirGRstep2[GenomicRanges::width(bidirGRstep2)>b]<-GenomicRanges::resize(
    bidirGRstep2[GenomicRanges::width(bidirGRstep2)>b], b, fix = "start")
  bidirGRstep2FLIP<-bidirGRstep2
  GenomicRanges::strand(bidirGRstep2FLIP) <- ifelse(GenomicRanges::strand(bidirGRstep2FLIP) == "+", "-", "+")
  step3hits<-GenomicRanges::findOverlaps(bidirGRstep2,bidirGRstep2FLIP)
  fiveprimeFROM<-get5prime(INPUT.GRANGES[step3hits@from])
  threeprimeFROM<-get3prime(INPUT.GRANGES[step3hits@from])
  fiveprimeTO<-get5prime(INPUT.GRANGES[step3hits@to])
  threeprimeTO<-get3prime(INPUT.GRANGES[step3hits@to])
  step0hits<-GenomicRanges::findOverlaps(INPUT.GRANGES,FLIP.INPUT.GRANGES)
  selectedNumbers<-step3hits@from[(abs(fiveprimeFROM-fiveprimeTO)<abs(fiveprimeFROM-threeprimeTO)) & (abs(fiveprimeFROM-fiveprimeTO)<abs(threeprimeFROM-fiveprimeTO))]
  selectedNumbers<-selectedNumbers[selectedNumbers %in%
                                     step0hits@from[(step0hits@from %in% as.integer(names(table(step0hits@from)[table(step0hits@from)<2]))) &
                                                      (step0hits@to %in% as.integer(names(table(step0hits@to)[table(step0hits@to)<2])))]]
  GenomicRanges::mcols(INPUT.GRANGES)[['bidirections_stranded']]=FALSE
  GenomicRanges::mcols(INPUT.GRANGES)[['bidirections_stranded']][selectedNumbers]=TRUE
  return(INPUT.GRANGES)
}


