#' Build piRNA seeds/cores/clusters from alignments
#'
#' @param IN.ALIGNMENTS list of alignments from PICBload
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH width of sliding window for unique mappers. 350 nt by default
#' @param UNIQUEMAPPERS.SLIDING.WINDOW.STEP step of sliding windows for unique mappers. width/10 by default
#' @param PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH width of sliding window for primary multimapping alignments. 350 nt by default
#' @param PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.STEP step of sliding windows for primary multimapping alignments. width/10 by default
#' @param SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH width of sliding window for secondary multimapping alignments. 1000 nt by default
#' @param SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.STEP step of sliding windows for secondary multimapping alignments. width/10 by default
#' @param LIBRARY.SIZE number of reads in the library. By default computed as number of unique mapping alignments + number of primary multimapping alignments.
#' @param MIN.UNIQUE.ALIGNMENTS.PER.WINDOW absolute number of unique mapping alignments per window to call it. By default computed as 2 FPKM.
#' @param MIN.UNIQUE.SEQUENCES.PER.WINDOW absolute number of unique mapping sequences per window to call it. By default computed as width/50.
#' @param MIN.PRIMARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW absolute number of primary multimapping alignments per window to call it. By default computed as 4 FPKM.
#' @param MIN.SECONDARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW absolute number of secondary multimapping alignments per window to call it. By default computed as 0.2 FPKM.
#' @param MIN.SEED.LENGTH minimum length of a seed. By default computed as 2x unique mapper window size + 100.
#' @param MIN.COVERED.SEED.LENGTH minimum number of seed nucleotides covered by unique mappers. 0 by default.
#' @param THRESHOLD.SEEDS.GAP minimum gap between seeds to not merge them. 0 by default.
#' @param THRESHOLD.CORES.GAP minimum gap between cores to not merge them. 0 by default.
#' @param THRESHOLD.CLUSTERS.GAP minimum gap between clusters to not merge them. 0 by default.
#' @param SEQ.LEVELS.STYLE style of chromosome names for BSgenome. "UCSC" by default.
#' @param MIN.OVERLAP minimum overlap between seeds and cores, as well as between cores and clusters 5 nt by default.
#' @param PROVIDE.NON.NORMALIZED include non-normalized to the library size statistics in the output annotations
#' @param COMPUTE.1U.10A.FRACTIONS for each locus and each alignments type (unique mapping, primary multimapping, secodnary multimapping) compute fraction 1U and 10A containing reads overlapping the locus. Default FALSE.
#' @param VERBOSITY verbosity level 0/1/2/3. 2 by default.
#'
#' @author Pavol Genzor
#' @author Daniel Stoyko
#' @author Aleksandr Friman
#' @return list of annotated Granges objects named "seeds" for seeds,
#' "cores" for cores,
#' "clusters" for clusters
#' @export
#'
#' @examples 
#' #'
#' @examples 
#' library(BSgenome.Dmelanogaster.UCSC.dm6)
#' myAlignmentsFromPICBload <- PICBload(
#'   BAMFILE=system.file("extdata", "Fly_Ov1_chr2L_20To21mb.bam", package = "PICB"), 
#'   REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6", 
#'   VERBOSE=FALSE
#' )
#' 
#' outputOfPICBbuild <- PICBbuild(
#'   IN.ALIGNMENTS = myAlignmentsFromPICBload,
#'   REFERENCE.GENOME = "BSgenome.Dmelanogaster.UCSC.dm6", 
#'   VERBOSITY=0
#' )
PICBbuild <-
  function(
    ## MAIN
    IN.ALIGNMENTS,
    REFERENCE.GENOME,
    ## OPTIONS
    UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH=350,
    UNIQUEMAPPERS.SLIDING.WINDOW.STEP=round(UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH/10,0),
    PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH=350,
    PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.STEP=round(PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH/10,0),
    SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH=1000,
    SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.STEP=round(SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH/10,0),
    LIBRARY.SIZE=length(IN.ALIGNMENTS$unique)+length(IN.ALIGNMENTS$multi.primary),
    MIN.UNIQUE.ALIGNMENTS.PER.WINDOW=2*(UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH/1e3)*(LIBRARY.SIZE/1e6),
    MIN.UNIQUE.SEQUENCES.PER.WINDOW=min(MIN.UNIQUE.ALIGNMENTS.PER.WINDOW, round(UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH/50,0)),
    MIN.PRIMARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW=4*(PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH/1e3)*(LIBRARY.SIZE/1e6),
    MIN.SECONDARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW=0.2*(SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH/1e3)*(LIBRARY.SIZE/1e6),
    MIN.SEED.LENGTH=2*UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH+100,
    MIN.COVERED.SEED.LENGTH=0,
    ## CLUSTER FILTERS
    THRESHOLD.SEEDS.GAP=0, THRESHOLD.CORES.GAP=0, THRESHOLD.CLUSTERS.GAP=0,
    ## EXTRA OPTIONS
    SEQ.LEVELS.STYLE="UCSC",
    MIN.OVERLAP=5,
    PROVIDE.NON.NORMALIZED = FALSE,
    COMPUTE.1U.10A.FRACTIONS = FALSE,
    VERBOSITY=2){


    ## Authors: Pavol Genzor, Alex Friman, Daniel Stoyko
    ## Usage: Build piRNA clusters using piRNA sequencing data

    ## OUTPUT FOR DEBUGGING PURPOSES
    if (VERBOSITY>2) {
      tmpENV<-as.list(environment())
      tmpENV$IN.ALIGNMENTS<-NULL
      message(tmpENV)
    }
    ## Check input

    if(is.null(IN.ALIGNMENTS)) stop("Please provide IN.ALIGNMENTS !")
    if(is.null(REFERENCE.GENOME)) stop("Please provide REFERENCE.GENOME !")
    #check if IN.ALIGNMENTS contains columns "unique" and "multi.primary"
    typeAlignments <- names(IN.ALIGNMENTS)
    if(!"unique" %in% typeAlignments) stop("IN.ALIGNMENTS must contain 'unique' column! Use PICBload() with IS.SECONDARY.ALIGNMENT = NA (default) to load piRNAs. \n")
    if(!"multi.primary" %in% typeAlignments) stop("IN.ALIGNMENTS must contain 'multi.primary' column! Use PICBload() with IS.SECONDARY.ALIGNMENT = NA (default) to load piRNAs. \n")
    if(!"multi.secondary" %in% typeAlignments) warning("IN.ALIGNMENTS does not contain secondary multimappers (column 'multi.secondary')! We recommend including secondary multimappers; please set IS.SECONDARY.ALIGNMENT = NA (default) in PICBload. \n\tContinuing without secondary multimappers...\n")

    for (columnName in typeAlignments){
      if(!"NH" %in% colnames(GenomicRanges::mcols(IN.ALIGNMENTS[[columnName]]))) stop("The IN.ALIGNMENTS must contain NH information !")
    }

    ##
    if (VERBOSITY>0) message("PICB v", utils::packageVersion("PICB")," Starting ... ")
    ##

    ## KEEP STANDARD CHROMOSOMES
    if (typeof(REFERENCE.GENOME)=="character"){
      SI<-PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)
    }else{
      SI<-REFERENCE.GENOME
    }
    ##
    ## Clean data
    if (VERBOSITY>1) message("\n\tKeeping standard linear chromosomes")

    for (columnName in typeAlignments){
      #trying to change the seqlevels style.
      #read https://github.com/Bioconductor/GenomeInfoDb/blob/devel/inst/extdata/dataFiles/README
      #for more info
      tryCatch( expr = {GenomeInfoDb::seqlevelsStyle(IN.ALIGNMENTS[[columnName]]) <- SEQ.LEVELS.STYLE},
                error = function(e){
                  message(e)
                  message("Failed to change SEQ.LEVELS.STYLE")
                  message("read https://github.com/Bioconductor/GenomeInfoDb/blob/devel/inst/extdata/dataFiles/README")
                  message("Continuing hoping for the best")
                  })

      KEEP.SEQLEVELS <- GenomeInfoDb::seqlevels(SI)
      KEEP.SEQLEVELS <- KEEP.SEQLEVELS[KEEP.SEQLEVELS %in% GenomeInfoDb::seqlevels(IN.ALIGNMENTS[[columnName]])]
      IN.ALIGNMENTS[[columnName]] <- GenomeInfoDb::keepSeqlevels(x = IN.ALIGNMENTS[[columnName]], value = KEEP.SEQLEVELS, pruning.mode = "coarse")
    }

    WGRU<-IN.ALIGNMENTS$unique
    WGRMP<-IN.ALIGNMENTS$multi.primary

    ##
    if (VERBOSITY>0) message("\nBuilding ... STEP 1... Searching using windows\n")
    ##

    ##
    ## Make genome sliding window for unique mappers

    if (VERBOSITY>1) message("\tSliding window analysis")
    if (VERBOSITY>1) message("\t\tUNIQUE MAPPERS\n\t\tWINDOW: ",UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH,"\n\t\tSTEP: ", UNIQUEMAPPERS.SLIDING.WINDOW.STEP)
    AG.gr <- GenomicRanges::GRanges(data.table::data.table("chr"=GenomicRanges::seqnames(SI),"start"=1,"end"=GenomeInfoDb::seqlengths(SI),"strand"="*"))
    AG.sw.uniq <- unlist(GenomicRanges::slidingWindows(x = AG.gr, width = UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH, step = UNIQUEMAPPERS.SLIDING.WINDOW.STEP))

    ## Make genome sliding window for multi mappers

    if (VERBOSITY>1) message("\t\tPRIMARY MULTI MAPPERS\n\t\tWINDOW: ",PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH,"\n\t\tSTEP: ", PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.STEP)
    #AG.gr <- GRanges(data.table("chr"=seqnames(SI),"start"=1,"end"=seqlengths(SI),"strand"="*"))
    AG.sw.primary.mult <- unlist(GenomicRanges::slidingWindows(x = AG.gr, width = PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH, step = PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.STEP))
    if (VERBOSITY>1 && "multi.secondary" %in% typeAlignments) message("\t\tSECONDARY MULTI MAPPERS\n\t\tWINDOW: ",SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH,"\n\t\tSTEP: ", SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.STEP)

    ##
    ## Counting piRNAs per window: Uniq and all multi
    GenomicRanges::mcols(AG.sw.uniq)[["uniq_piRNA_plus"]] <- GenomicRanges::countOverlaps(query = AG.sw.uniq, subject = WGRU[GenomicRanges::strand(WGRU) == "+"])
    GenomicRanges::mcols(AG.sw.uniq)[["uniq_piRNA_minus"]] <- GenomicRanges::countOverlaps(query = AG.sw.uniq, subject = WGRU[GenomicRanges::strand(WGRU) == "-"])
    GenomicRanges::mcols(AG.sw.uniq)[["uniq_intervals_plus"]] <- GenomicRanges::countOverlaps(query = AG.sw.uniq, subject = unique(WGRU[GenomicRanges::strand(WGRU) == "+"]))
    GenomicRanges::mcols(AG.sw.uniq)[["uniq_intervals_minus"]] <- GenomicRanges::countOverlaps(query = AG.sw.uniq, subject = unique(WGRU[GenomicRanges::strand(WGRU) == "-"]))
    GenomicRanges::mcols(AG.sw.primary.mult)[["primary_mult_piRNA_plus"]] <- GenomicRanges::countOverlaps(query = AG.sw.primary.mult, subject = WGRMP[GenomicRanges::strand(WGRMP) == "+"])
    GenomicRanges::mcols(AG.sw.primary.mult)[["primary_mult_piRNA_minus"]] <- GenomicRanges::countOverlaps(query = AG.sw.primary.mult, subject = WGRMP[GenomicRanges::strand(WGRMP) == "-"])
    if ("multi.secondary" %in% typeAlignments){
      AG.sw.secondary.mult <- unlist(GenomicRanges::slidingWindows(x = AG.gr, width = SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH, step = SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.STEP))
      GenomicRanges::mcols(AG.sw.secondary.mult)[["secondary_mult_piRNA_plus"]] <- GenomicRanges::countOverlaps(query = AG.sw.secondary.mult, subject = IN.ALIGNMENTS$multi.secondary[GenomicRanges::strand(IN.ALIGNMENTS$multi.secondary) == "+"])
      GenomicRanges::mcols(AG.sw.secondary.mult)[["secondary_mult_piRNA_minus"]] <- GenomicRanges::countOverlaps(query = AG.sw.secondary.mult, subject = IN.ALIGNMENTS$multi.secondary[GenomicRanges::strand(IN.ALIGNMENTS$multi.secondary) == "-"])
    }

    ##
    ## Filter by minimum reads per window and reduce to new range
    if (VERBOSITY>1) message("\t\tMIN.UNIQUE.ALIGNMENTS.PER.WINDOW >= ",MIN.UNIQUE.ALIGNMENTS.PER.WINDOW)
    sw.UNIQ.PLUS <- AG.sw.uniq[(GenomicRanges::mcols(AG.sw.uniq)[["uniq_piRNA_plus"]] >= MIN.UNIQUE.ALIGNMENTS.PER.WINDOW) & (GenomicRanges::mcols(AG.sw.uniq)[["uniq_intervals_plus"]] >= MIN.UNIQUE.SEQUENCES.PER.WINDOW)]
    sw.UNIQ.MINUS <- AG.sw.uniq[(GenomicRanges::mcols(AG.sw.uniq)[["uniq_piRNA_minus"]] >= MIN.UNIQUE.ALIGNMENTS.PER.WINDOW) & (GenomicRanges::mcols(AG.sw.uniq)[["uniq_intervals_minus"]] >= MIN.UNIQUE.SEQUENCES.PER.WINDOW)]
    AG.sw.uniq<-NULL # clearing memory
    GenomicRanges::strand(sw.UNIQ.PLUS) <- "+"; GenomicRanges::strand(sw.UNIQ.MINUS) <- "-"
    sw.UNIQ.RED <- GenomicRanges::sort.GenomicRanges(GenomicRanges::reduce(c(sw.UNIQ.PLUS,sw.UNIQ.MINUS)))
    sw.UNIQ.PLUS<-NULL #clearing memory
    ## Removing gaps between SEEDS
    if (THRESHOLD.SEEDS.GAP>0){
      GR.GAPS <- GenomicRanges::gaps(sw.UNIQ.RED)
      GR.GAPS <- GenomicRanges::sort.GenomicRanges(GR.GAPS)
      GR.GAPS.WIDER.THAN.MIN <- GR.GAPS[GenomicRanges::width(GR.GAPS) >= THRESHOLD.SEEDS.GAP]
      sw.UNIQ.RED <- GenomicRanges::gaps(GR.GAPS.WIDER.THAN.MIN) #gaps of gaps
      sw.UNIQ.RED <- GenomicRanges::sort.GenomicRanges(sw.UNIQ.RED)
    }
    sw.UNIQ.RED<-GenomicRanges::reduce(sw.UNIQ.RED[GenomicRanges::width(sw.UNIQ.RED)>=MIN.SEED.LENGTH])
    if (VERBOSITY>1) message("Annotating seeds")
    BM.SEEDS<-PICBannotate(sw.UNIQ.RED,IN.ALIGNMENTS, REFERENCE.GENOME = REFERENCE.GENOME, LIBRARY.SIZE = LIBRARY.SIZE,
                           SEQ.LEVELS.STYLE = SEQ.LEVELS.STYLE, PROVIDE.NON.NORMALIZED = TRUE)
    if (VERBOSITY>1) message("Removing seed with unique mapping coverage less than", MIN.COVERED.SEED.LENGTH, "nt")
    BM.SEEDS<-BM.SEEDS[GenomicRanges::mcols(BM.SEEDS)[['width_covered_by_unique_alignments']]>=MIN.COVERED.SEED.LENGTH]
    outputList<-list()
    outputList[[uniqueonly]]<-BM.SEEDS

    if (VERBOSITY>1) message("\t\tMIN.PRIMARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW >= ",MIN.PRIMARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW)
    sw.PRIMARY.MULT.PLUS <- AG.sw.primary.mult[GenomicRanges::mcols(AG.sw.primary.mult)[["primary_mult_piRNA_plus"]] >= MIN.PRIMARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW]
    sw.PRIMARY.MULT.MINUS <- AG.sw.primary.mult[GenomicRanges::mcols(AG.sw.primary.mult)[["primary_mult_piRNA_minus"]] >= MIN.PRIMARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW]
    GenomicRanges::strand(sw.PRIMARY.MULT.PLUS) <- "+"; GenomicRanges::strand(sw.PRIMARY.MULT.MINUS) <- "-"
    sw.PRIMARY.MULT.PLUS.RED<-GenomicRanges::reduce(sw.PRIMARY.MULT.PLUS)
    sw.PRIMARY.MULT.MINUS.RED<-GenomicRanges::reduce(sw.PRIMARY.MULT.MINUS)
    sw.PRIMARY.MULT.PLUS.ANCHORED <- IRanges::subsetByOverlaps(x = sw.PRIMARY.MULT.PLUS.RED, ranges = BM.SEEDS, minoverlap = MIN.OVERLAP)
    sw.PRIMARY.MULT.MINUS.ANCHORED <- IRanges::subsetByOverlaps(x = sw.PRIMARY.MULT.MINUS.RED, ranges = BM.SEEDS, minoverlap = MIN.OVERLAP)

    sw.PRIMARY.MULT.RED.ANCHORED <- GenomicRanges::sort.GenomicRanges(GenomicRanges::reduce(c(sw.PRIMARY.MULT.PLUS.ANCHORED,sw.PRIMARY.MULT.MINUS.ANCHORED)))
    sw.PRIMARY.MULT.RED.ANCHORED <- GenomicRanges::sort.GenomicRanges(GenomicRanges::reduce(c(sw.PRIMARY.MULT.RED.ANCHORED,BM.SEEDS)))

    ONE.WINDOW.GR <- c(BM.SEEDS,sw.PRIMARY.MULT.RED.ANCHORED)
    BM.CORES <- GenomicRanges::reduce(ONE.WINDOW.GR)

    ##
    ## Removing gaps between cores
    if (THRESHOLD.CORES.GAP>0){
      GR.GAPS <- GenomicRanges::gaps(BM.CORES)
      GR.GAPS <- GenomicRanges::sort.GenomicRanges(GR.GAPS)
      GR.GAPS.WIDER.THAN.MIN <- GR.GAPS[GenomicRanges::width(GR.GAPS) >= THRESHOLD.CORES.GAP]
      BM.CORES <- GenomicRanges::gaps(GR.GAPS.WIDER.THAN.MIN) #gaps of gaps
      BM.CORES <- GenomicRanges::sort.GenomicRanges(BM.CORES)
    }

    ##Regions
    if ("multi.secondary" %in% typeAlignments){
      if (VERBOSITY>1) message("\t\tMIN.SECONDARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW >= ",MIN.SECONDARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW)

      sw.SECONDARY.MULT.PLUS <- AG.sw.secondary.mult[GenomicRanges::mcols(AG.sw.secondary.mult)[["secondary_mult_piRNA_plus"]] >= MIN.SECONDARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW]
      sw.SECONDARY.MULT.MINUS <- AG.sw.secondary.mult[GenomicRanges::mcols(AG.sw.secondary.mult)[["secondary_mult_piRNA_minus"]] >= MIN.SECONDARY.MULTIMAPPING.ALIGNMENTS.PER.WINDOW]
      GenomicRanges::strand(sw.SECONDARY.MULT.PLUS) <- "+"; GenomicRanges::strand(sw.SECONDARY.MULT.MINUS) <- "-"
      sw.SECONDARY.MULT.PLUS.RED<-GenomicRanges::reduce(sw.SECONDARY.MULT.PLUS)
      sw.SECONDARY.MULT.MINUS.RED<-GenomicRanges::reduce(sw.SECONDARY.MULT.MINUS)
      sw.SECONDARY.MULT.PLUS.ANCHORED <- IRanges::subsetByOverlaps(x = sw.SECONDARY.MULT.PLUS.RED, ranges = BM.CORES, minoverlap = MIN.OVERLAP)
      sw.SECONDARY.MULT.MINUS.ANCHORED <- IRanges::subsetByOverlaps(x = sw.SECONDARY.MULT.MINUS.RED, ranges = BM.CORES, minoverlap = MIN.OVERLAP)

      sw.SECONDARY.MULT.RED.ANCHORED <- GenomicRanges::sort.GenomicRanges(GenomicRanges::reduce(c(sw.SECONDARY.MULT.PLUS.ANCHORED,sw.SECONDARY.MULT.MINUS.ANCHORED)))
      sw.SECONDARY.MULT.RED.ANCHORED <- GenomicRanges::sort.GenomicRanges(GenomicRanges::reduce(c(sw.SECONDARY.MULT.RED.ANCHORED,BM.CORES)))

      CLUSTERS.GR <- c(BM.CORES, sw.SECONDARY.MULT.RED.ANCHORED)
    } else {
      CLUSTERS.GR <- BM.CORES
    }

    CLUSTERS <- GenomicRanges::reduce(CLUSTERS.GR)
    ## Removing gaps between CLUSTERS
    if (THRESHOLD.CLUSTERS.GAP>0){
      GR.GAPS <- GenomicRanges::gaps(CLUSTERS)
      GR.GAPS <- GenomicRanges::sort.GenomicRanges(GR.GAPS)
      GR.GAPS.WIDER.THAN.MIN <- GR.GAPS[GenomicRanges::width(GR.GAPS) >= THRESHOLD.CLUSTERS.GAP]
      CLUSTERS <- GenomicRanges::gaps(GR.GAPS.WIDER.THAN.MIN) #gaps of gaps
      CLUSTERS <- GenomicRanges::sort.GenomicRanges(CLUSTERS)
    }
    ##
    if (VERBOSITY>0) message("\nBuilding ... STEP 2 ... Annotating & filtering\n")
    ##

    ##
    ## Annotate the cores

    if (VERBOSITY>1) message("\tAnnotating cores and clusters")
    outputList[[uniqueandprimary]]<-BM.CORES
    outputList[[allalignments]]<-CLUSTERS
    outputList<-PICBannotate(outputList, IN.ALIGNMENTS, , REFERENCE.GENOME = REFERENCE.GENOME,
                             SEQ.LEVELS.STYLE = SEQ.LEVELS.STYLE, PROVIDE.NON.NORMALIZED = TRUE, COMPUTE.1U.10A.FRACTIONS = COMPUTE.1U.10A.FRACTIONS, LIBRARY.SIZE = LIBRARY.SIZE)

    ##
    ## REPORT

    if (VERBOSITY>1) message("\n\tClusters stats: ")
    used_uniq <- sum(GenomicRanges::mcols(outputList[[allalignments]])[["uniq_reads"]])
    used_uniq_piRNA <- round(((used_uniq/length(IN.ALIGNMENTS$unique))*100),3)
    used_mult <- sum(GenomicRanges::mcols(outputList[[allalignments]])[["multimapping_reads_primary_alignments"]])
    used_mult_primary_piRNA <- round(((used_mult/length(WGRMP))*100),3)
    if (VERBOSITY>1) {
      message("\t\tAccomodated UNIQUE MAPPERS: ",used_uniq," (",used_uniq_piRNA," %)")
      message("\t\tAccomodated PRIMARY MULTI MAPPERS: ",used_mult," (",used_mult_primary_piRNA," %)")
    }

    if (! PROVIDE.NON.NORMALIZED){ #removing stats hard to understand
      for ( t in c(uniqueonly, uniqueandprimary, allalignments)){
        GenomicRanges::mcols(outputList[[t]])[["uniq_reads"]] <- NULL
        GenomicRanges::mcols(outputList[[t]])[["uniq_sequences"]]  <- NULL
        GenomicRanges::mcols(outputList[[t]])[["multimapping_reads_primary_alignments"]] <- NULL
        GenomicRanges::mcols(outputList[[t]])[["all_reads_primary_alignments"]] <- NULL
        GenomicRanges::mcols(outputList[[t]])[["width_covered_by_unique_alignments"]] <- NULL
      }
    }
    if (VERBOSITY>0) message("\nDone!")
    return(outputList)

  }
