#' Runs PICBbuild multiple times with provided parameters and returns optimization data frame.
#'
#' @param IN.ALIGNMENTS list of alignments from PICBload
#' @param REFERENCE.GENOME name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6"
#' @param VERBOSITY verbosity level 0/1/2/3. 2 by default.
#' @param ... rest of the parameters used by PICBbuild and provided as iteratable vectors
#'
#' @author Aleksandr Friman
#' @return Optimization values dataframe
#' @export
#'
#' @examples PICBoptimize(IN.ALIGNMENTS = myAlignmentsFromPIBCload, REFERENCE.GENOME="BSgenome.Dmelanogaster.UCSC.dm6", MIN.UNIQUE.ALIGNMENTS.PER.WINDOW=c(1,2,3,4,5))
PICBoptimize <- function(
    IN.ALIGNMENTS,
    REFERENCE.GENOME,
    VERBOSITY=2,
    SEQ.LEVELS.STYLE="UCSC",
    ...
){
  numberOfallReadsExplained<-function(gr, alignments){
    foundnames<-c()
    if (typeof(alignments)=="list"){
      for (clmn in c("unique", "multi.primary", "multi.secondary")){
        usedAlignments<-IRanges::subsetByOverlaps(alignments[[clmn]], gr)
        foundnames<-c(foundnames, names(usedAlignments))
      }
    }else{
      usedAlignments<-IRanges::subsetByOverlaps(alignments, gr)
      foundnames<-c(foundnames, names(usedAlignments))
    }
    numbr=length(unique(foundnames))
    return(numbr)
  }
  
  inputArgs=as.list(match.call(expand.dots=TRUE))
  PICBbuildArgs=formals(PICBbuild)
  ArgsToIterate=NULL
  for (argName in names(PICBbuildArgs)){
    #checking if the argument need to be iterated over
    if (is.numeric(PICBbuildArgs[[argName]]) ||
        (is.language(PICBbuildArgs[[argName]]) && PICBbuildArgs[[argName]] != "") ){
      # iteration is possible
      if (argName %in% names(inputArgs)){
        #need to iterate
        if (length(inputArgs[[argName]]) > 1){
          if (VERBOSITY >0 ){
            message(paste("Adding", argName, "to iteration targets"))
          }
          if (is.null(ArgsToIterate)){
            ArgsToIterate=argName
          }else{
            ArgsToIterate=c(ArgsToIterate, argName)
          }
        }
      }
    }
  }
  #checking the inputs
  if (VERBOSITY >0 ){
    message("Checking the inputs")
  }
  if (is.null(ArgsToIterate)){
    stop("Provide arguments to iterave over. See example.")
  }
  totalReads=length(IN.ALIGNMENTS$unique) + length(IN.ALIGNMENTS$multi.primary)
  if (length(unique(c(names(IN.ALIGNMENTS$unique),
                      names(IN.ALIGNMENTS$multi.primary)))
  ) != totalReads){
    warning("Total number of primary aligments is not equal to total number of read names. Numbers of explained reads may be wrong.")   
  }
  #building the parameter combinations
  if (VERBOSITY >0 ){
    message("Building the parameter combinations")
  }
  listOfArgLists=list()
  for (argName in ArgsToIterate){
    firstArg= length(listOfArgLists)==0
    argVals=eval(inputArgs[[argName]])
    if (firstArg){ #first arg
      for (curVal in argVals){
        tmpList=list()
        tmpList[[argName]]=curVal
        listOfArgLists[[length(listOfArgLists)+1]]=tmpList
      }
    }else{
      #add the argName to each combination
      newlistOfArgLists=list()
      for (combList in listOfArgLists){
        for (curVal in argVals){
          combList[[argName]]=curVal
          newlistOfArgLists[[length(newlistOfArgLists)+1]]=combList
        }
      }
      listOfArgLists=newlistOfArgLists
    }
    
  }
  IN.ALIGNMENTS.ALL=c(IN.ALIGNMENTS$unique, IN.ALIGNMENTS$multi.primary,
                      IN.ALIGNMENTS$multi.secondary)
  #preparing output dataframe
  outDF=data.frame(combinationNumber=1:length(listOfArgLists))
  for (argName in names(listOfArgLists[[1]])){
    #adding columns we iterate over
    outDF[[argName]]=NA
  }
  for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
    outDF[[paste0("number.of.", locustype)]] = NA
  }
  for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
    outDF[[paste0(locustype, ".total.width")]] = NA
  }
  for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
    outDF[[paste0("reads.explained.by.", locustype)]] = NA
  }
  #running PICBbuild for every possible combation
  i=1
  for (combList in listOfArgLists){
    if (VERBOSITY >0 ){
      message(paste("Iteration", i, "out of", length(listOfArgLists)))
    }
    for (argName in names(combList)){
      #adding values we iterate over
      outDF[[argName]][i]=combList[[argName]]
    }
    listToRun=combList
    listToRun[["IN.ALIGNMENTS"]]=IN.ALIGNMENTS
    listToRun[["REFERENCE.GENOME"]]=REFERENCE.GENOME
    listToRun[["VERBOSITY"]]=VERBOSITY
    if (VERBOSITY >1 ){
      message("Running PICBbuild")
    }
    tmrLoci=do.call(what = PICBbuild, args = listToRun)
    if (VERBOSITY >1 ){
      message("PICBbuild done")
    }
    for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
      outDF[[paste0("number.of.", locustype)]][i]= length(tmrLoci[[locustype]])
      outDF[[paste0(locustype, ".total.width")]][i]= as.numeric(sum(GenomicRanges::width(tmrLoci[[locustype]])))
      outDF[[paste0("reads.explained.by.", locustype)]][i]= as.numeric(numberOfallReadsExplained(tmrLoci[[locustype]], IN.ALIGNMENTS.ALL))
      
    }
    i=i+1
  }
  outDF$combinationNumber=NULL
  #adding normaized values
  if (VERBOSITY >0 ){
    message("Values normalization")
  }

  for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
    outDF[[paste0("fraction.of.library.explained.by.", locustype)]]= outDF[[paste0("reads.explained.by.", locustype)]]/totalReads
  }
  for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
    outDF[[paste0(locustype,".mean.RPKM")]]= outDF[[paste0("reads.explained.by.", locustype)]]*1e9/(outDF[[paste0(locustype, ".total.width")]]*totalReads)
    outDF[[paste0(locustype,".mean.RPKM")]][ is.na(outDF[[paste0(locustype,".mean.RPKM")]])] = 0
  }
  if (typeof(REFERENCE.GENOME)=="character"){
    SI<-PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE)
  }else{
    SI<-REFERENCE.GENOME
  }
  GenomeSize=sum(GenomeInfoDb::seqlengths(SI))
  for (locustype in c(uniqueonly,uniqueandprimary,allalignments )){
    outDF[[paste0(locustype, ".fraction.of.genome.space")]]=outDF[[paste0(locustype, ".total.width")]]/(2*GenomeSize)
  }
  return(outDF)
}