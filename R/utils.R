


#' @importFrom methods is
.validateGuideSet <- function(obj,
                              eMessage=NULL
){
    if (is.null(eMessage)){
        eMessage <- "guideSet argument must be a GuideSet object."
    }
    isGuideSet <- methods::is(obj, "GuideSet")
    if (!isGuideSet){
        stop(eMessage)
    }
    return(obj)
}






.validateGuideSetList <- function(obj,
                                  eMessage=NULL
){
    if (is.null(eMessage)){
        eMessage <- "guideSets argument must be a list of GuideSets."
    }
    if (!is.list(obj)){
        stop(eMessage)
    }
    lapply(obj, .validateGuideSet, eMessage=eMessage)
    return(obj)
}








#' @importFrom GenomeInfoDb seqnames
.validateSingleChr <- function(obj
){
    chr <- GenomeInfoDb::seqnames(obj)
    chr <- as.character(unique(chr))
    stopifnot("Only GuideSets targeting a single chromosome may be plotted." = {
        length(chr) == 1
    })
    return(chr)
}



#' @importFrom GenomeInfoDb seqnames
.validateSingleChrList <- function(obj
){
    eMessage <- "All GuideSets must target the same chromosome."
    if (!is.list(obj)){
        stop(eMessage)
    }
    chrs <- vapply(obj, function(x){
        chr <- GenomeInfoDb::seqnames(x)
        as.character(unique(chr))
    }, FUN.VALUE=character(1))
    chr <- unique(chrs)
    if (length(chr) != 1){
        stop(eMessage)
    }
    return(chr)
}






#' @importFrom methods is
#' @importFrom crisprDesign TxDb2GRangesList
.validateGRangesList <- function(obj
){
    if (methods::is(obj, "TxDb")){
        obj <- crisprDesign::TxDb2GRangesList(obj)
        return(obj)
    }
    stopifnot("Object must be a TxDb or GRangesList." = {
        methods::is(obj, "GRangesList")
    })
    fields <- c("transcripts",
                "exons",
                "cds",
                "fiveUTRs",
                "threeUTRs",
                "introns",
                "tss")
    if (!all(fields %in% names(obj))){
        stop("GRangesList is missing required genomic regions. ",
             "Use crisprDesign::TxDb2GRangesList() to get properly ",
             "formatted GRangesList.")
    }
    return(obj)
}






.validateOnTargetScore <- function(onTargetScore,
                                   guideSet
){
    if (!is.null(onTargetScore)){
        stopifnot("onTargetScore must be a character string" = {
            length(onTargetScore) == 1 &&
                is.vector(onTargetScore, mode="character")
        })
        mcolnames <- names(S4Vectors::mcols(guideSet))
        if (!onTargetScore %in% mcolnames){
            warning("Could not find '", onTargetScore,
                    "' column in guideSet, ignoring.",
                    immediate.=TRUE)
            onTargetScore <- NULL
        }
    }
    return(onTargetScore)
}

