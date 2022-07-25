#' @title Plotting a GuideSet and other genomic annotation
#' 
#' @description Function to plot guide targets stored in multiple
#'     \linkS4class{GuideSet} objects in a gene browser view supported by
#'     \code{Gviz}. Target gene isoforms and other genomic annotation,
#'     along with the target chromosome ideogram and sequence, may also be
#'     added, permitting a comprehensive visualization of the genomic context
#'     targeted by spacers in \linkS4class{GuideSet} objects.
#' 
#' @param x A named (optional) list of \linkS4class{GuideSet} objects. List
#'     names are displayed as track titles.
#' @param geneModel A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \link[crisprDesign]{TxDb2GRangesList}.
#' @param targetGene A character string giving the gene symbol or Ensembl ID of
#'     the gene to plot.
#' @param annotations A named (optional) list of genomic annotations as
#'     \linkS4class{GRanges} to plot. Provided names are displayed as track
#'     titles.
#' @param from,to Numeric value giving the genomic coordinate range to plot;
#'     see \link[Gviz]{plotTracks}.
#' @param extend.left,extend.right  Numeric value giving the length in bases to
#'     extend the plotting range; see \link[Gviz]{plotTracks}.
#' @param margin A numeric value that sets the margin of the plotting range
#'     with respect to the range of \code{x}. The value is a ratio of the width
#'     of \code{x} considered as a single range, with a value of \code{1}
#'     (default) extending the plotting window by the same distance as that
#'     between \code{min(start(x))} and \code{max(end(x))}. Only works for the
#'     respective side of the plot if either \code{from} or \code{to} are
#'     \code{NULL}.
#' @param includeIdeogram Logical; whether to include an
#'     \link[Gviz]{IdeogramTrack} in the plot.
#' @param bands An optional \code{data.frame} of cytoband information for the 
#'     target genome; see \link[Gviz]{IdeogramTrack}.
#' @param bsgenome A \linkS4class{BSgenome} object; used to generate
#'     \link[Gviz]{SequenceTrack} and GC content \link[Gviz]{DataTrack}.
#' @param onTargetScores Optional list of column names for each element in
#'     \code{x} of on-target scores. Applies a color scheme to the guide track
#'     based on on-target scores, with light gray corresponding to 0 and dark
#'     blue corresponding to 1.
#' @param includeSNPTrack Logical; whether to include an
#'     \link[Gviz]{AnnotationTrack} for SNPs if such annotation exists in
#'     \code{x}.
#' @param gcWindow If not \code{NULL}, a numeric value specifying the distance
#'     from a given base for which to establish a window for calculating GC
#'     content at that base. These values are then added to the plot in a
#'     \link[Gviz]{DataTrack}.
#' 
#' @return A Gviz plot; see \code{\link[Gviz]{plotTracks}}.
#' 
#' @author Luke Hoberecht, Jean-Philippe Fortin
#' 
#' @seealso \code{\link{plotGuideSet}}
#' 
#' @examples
#' if (interactive()){
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' data(txdb_human, package="crisprDesignData")
#' data(cas9GuideSet, package="crisprViz")
#' data(cas12aGuideSet, package="crisprViz")
#' plotMultipleGuideSets(list(SpCas9=cas9GuideSet, AsCas12a=cas12aGuideSet),
#'                       geneModel=txdb_human,
#'                       targetGene="LTN1",
#'                       bsgenome=BSgenome.Hsapiens.UCSC.hg38,
#'                       margin=0.2,
#'                       gcWindow=10)
#' }
#' 
#' @export
plotMultipleGuideSets <- function(x,
                                  geneModel=NULL,
                                  targetGene=NULL,
                                  annotations=list(),
                                  from=NULL,
                                  to=NULL,
                                  extend.left=0,
                                  extend.right=0,
                                  margin=1,
                                  includeIdeogram=TRUE,
                                  bands=NULL,
                                  bsgenome=NULL,
                                  onTargetScores=NULL,
                                  includeSNPTrack=TRUE,
                                  gcWindow=NULL
){
    ## check inputs
    guideSets <- .validateGuideSetList(x)
    chr <- .validateSingleChrList(guideSets)
    geneModel <- .validateGRangesList(geneModel)
    stopifnot("targetGene must be a character vector or NULL" = {
        is.vector(targetGene, mode="character") || is.null(targetGene)
    })
    stopifnot("margin must be a non-negative numeric value" = {
        is.vector(margin, mode="numeric") &&
            length(margin) == 1 &&
            margin >= 0
    })
    
    ## set tracks
    genome <- lapply(guideSets, function(x){
        unique(GenomeInfoDb::genome(x))
    })
    genome <- unique(unlist(genome))
    ideogramTrack <- .getIdeogramTrack(
        includeIdeogram=includeIdeogram,
        chr=chr,
        genome=genome,
        bands=bands
    )
    genomeAxisTrack <- .getGenomeAxisTrack()
    guideTrack <- .getGuideTrackList(
        guideSets=guideSets,
        onTargetScores=onTargetScores
    )
    
    ## set plot range
    plotWindowMin <- lapply(guideTrack, min)
    plotWindowMin <- min(unlist(plotWindowMin))
    plotWindowMax <- lapply(guideTrack, max)
    plotWindowMax <- max(unlist(plotWindowMax))
    minimumWindowMargin <- 20
    plotWindowSize <- max(plotWindowMax - plotWindowMin,
                          minimumWindowMargin)
    if (is.null(from)){
        from <- plotWindowMin - margin * plotWindowSize
    }
    if (is.null(to)){
        to <- plotWindowMax + margin * plotWindowSize
    }
    
    ## have geneTrack adjust to plot window
    geneTrack <- .getGeneTrack(geneModel=geneModel,
                               targetGene=targetGene,
                               from=from - extend.left,
                               to=to + extend.right)
    
    annotationTrack <- .getAnnotationTrack(
        annotations=annotations,
        chr=chr,
        from=from,
        to=to
    )
    
    if (includeSNPTrack){
        snpTrack <- .getSnpTrackFromList(
            guideSets=guideSets,
            chr=chr
        )
    } else {
        snpTrack <- list()
    }
    
    if (!is.null(gcWindow)){
        stopifnot("gcWindow must be a positive integer" = {
            is.vector(gcWindow, "numeric") &&
                length(gcWindow) == 1 &&
                gcWindow == round(gcWindow)
        })
        gcTrack <- .getGcTrack(
            bsgenome=bsgenome,
            chr=chr,
            from=from - extend.left,
            to=to + extend.right,
            gcWindow=gcWindow
        )
    } else {
        gcTrack <- NULL
    }
    
    tracks <- c(list(ideogramTrack, genomeAxisTrack), geneTrack, annotationTrack, snpTrack, gcTrack, guideTrack)
    tracks <- tracks[vapply(tracks, function(x){
        !is.null(x)
    }, FUN.VALUE=logical(1))]
    if (!is.null(bsgenome)){
        tracks <- .addSequenceTrack(
            tracks=tracks,
            chr=chr,
            bsgenome=bsgenome
        )
    }
    
    Gviz::plotTracks(
        tracks,
        chromosome=chr,
        from=from,
        to=to,
        just.group="above",
        groupAnnotation="group",
        extend.left=extend.left,
        extend.right=extend.right,
        background.title="darkblue",
        background.panel="#F5F9FC"
    )
}







## tracks =====================================================================







#' @importFrom crisprDesign snps
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom Gviz AnnotationTrack
.getSnpTrackFromList <- function(guideSets,
                                 chr
){
    snps <- lapply(guideSets, crisprDesign::snps)
    snps <- Reduce(c, snps)
    if (is.null(snps) || nrow(snps) == 0){
        return(list())
    }
    snps <- snps[!duplicated(snps$rs), , drop=FALSE]
    
    snpRanges <- GenomicRanges::GRanges(
        seqnames=chr,
        ranges=IRanges::IRanges(start=snps$rs_site, width=snps$length)
    )
    names(snpRanges) <- snps$rs
    
    snpTrack <- Gviz::AnnotationTrack(
        snpRanges,
        chromosome=chr,
        shape="box",
        fill="#CDCDCD",
        group=names(snpRanges),
        name="SNP"
    )
    return(snpTrack)
}







#' @importFrom GenomicFeatures transcripts
#' @importFrom Gviz GeneRegionTrack displayPars<- transcript strand
.getGuideTrackList <- function(guideSets,
                               onTargetScores
){
    lapply(seq_along(guideSets), function(x){
        guideRanges <- guideSets[[x]]
        title <- names(guideSets)[x]
        onTargetScore <- onTargetScores[[x]]
        
        colors <- .getScoreColors(
            guideSet=guideRanges,
            onTargetScore=onTargetScore
        )
        
        track <- Gviz::GeneRegionTrack(
            range=guideRanges,
            name=title,
            stacking="dense",
            fontcolor.group="black",
            col.line="black"
        )
        Gviz::group(track) <- names(guideRanges)
        Gviz::displayPars(track)[["fill"]] <- colors[names(guideRanges)]
        track
    })
}






