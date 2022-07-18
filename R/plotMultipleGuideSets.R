#' @title Plotting a GuideSet and other genomic annotation
#' 
#' @description Function to plot guide targets stored in a
#'     \linkS4class{GuideSet} object in a gene browser view supported by
#'     \code{Gviz}. Target gene isoforms and other genomic annotation,
#'     along with the target chromosome ideogram and sequence, may also be
#'     added, permitting a comprehensive visualization of the genomic context
#'     targeted by spacers in the \linkS4Class{GuideSet}.
#' 
#' @param x named list of GuideSets
#' @param geneModel tbd
#' @param targetGene gene symbol or id
#' @param annotations named list of annotations to plot (GRanges)
#' @param from,to see ?plotTracks
#' @param extend.left,extend.right see ?plotTracks
#' @param includeIdeogram boolean
#' @param bands see ?IdeogramTrack
#' @param bsgenome used in SequenceTrack
#' @param onTargetScores list of on-target scores
#' @param includeSNPTrack tbd
#' @param displayPars tbd (possible to implement?)
#' 
#' @return A Gviz plot... see ?Gviz::plotTracks
#' 
#' @author Luke Hoberecht, Jean-Philippe Fortin
#' 
#' @seealso \code{\link{plotMultipleGuideSets}}
#' 
#' @examples
#' \dontrun{
#' ## WIP =====================================================================
#' library(crisprDesign)
#' data(guideSetExample)
#' library(crisprDesignData)
#' data(txdb_human)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' gs <- guideSetExample[1:10]
#' targetGene <- "IQSEC3"
#' geneModel <- txdb_human
#' plotGuideSet(gs,
#'              geneModel=geneModel,
#'              targetGene=targetGene)
#'
#' kras <- queryTxObject(txdb_human, "cds", "gene_symbol", "KRAS")
#' kras <- kras[1]
#' gs <- findSpacers(kras, crisprNuclease=SpCas9, bsgenome=BSgenome.Hsapiens.UCSC.hg38)
#' ideogram <- readRDS("inst/shiny/data/ideogram_hg38.rds")
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
                                  includeIdeogram=TRUE,
                                  bands=NULL,
                                  bsgenome=NULL,
                                  onTargetScores=NULL,
                                  includeSNPTrack=TRUE,
                                  displayPars=NULL
){
    ## check inputs
    guideSets <- .validateGuideSetList(x)
    chr <- .validateSingleChrList(guideSets)
    geneModel <- .validateGRangesList(geneModel)
    stopifnot("targetGene must be a character vector or NULL" = {
        is.vector(targetGene, mode="character") || is.null(targetGene)
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
        from <- plotWindowMin - plotWindowSize
    }
    if (is.null(to)){
        to <- plotWindowMax + plotWindowSize
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
    
    tracks <- c(list(ideogramTrack, genomeAxisTrack), geneTrack, annotationTrack, snpTrack, guideTrack)
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






