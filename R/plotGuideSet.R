#' @title Plotting a GuideSet and other genomic annotation
#' 
#' @description Function to plot guide targets stored in a
#'     \linkS4class{GuideSet} object in a gene browser view supported by
#'     \code{Gviz}. Target gene isoforms and other genomic annotation,
#'     along with the target chromosome ideogram and sequence, may also be
#'     added, permitting a comprehensive visualization of the genomic context
#'     targeted by spacers in the \linkS4Class{GuideSet}.
#' 
#' @param x tbd
#' @param geneModel tbd
#' @param targetGene gene symbol or Ensembl ID
#' @param annotations  named list of annotations to plot (GRanges)
#' @param from,to  see ?plotTracks
#' @param extend.left,extend.right  see ?plotTracks
#' @param bands see ?IdeogramTrack
#' @param guideStacking  how to stack all guides on a single track (squish, dense, or hide) or NA to have each guide occupy a separate track
#' @param bsgenome used in SequenceTrack
#' @param pamSiteOnly whether to plot only the PAM site for spacers or plot the spacer and PAM sequence (default)
#' @param showGuideLabels tbd
#' @param onTargetScore  color coding guides...need to add legend
#' @param includeSNPTrack tbd
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
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom Gviz plotTracks
plotGuideSet <- function(x,
                         geneModel=NULL,
                         targetGene=NULL,
                         annotations=list(),
                         from=NULL,
                         to=NULL,
                         extend.left=0,
                         extend.right=0,
                         includeIdeogram=TRUE,
                         bands=NULL,
                         guideStacking="squish",
                         bsgenome=NULL,
                         pamSiteOnly=FALSE,
                         showGuideLabels=TRUE,
                         onTargetScore=NULL,
                         includeSNPTrack=TRUE
){
    ## check inputs
    guideSet <- .validateGuideSet(x)
    chr <- .validateSingleChr(guideSet)
    geneModel <- .validateGRangesList(geneModel)
    stopifnot("targetGene must be a character vector or NULL" = {
        is.vector(targetGene, mode="character") || is.null(targetGene)
    })
    stopifnot("guideStacking must be NA or one of 'squish', 'dense', 'hide', or 'none'" = {
        length(guideStacking) == 1 &&
            (
                is.na(guideStacking) ||
                    guideStacking %in% c("squish", "dense", "hide", "none")
            )
    })
    stopifnot("pamSiteOnly must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(pamSiteOnly)
    })
    
    ## set tracks
    ideogramTrack <- .getIdeogramTrack(
        includeIdeogram=includeIdeogram,
        chr=chr,
        genome=unique(GenomeInfoDb::genome(guideSet)),
        bands=bands
    )
    genomeAxisTrack <- .getGenomeAxisTrack()
    guideTrack <- .getGuideTrack(
        guideSet=guideSet,
        guideStacking=guideStacking,
        pamSiteOnly=pamSiteOnly,
        showGuideLabels=showGuideLabels,
        onTargetScore=onTargetScore
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
        snpTrack <- .getSnpTrack(
            guideSet=guideSet,
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



#' @importFrom Gviz IdeogramTrack
.getIdeogramTrack <- function(includeIdeogram,
                              chr,
                              genome,
                              bands
){
    if (!includeIdeogram || (is.null(genome) && is.null(bands))){
        return(NULL)
    } else {
        Gviz::IdeogramTrack(
            chromosome=chr,
            genome=genome,
            bands=bands,
            fontcolor="black",
            cex=1
        )
    }
}



#' @importFrom Gviz GenomeAxisTrack
.getGenomeAxisTrack <- function(
){
    Gviz::GenomeAxisTrack(
        add35=TRUE,
        add53=TRUE,
        cex=1,
        fontcolor="black"
    )
}









#' @importFrom crisprDesign queryTxObject
#' @importFrom BiocGenerics start end strand
#' @importFrom IRanges restrict
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames seqinfo seqlengths
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicFeatures makeTxDb
#' @importFrom Gviz GeneRegionTrack transcript strand
.getGeneTrack <- function(geneModel,
                          targetGene,
                          from,
                          to
                        
){
    if (is.null(geneModel) || is.null(targetGene)){
        return()
    }
    
    transcripts <- crisprDesign::queryTxObject(
        txObject=geneModel,
        featureType="transcripts",
        queryColumn="gene_symbol",
        queryValue=targetGene
    )
    if (length(transcripts) == 0){
        warning("targetGene not found in geneModel, ignoring.",
                immediate.=TRUE)
        return()
    }
    ## get tx_id whose start/end exceed limits (for each limit)
    txLeft <- transcripts$tx_id[BiocGenerics::start(transcripts) < from &
                                    BiocGenerics::end(transcripts) > from]
    txRight <- transcripts$tx_id[BiocGenerics::end(transcripts) > to &
                                     BiocGenerics::start(transcripts) < to]
    extendDistance <- 0.1 * (to - from)
    transcripts <- IRanges::restrict(transcripts,
                                     start=from,
                                     end=to)
    splicings <- crisprDesign::queryTxObject(
        txObject=geneModel,
        featureType="exons",
        queryColumn="gene_symbol",
        queryValue=targetGene
    )
    splicings <- IRanges::restrict(splicings,
                                   start=from,
                                   end=to)
    splicings$exon_rank <- splicings$exon_rank + 1 # for upstream shift, if necessary
    
    ## for tx_id above that exceed limit, append "exon" beyond limit
    exonLeft <- lapply(txLeft, function(x){
        strand <- unique(BiocGenerics::strand(splicings))
        currentRanks <- splicings$exon_rank[splicings$tx_id == x]
        if (length(currentRanks) == 0){
            exon_rank <- ifelse(strand == "+", 1, 2)
        } else if (strand == "+"){
            exon_rank <- min(currentRanks) - 1
        } else {
            exon_rank <- max(currentRanks) + 1
        }
        GenomicRanges::GRanges(
            seqnames=unique(GenomeInfoDb::seqnames(splicings)),
            ranges=IRanges::IRanges(start=from-extendDistance, width=1),
            strand=strand,
            tx_id=x,
            gene_id=unique(splicings$gene_id),
            protein_id=NA,
            tx_type=NA,
            gene_symbol=NA,
            exon_id=0,
            exon_rank=exon_rank,
            cds_start=NA,
            cds_end=NA,
            tx_start=NA,
            tx_end=NA,
            cds_len=NA,
            exon_start=NA,
            exon_end=NA
        )
    })
    exonLeft <- Reduce(c, exonLeft)
    
    exonRight <- lapply(txRight, function(x){
        strand <- unique(BiocGenerics::strand(splicings))
        currentRanks <- splicings$exon_rank[splicings$tx_id == x]
        if (length(currentRanks) == 0){
            exon_rank <- ifelse(strand == "+", 2, 1)
        } else if (strand == "+"){
            exon_rank <- max(currentRanks) + 1
        } else {
            exon_rank <- min(currentRanks) - 1
        }
        GenomicRanges::GRanges(
            seqnames=unique(GenomeInfoDb::seqnames(splicings)),
            ranges=IRanges::IRanges(start=to+extendDistance, width=1),
            strand=strand,
            tx_id=x,
            gene_id=unique(splicings$gene_id),
            protein_id=NA,
            tx_type=NA,
            gene_symbol=NA,
            exon_id=999999999,
            exon_rank=exon_rank,
            cds_start=NA,
            cds_end=NA,
            tx_start=NA,
            tx_end=NA,
            cds_len=NA,
            exon_start=NA,
            exon_end=NA
        )
    })
    exonRight <- Reduce(c, exonRight)
    
    splicings <- c(splicings, exonLeft, exonRight)
    
    ## adjust start/end of cds and exons
    for (i in seq_along(splicings)){
        tx_id <- splicings$tx_id[i]
        tx_index <- transcripts$tx_id == tx_id
        splicings$cds_start[i] <- max(BiocGenerics::start(splicings[i]), splicings$cds_start[i])
        splicings$cds_start[i] <- min(splicings$cds_start[i],
                                      BiocGenerics::end(transcripts[tx_index]))
        splicings$cds_end[i] <- min(BiocGenerics::end(splicings[i]), splicings$cds_end[i])
        splicings$cds_end[i] <- max(splicings$cds_end[i],
                                      BiocGenerics::start(transcripts[tx_index]))
    }
    
    chrominfo <- GenomeInfoDb::seqinfo(geneModel[["transcripts"]])
    chrominfo <- data.frame(
        chrom=GenomeInfoDb::seqnames(chrominfo),
        length=GenomeInfoDb::seqlengths(chrominfo)
    )
    
    .id2Integer <- function(x){
        as.integer(gsub('[A-Z]', '', x))
    }
    
    transcripts <- data.frame(
        tx_id=.id2Integer(transcripts$tx_id),
        tx_chrom=as.character(GenomeInfoDb::seqnames(transcripts)),
        tx_strand=as.character(BiocGenerics::strand(transcripts)),
        tx_start=as.integer(BiocGenerics::start(transcripts)),
        tx_end=as.integer(BiocGenerics::end(transcripts)),
        tx_name=transcripts$tx_id,
        gene_id=transcripts$gene_id)
    splicings <- data.frame(
        tx_id=.id2Integer(splicings$tx_id),
        exon_rank=splicings$exon_rank,
        exon_id=.id2Integer(splicings$exon_id),
        exon_start=as.integer(BiocGenerics::start(splicings)),
        exon_end=as.integer(BiocGenerics::end(splicings)),
        cds_start=splicings$cds_start,
        cds_end=splicings$cds_end)
    
    txdb <- suppressWarnings(
        GenomicFeatures::makeTxDb(
            transcripts=transcripts,
            splicings=splicings,
            chrominfo=chrominfo
        )
    )
    
    txTrack <- Gviz::GeneRegionTrack(
        range=txdb,
        name=targetGene,
        showID=TRUE,
        # fill="darkorange",
        geneSymbol=TRUE,
        transcriptAnnotation="transcript",
        fontcolor.group="black",
        col.line="black"
    )
    
    
    ## separate function to add direction to labels
    txTrack <- .addStrandToLabel(
        track=txTrack,
        label=paste0(targetGene, " (", Gviz::transcript(txTrack), ")"),
        strand=Gviz::strand(txTrack)
    )
    return(txTrack)
}





#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom Gviz AnnotationTrack
.getAnnotationTrack <- function(annotations,
                                chr,
                                from,
                                to
){
    windowRange <- GenomicRanges::GRanges(
        seqnames=chr,
        ranges=IRanges::IRanges(start=from, end=to)
    )
    annTracks <- lapply(seq_along(annotations), function(x){
        overlaps <- suppressWarnings(GenomicRanges::findOverlaps(
            annotations[[x]],
            windowRange,
            ignore.strand=TRUE
        ))
        hits <- S4Vectors::queryHits(overlaps)
        ann <- annotations[[x]][hits]
        Gviz::AnnotationTrack(
            ann,
            shape="arrow",
            fill="#F5CBD0",
            showId=FALSE,
            name=names(annotations)[x]
        )
    })
    return(annTracks)
}





#' @importFrom crisprDesign snps
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom Gviz AnnotationTrack
.getSnpTrack <- function(guideSet,
                         chr
){
    snps <- crisprDesign::snps(guideSet)
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
.getGuideTrack <- function(guideSet,
                           guideStacking,
                           pamSiteOnly,
                           showGuideLabels,
                           onTargetScore
){
    # if (guideStacking == "none"){
    #     return()
    # }
    if (pamSiteOnly){
        guideRanges <- guideSet
        if (is.na(guideStacking)){
            guideRanges <- split(guideRanges, f=names(guideRanges))
        } else {
            guideRanges <- list(guideRanges)
        }
    } else {
        guideRanges <- .extractGuideRanges(
            guideSet=guideSet,
            guideStacking=guideStacking)
    }
    
    if (is.na(guideStacking)){
        guideStacking <- "squish"
        name <- ""
    } else {
        name <- "gRNAs"
    }
    
    colors <- .getScoreColors(
        guideSet=guideSet,
        names=names,
        onTargetScore=onTargetScore
    )
    
    guideTracks <- lapply(guideRanges, function(x){
        if (methods::is(x, "TxDb")){
            names <- GenomicFeatures::transcripts(x)$tx_name
        } else {
            names <- names(x)
        }
        track <- Gviz::GeneRegionTrack(
            range=x,
            name=name,
            symbol=names,
            showID=TRUE,
            geneSymbol=TRUE,
            stacking=guideStacking,
            fontcolor.group="black",
            col.line="black"
        )
        if (methods::is(x, "GuideSet")){
            Gviz::group(track) <- names
        }
        groups <- Gviz::group(track)
        Gviz::displayPars(track)[["fill"]] <- colors[groups]
        if (!methods::is(x, "TxDb")){
            Gviz::transcript(track) <- names
        }
        if (showGuideLabels){
            Gviz::displayPars(track) <- list(transcriptAnnotation="transcript")
            track <- .addStrandToLabel(
                track=track,
                label=Gviz::transcript(track),
                strand=Gviz::strand(track)
            )
        }
        track
    })
    return(guideTracks)
}




#' @importFrom crisprBase getProtospacerRanges getPamRanges
#' @importFrom crisprDesign crisprNuclease
#' @importFrom GenomicRanges punion
#' @importFrom methods as
#' @importFrom S4Vectors mcols<-
#' @importFrom BiocGenerics start end strand
#' @importFrom GenomeInfoDb seqnames seqinfo seqlengths
#' @importFrom GenomicFeatures makeTxDb
.extractGuideRanges <- function(guideSet,
                                guideStacking
){
    protospacerRanges <- crisprBase::getProtospacerRanges(
        gr=guideSet,
        nuclease=crisprDesign::crisprNuclease(guideSet),
    )
    pamRanges <- crisprBase::getPamRanges(
        gr=guideSet,
        nuclease=crisprDesign::crisprNuclease(guideSet)
    )
    transcripts <- GenomicRanges::punion(
        methods::as(protospacerRanges, "GRanges"),
        methods::as(pamRanges, "GRanges")
    )
    
    S4Vectors::mcols(transcripts)[["tx_id"]] <- seq_along(transcripts)
    S4Vectors::mcols(transcripts)[["exon_rank"]] <- 1
    S4Vectors::mcols(transcripts)[["exon_id"]] <- as.integer(seq_along(transcripts))
    S4Vectors::mcols(transcripts)[["cds_start"]] <- BiocGenerics::start(protospacerRanges)
    S4Vectors::mcols(transcripts)[["cds_end"]] <- BiocGenerics::end(protospacerRanges)
    
    tx <- data.frame(
        tx_id=transcripts$tx_id,
        tx_chrom=as.character(GenomeInfoDb::seqnames(transcripts)),
        tx_strand=as.character(BiocGenerics::strand(transcripts)),
        tx_start=as.integer(BiocGenerics::start(transcripts)),
        tx_end=as.integer(BiocGenerics::end(transcripts)),
        tx_name=names(guideSet),
        gene_id="gRNAs")
    splicings <- data.frame(
        tx_id=transcripts$tx_id,
        tx_strand=as.character(BiocGenerics::strand(transcripts)),
        exon_rank=transcripts$exon_rank,
        exon_id=transcripts$exon_id,
        exon_start=as.integer(BiocGenerics::start(transcripts)),
        exon_end=as.integer(BiocGenerics::end(transcripts)),
        cds_start=transcripts$cds_start,
        cds_end=transcripts$cds_end)
    chrominfo <- GenomeInfoDb::seqinfo(guideSet)
    chrominfo <- data.frame(
        chrom=GenomeInfoDb::seqnames(chrominfo),
        length=GenomeInfoDb::seqlengths(chrominfo)
    )
    
    if (is.na(guideStacking)){
        txdb <- lapply(seq_len(nrow(tx)), function(x){
            suppressWarnings(
                GenomicFeatures::makeTxDb(
                    transcripts=tx[x, , drop=FALSE],
                    splicings=splicings[x, , drop=FALSE],
                    chrominfo=chrominfo
                )
            )
        })
    } else {
        txdb <- list(
            suppressWarnings(
                GenomicFeatures::makeTxDb(
                    transcripts=tx,
                    splicings=splicings,
                    chrominfo=chrominfo
                )
            )
        )
    }
    return(txdb)
}











#' @importFrom methods is
#' @importFrom Gviz SequenceTrack
.addSequenceTrack <- function(tracks,
                              chr,
                              bsgenome
){
    stopifnot("bsgenome must be a BSgenome object." = {
        methods::is(bsgenome, "BSgenome")
    })
    seqTrack <- Gviz::SequenceTrack(
        sequence=bsgenome,
        chromosome=chr,
        add53=TRUE,
        noLetters=FALSE
    )
    tracks <- c(tracks, seqTrack)
    return(tracks)
}






## other ======================================================================



#' @importFrom Gviz transcript
.addStrandToLabel <- function(track,
                              label,
                              strand
){
    forwardStrand <- strand == "+"
    label[forwardStrand] <- paste(label[forwardStrand], "->", recycle0=TRUE)
    label[!forwardStrand] <- paste("<-", label[!forwardStrand], recycle0=TRUE)
    Gviz::transcript(track) <- label
    return(track)
}










#' @importFrom S4Vectors mcols
#' @importFrom grDevices colorRampPalette
.getScoreColors <- function(guideSet,
                            names,
                            onTargetScore,
                            score0="#D2D2D2",
                            score1="#000080"
){
    onTargetScore <- .validateOnTargetScore(
        onTargetScore,
        guideSet=guideSet
    )
    if (is.null(onTargetScore)){
        return("lightblue") # default color
    }
    scores <- S4Vectors::mcols(guideSet)[[onTargetScore]]
    colorFunction <- grDevices::colorRampPalette(c(score0, score1))
    colors <- colorFunction(100)
    intervals <- cut(scores, breaks=seq(0, 1, by=0.01))
    intervals <- as.numeric(intervals)
    scoreColors <- colors[intervals]
    names(scoreColors) <- names(guideSet)
    return(scoreColors)
}
