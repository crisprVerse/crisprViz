
# library(crisprDesign)
# library(crisprDesignData)
# library(crisprViz)
# data(guideSetExample)
# data(txdb_human)
# gs <- guideSetExample[1:10]
# targetGene <- "IQSEC3"
# geneModel <- txdb_human
# plotGuideSets(list(gs),
#               ideogram=Ideogram_GRCh38,
#               geneModel=geneModel,
#               targetGene=targetGene,
#               from=66500,
#               to=67500,
#               genome="hg38")

# kras <- queryTxObject(txdb_human, "cds", "gene_symbol", "KRAS")
# kras <- kras[1]
# gs <- findSpacers(kras, crisprNuclease=SpCas9, bsgenome=BSgenome.Hsapiens.UCSC.hg38)
# ideogram <- readRDS("inst/shiny/data/ideogram_hg38.rds")



plotGuideSet <- function(guideSet,
                         bands=NULL, # see ?IdeogramTrack
                         bsgenome=NULL, # used in SequenceTrack
                         geneModel=NULL,
                         targetGene=NULL, # gene symbol or id
                         from=NULL, # see plotTracks
                         to=NULL, # see plotTracks
                         extend.left=0, # see plotTracks
                         extend.right=0, # see plotTracks
                         includeSequenceTrack=FALSE,
                         legendMetric=NULL#, # color coding guides...need to add legend
                         # displayOptions=NULL # named list of options that can be passed to subfunctions
){
    ## check inputs
    stopifnot("guideSet must be a GuideSet object" = {
        methods::is(guideSet, "GuideSet")
    })
    chr <- .validateSingleChr(guideSet)
    stopifnot("bsgenome must be a BSgenome object or NULL" = {
        methods::is(bsgenome, "BSgenome") || is.null(bsgenome)
    })
    stopifnot("targetGene must be a character vector or NULL" = {
        is.vector(targetGene, mode="character") || is.null(targetGene)
    })
    stopifnot(S4Vectors::isTRUEorFALSE(includeSequenceTrack))
    
    
    
    
    
    
    
    ## set tracks
    ideogramTrack <- .getIdeogramTrack(chr=chr,
                                       genome=unique(GenomeInfoDb::genome(guideSet)),
                                       bands=bands)
    genomeAxisTrack <- .getGenomeAxisTrack()
    guideTrack <- .getGuideTrack(guideSet)
    
    ## set plot range
    plotWindowMin <- min(guideTrack)
    plotWindowMax <- max(guideTrack)
    plotWindowSize <- plotWindowMax - plotWindowMin
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
    
    tracks <- c(list(ideogramTrack, genomeAxisTrack), geneTrack, guideTrack)
    tracks <- tracks[vapply(tracks, function(x){
        !is.null(x)
    }, FUN.VALUE=logical(1))]
    if (includeSequenceTrack){
        tracks <- c(tracks,
                    .getSequenceTrack(
                        chr=chr,
                        bsgenome=bsgenome
                    )
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
        background.panel="#F5F9FC")
}







## input validation ===========================================================



#' @importFrom GenomeInfoDb seqnames
.validateSingleChr <- function(guideSet
){
    chr <- GenomeInfoDb::seqnames(guideSet)
    chr <- as.character(unique(chr))
    stopifnot("Can only plot a single chromosome" = {
        length(chr) == 1
    })
    return(chr)
}
















## tracks =====================================================================



.getIdeogramTrack <- function(chr,
                              genome,
                              bands
){
    if (is.null(genome) && is.null(bands)){
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



.getGenomeAxisTrack <- function(
){
    Gviz::GenomeAxisTrack(
        add35=TRUE,
        add53=TRUE,
        cex=1,
        fontcolor="black"
    )
}






## part of txTracks
# .getPromoterTrack <- function(tssModel,
#                               targetGene
#                               
# ){
#     tss <- crisprDesign::queryTss(
#         tssObject=tssModel,
#         queryColumn="gene_symbol",
#         queryValue=targetGene
#     )
#     if (length(tss) == 0){
#         return(NULL)
#     }
#     
#     promoterTrack <- lapply(seq_along(tss), function(x){
#         id <- S4Vectors::mcols(tss)[["ID"]][x]
#         strand <- as.character(BiocGenerics::strand(tss)[x])
#         directionalId <- ifelse(strand == "+",
#                                 paste(id, "->"),
#                                 paste("<-", id))
#         Gviz::AnnotationTrack(
#             start=BiocGenerics::start(tss)[x],
#             end=BiocGenerics::end(tss)[x],
#             chromosome=GenomeInfoDb::seqnames(tss)[x],
#             strand=strand,
#             id=directionalId,
#             name=id,
#             group=directionalId,
#             showID=TRUE,
#             stacking="full",
#             collapse=FALSE,
#             col="purple",
#             fill="orchid",
#             shape="box",
#             fontcolor.group="black",
#             fontsize.group=18,
#             cex.title=0.8
#         )
#     })
#     return(promoterTrack)
# }





.getGeneTrack <- function(geneModel,
                          targetGene,
                          from,
                          to
                        
){
    if (is.null(geneModel) || is.null(targetGene)){ # need more checks (that return errors)
        return()
    }
    
    transcripts <- crisprDesign::queryTxObject(
        txObject=geneModel,
        featureType="transcripts",
        queryColumn="gene_symbol",
        queryValue=targetGene
    )
    if (length(transcripts) == 0){
        warning("targetGene not found in geneModel")
        return()
    }
    ## get tx_id whose start/end exceed limits (for each limit)
    txLeft <- transcripts$tx_id[BiocGenerics::start(transcripts) < from]
    txRight <- transcripts$tx_id[BiocGenerics::end(transcripts) > to]
    extendDistance <- to - from
    
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
    
    ## for tx_id above that exceed limit, append 1bp exon beyond limit
    exonLeft <- lapply(txLeft, function(x){
        strand <- unique(BiocGenerics::strand(splicings))
        currentRanks <- splicings$exon_rank[splicings$tx_id == x]
        if (strand == "+"){
            exon_rank <- max(1, min(currentRanks) - 1)
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
        if (strand == "+"){
            exon_rank <- max(currentRanks) + 1
        } else {
            exon_rank <- max(1, min(currentRanks) - 1)
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
        splicings$cds_start[i] <- max(BiocGenerics::start(splicings[i]), splicings$cds_start[i])
        splicings$cds_end[i] <- min(BiocGenerics::end(splicings[i]), splicings$cds_end[i])
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
        fill="darkorange",
        geneSymbol=TRUE,
        cex=2,
        cex.title=1,
        transcriptAnnotation="transcript",
        fontcolor.group="black",
        col.line="black",
        fontsize.group=20,
        min.height=12
    )
    
    txLabel <- paste0(targetGene, " (", Gviz::transcript(txTrack), ")")
    txStrand <- unique(transcripts$tx_strand)
    if (txStrand == "+"){
        txLabel <- paste(txLabel, "->")
    } else {
        txLabel <- paste("<-", txLabel)
    }
    Gviz::transcript(txTrack) <- txLabel
    
    return(txTrack)
}








#' @importFrom crisprBase getProtospacerRanges getPamRanges
#' @importFrom crisprDesign crisprNuclease
#' @importFrom GenomicRanges punion
#' @importFrom methods as
#' @importFrom S4Vectors mcols<-
#' @importFrom BiocGenerics start end strand
#' @importFrom GenomeInfoDb seqnames seqinfo seqlengths
#' @importFrom GenomicFeatures makeTxDb
#' @importFrom Gviz GeneRegionTrack
.getGuideTrack <- function(guideSet,
                           stacking="squish"
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
    txdb <- suppressWarnings(
        GenomicFeatures::makeTxDb(
            transcripts=tx,
            splicings=splicings,
            chrominfo=chrominfo
        )
    )
    
    Gviz::GeneRegionTrack(
        range=txdb,
        name="gRNAs",
        showID=TRUE,
        fill="lightblue",
        geneSymbol=TRUE,
        cex=2,
        cex.title=1,
        transcriptAnnotation="transcript",
        fontcolor.group="black",
        col.line="black",
        fontsize.group=20,
        min.height=12
    )
}







.getSequenceTrack <- function(chr,
                              bsgenome
){
    Gviz::SequenceTrack(
        sequence=bsgenome,
        chromosome=chr,
        add53=TRUE
    )
}









## other ======================================================================



.getPlotExtend <- function(extend,
                           tracks
){
    if (is.null(extend)){
        return(0)
    }
    stopifnot("extend must be a single, non-negative number" = {
        is.numeric(extend) &&
            length(extend) == 1 &&
            extend >= 0
    })
    if (extend != round(extend)){
        plotMin <- lapply(tracks, BiocGenerics::start)
        plotMin <- min(unlist(plotMin))
        plotMax <- lapply(tracks, BiocGenerics::end)
        plotMax <- max(unlist(plotMax))
        extend <- extend * (plotMax - plotMin)
        extend <- round(extend)
    }
    return(extend)
}
