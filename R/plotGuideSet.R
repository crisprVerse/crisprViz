
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






plotGuideSet <- function(guideSet, # change to guideSets...force guideStacking to be dense?
                         geneModel=NULL,
                         targetGene=NULL, # gene symbol or id
                         annotations=list(), # named list of annotations to plot (GRanges)
                         from=NULL, # see ?plotTracks
                         to=NULL, # see ?plotTracks
                         extend.left=0, # see ?plotTracks
                         extend.right=0, # see ?plotTracks
                         bands=NULL, # see ?IdeogramTrack
                         guideStacking="squish", # how to stack all guides on a single track (squish, dense, or hide) or NA to have each guide occupy a separate track
                         bsgenome=NULL, # used in SequenceTrack
                         pamSiteOnly=FALSE, # whether to plot only the PAM site for spacers or plot the spacer and PAM sequence (default)
                         onTargetScore=NULL, # color coding guides...need to add legend
                         includeSNPTrack=TRUE
                         # displayOptions=NULL # named list of options that can be passed to subfunctions
){
    ## check inputs
    guideSet <- .validateGuideSet(guideSet)
    chr <- .validateSingleChr(guideSet)
    geneModel <- .validateGRangesList(geneModel)
    stopifnot("targetGene must be a character vector or NULL" = {
        is.vector(targetGene, mode="character") || is.null(targetGene)
    })
    stopifnot("guideStacking must be NA or one of 'squish', 'dense', or 'hide'" = {
        length(guideStacking) == 1 &&
            (
                is.na(guideStacking) ||
                    guideStacking %in% c("squish", "dense", "hide")
            )
    })
    stopifnot("pamSiteOnly must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(pamSiteOnly)
    })
    
    
    
    
    
    
    
    
    
    ## set tracks
    ideogramTrack <- .getIdeogramTrack(chr=chr,
                                       genome=unique(GenomeInfoDb::genome(guideSet)),
                                       bands=bands)
    genomeAxisTrack <- .getGenomeAxisTrack()
    guideTrack <- .getGuideTrack(
        guideSet=guideSet,
        guideStacking=guideStacking,
        pamSiteOnly=pamSiteOnly,
        onTargetScore=onTargetScore
    )
    
    ## set plot range
    plotWindowMin <- lapply(guideTrack, min)
    plotWindowMin <- min(unlist(plotWindowMin))
    plotWindowMax <- lapply(guideTrack, max)
    plotWindowMax <- max(unlist(plotWindowMax))
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
    txLeft <- transcripts$tx_id[BiocGenerics::start(transcripts) < from &
                                    BiocGenerics::end(transcripts) > from]
    txRight <- transcripts$tx_id[BiocGenerics::end(transcripts) > to &
                                     BiocGenerics::start(transcripts) < to]
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
        # fill="darkorange",
        geneSymbol=TRUE,
        # cex=2,
        # cex.title=1,
        transcriptAnnotation="transcript",
        fontcolor.group="black",
        col.line="black",
        # fontsize.group=20,
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
    lapply(seq_along(annotations), function(x){
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
}





#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom crisprDesign snps
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
    
    Gviz::AnnotationTrack(
        snpRanges,
        chromosome=chr,
        shape="box",
        fill="#CDCDCD",
        group=names(snpRanges),
        name="SNP"
    )
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
                           guideStacking,
                           pamSiteOnly,
                           onTargetScore
){
    if (pamSiteOnly){
        # function to extract pam sites as txdb coords
    } else {
        # current function
    }
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
        guideStacking <- "squish"
        name <- ""
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
        name <- "gRNAs"
    }
    
    colors <- .getScoreColor(
        guideSet=guideSet,
        onTargetScore=onTargetScore
    )
    
    guideTracks <- lapply(txdb, function(x){
        Gviz::GeneRegionTrack(
            range=x,
            name=name,
            showID=TRUE,
            fill=colors,
            geneSymbol=TRUE,
            stacking=guideStacking,
            # cex=2,
            # cex.title=1,
            transcriptAnnotation="transcript",
            fontcolor.group="black",
            col.line="black"#,
            # fontsize.group=20,
            # min.height=12
        )
    })
    
    return(guideTracks)
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


.getScoreColor <- function(guideSet,
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
    colorFunction <- colorRampPalette(c(score0, score1))
    colors <- colorFunction(100)
    intervals <- cut(scores, breaks=seq(0, 1, by=0.01))
    intervals <- as.numeric(intervals)
    return(colors[intervals])
}


