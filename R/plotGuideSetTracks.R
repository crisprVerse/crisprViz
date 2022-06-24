





## TODO:
#   additional tracks if SNP, repeat elements detected in gs, or provided?






plotGuideSets <- function(guideSets, # list of guideSets (names optional)
                          ideogram=NULL, # genome identifier or bands
                          geneModel=NULL,
                          targetGene=NULL, # either a gene symbol/id, or genomic coordinates? ...how to implement?
                          from=NULL, # see plotTracks
                          to=NULL, # see plotTracks
                          extend.left=0, # see plotTracks
                          extend.right=0, # see plotTracks
                          legendMetric=NULL, # color coding guides...need to add legend
                          ... # options that can be passed to subfunctions?
){
    ## check inputs
    if (!is.list(guideSets)){
        stop("not a list")
    }
    lapply(guideSets, function(x){
        stopifnot("guideSets must be a list of GuideSet objects" = {
            methods::is(x, "GuideSet")
        })
    })
    chr <- lapply(guideSets, function(x){
        as.character(unique(GenomeInfoDb::seqnames(x)))
    })
    chr <- unique(unlist(chr))
    stopifnot("Can only plot a single chromosome" = {
        length(chr) == 1
    })
    
    
    ## get tracks
    ideogramTrack <- .getIdeogramTrack(chr=chr,
                                       ideogram=ideogram)
    genomeAxisTrack <- .getGenomeAxisTrack()
    guideTracks <- .getGuideTrack(guideSets)
    txTrack <- .getTxTrack(geneModel=geneModel,
                           targetGene=targetGene)
    
    
    
    ## arrange tracks
    ## get plotting range
    
    ## plot
    tracks <- c(list(ideogramTrack, genomeAxisTrack), guideTracks, txTrack)
    tracks <- tracks[vapply(tracks, function(x){
        !is.null(x)
    }, FUN.VALUE=logical(1))]
    
    # extend.left <- .getPlotExtend(extend=extend.left,
    #                               tracks=guideTracks) # add tx/tss tracks, if any
    # extend.right <- .getPlotExtend(extend=extend.right,
    #                                tracks=guideTracks) # add tx/tss tracks, if any
    
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




## tracks =====================================================================




.getIdeogramTrack <- function(chr,
                              ideogram
){
    ## check whether ideogram is genome or bands
    if (is.null(ideogram)){
        return(NULL)
    } else if (is.data.frame(ideogram)){
        ## assume is bands -- need test
        Gviz::IdeogramTrack(
            chromosome=chr,
            bands=ideogram,
            fontcolor="black",
            cex=1
        )
    } else if (is.character(ideogram) && length(ideogram) == 1){
        Gviz::IdeogramTrack(
            chromosome=chr,
            genome=ideogram,
            fontcolor="black",
            cex=1
        )
    } else {
        stop("Can't recognize ideogram arg")
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



.getGuideTrack <- function(guideSets
){
    if (is.null(names(guideSets))){
        names(guideSets) <- paste0("GuideSet_", seq_along(guideSets))
    }
    lapply(seq_along(guideSets), function(x){
        name <- names(guideSets)[x]
        protospacerRanges <- crisprBase::getProtospacerRanges(
            gr=guideSets[[x]],
            nuclease=crisprDesign::crisprNuclease(guideSets[[x]]),
        )
        Gviz::AnnotationTrack(
            range=protospacerRanges,
            chromosome=GenomeInfoDb::seqnames(protospacerRanges),
            strand=BiocGenerics::strand(protospacerRanges),
            name=name,
            id="id",
            group=names(protospacerRanges),
            showID=TRUE,
            stacking="squish",
            collapse=FALSE,
            col="darkblue",
            shape="box",
            fontcolor.group="black",
            fontsize.group=18,
            cex.title=0.8
        )
    })
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





.getTxTrack <- function(geneModel,
                        targetGene
                        
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
    splicings <- crisprDesign::queryTxObject(
        txObject=geneModel,
        featureType="exons",
        queryColumn="gene_symbol",
        queryValue=targetGene
    )
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
