#' Example GuideSet targeting the human KRAS gene.
#'
#' An example \linkS4class{GuideSet} object targeting the CDS of the human KRAS
#'     gene. Generated using the SpCas9 CRISPR nuclease, the genome sequence
#'     from the \code{BSgenome.Hsapiens.UCSC.hg38} package, and the gene model
#'     from \code{txdb_human} in the \code{crisprDesignData} package (Ensembl
#'     release 104).
#' 
#' @format A \linkS4class{GuideSet} object.
#'     
#' @usage data(krasGuideSet, package="crisprViz")
"krasGuideSet"



#' CompressedGRangesList describing the human KRAS gene.
#'
#' A \linkS4class{CompressedGRangesList} object describing the gene model of
#'     the human KRAS gene. Coordinates were subset from \code{txdb_human} in
#'     the \code{crisprDesignData} package (Ensembl release 104).
#' 
#' @format A \linkS4class{CompressedGRangesList} object.
#'     
#' @usage data(krasGeneModel, package="crisprViz")
"krasGeneModel"



#' Example GuideSet targeting the human GPR21 gene.
#'
#' An example \linkS4class{GuideSet} object targeting a section of the CDS of
#'     the human GPR21 gene. Generated using the BE4max CRISPR nuclease, the
#'     genome sequence from the \code{BSgenome.Hsapiens.UCSC.hg38} package, and
#'     the gene model from \code{txdb_human} in the \code{crisprDesignData}
#'     package (Ensembl release 104).
#' 
#' @format A \linkS4class{GuideSet} object.
#'     
#' @usage data(gpr21GuideSet, package="crisprViz")
"gpr21GuideSet"



#' CompressedGRangesList describing the human GPR21 gene.
#'
#' A \linkS4class{CompressedGRangesList} object describing the gene model of
#'     the human GPR21 gene. Coordinates were subset from \code{txdb_human} in
#'     the \code{crisprDesignData} package (Ensembl release 104).
#' 
#' @format A \linkS4class{CompressedGRangesList} object.
#'     
#' @usage data(gpr21GeneModel, package="crisprViz")
"gpr21GeneModel"



#' Example GuideSet targeting the human MMP7 gene.
#'
#' An example \linkS4class{GuideSet} object targeting the promoter region of
#'     the human MMP7 gene. Generated using the SpCas9 CRISPR nuclease, the
#'     genome sequence from the \code{BSgenome.Hsapiens.UCSC.hg38} package, and
#'     the gene model from \code{tss_human} in the \code{crisprDesignData}
#'     package (Ensembl release 104).
#' 
#' @format A \linkS4class{GuideSet} object.
#'     
#' @usage data(mmp7GuideSet, package="crisprViz")
"mmp7GuideSet"



#' CompressedGRangesList describing the human MMP7 gene.
#'
#' A \linkS4class{CompressedGRangesList} object describing the gene model of
#'     the human MMP7 gene. Coordinates were subset from \code{txdb_human} in
#'     the \code{crisprDesignData} package (Ensembl release 104).
#' 
#' @format A \linkS4class{CompressedGRangesList} object.
#'     
#' @usage data(mmp7GeneModel, package="crisprViz")
"mmp7GeneModel"



#' Example GuideSet targeting the human LTN1 gene.
#'
#' An example \linkS4class{GuideSet} object targeting the CDS of the human LTN1
#'     gene. Generated using the SpCas9 CRISPR nuclease, the genome sequence
#'     from the \code{BSgenome.Hsapiens.UCSC.hg38} package, and the gene model
#'     from \code{txdb_human} in the \code{crisprDesignData} package (Ensembl
#'     release 104).
#' 
#' @format A \linkS4class{GuideSet} object.
#'     
#' @usage data(cas9GuideSet, package="crisprViz")
"cas9GuideSet"



#' Example GuideSet targeting the human LTN1 gene.
#'
#' An example \linkS4class{GuideSet} object targeting the CDS of the human LTN1
#'     gene. Generated using the AsCas12a CRISPR nuclease, the genome sequence
#'     from the \code{BSgenome.Hsapiens.UCSC.hg38} package, and the gene model
#'     from \code{txdb_human} in the \code{crisprDesignData} package (Ensembl
#'     release 104).
#' 
#' @format A \linkS4class{GuideSet} object.
#'     
#' @usage data(cas12aGuideSet, package="crisprViz")
"cas12aGuideSet"



#' CompressedGRangesList describing the human LTN1 gene.
#'
#' A \linkS4class{CompressedGRangesList} object describing the gene model of
#'     the human LTN1 gene. Coordinates were subset from \code{txdb_human} in
#'     the \code{crisprDesignData} package (Ensembl release 104).
#' 
#' @format A \linkS4class{CompressedGRangesList} object.
#'     
#' @usage data(ltn1GeneModel, package="crisprViz")
"ltn1GeneModel"



#' Subset of repeat elements for hg38.
#'
#' A \linkS4class{GRanges} object describing a subset of repeat elements for
#'     the hg38 genome. Coordinates were subset from \code{gr.repeats.hg38} in
#'     the \code{crisprDesignData} package.
#' 
#' @format A \linkS4class{GRanges} object.
#'     
#' @usage data(repeats, package="crisprViz")
"repeats"



#' CAGE peak annotation from AnnotationHub
#'
#' A sample of CAGE peak annotation from AnnotationHub (ascension number
#'     AH5084), lifted over to hg38.
#' 
#' @format A \linkS4class{GRanges} object.
#'     
#' @usage data(cage, package="crisprViz")
"cage"



#' DNase I hypersensitive site annotation from AnnotationHub
#'
#' A sample of DNase I hypersensitive site annotation from AnnotationHub
#'     (ascension number AH30743), lifted over to hg38.
#' 
#' @format A \linkS4class{GRanges} object.
#'     
#' @usage data(dnase, package="crisprViz")
"dnase"
