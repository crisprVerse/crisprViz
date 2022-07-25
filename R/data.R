#' Data.frame of human (ideogram) bands
#'
#' A data.frame of human (ideogram) bands (hg38)
#' 
#' @format Data.frame with 5 columns:
#'     \code{chrom}, \code{chromStart}, \code{chromEnd}, \code{name},
#'     and \code{gieStain}.
#'     
#' @usage data(Ideogram_GRCh38, package="crisprViz")
"Ideogram_GRCh38"



#' Data.frame of mouse (ideogram) bands
#'
#' A data.frame of mouse (ideogram) bands (hg38)
#' 
#' @format Data.frame with 5 columns:
#'     \code{chrom}, \code{chromStart}, \code{chromEnd}, \code{name},
#'     and \code{gieStain}.
#'     
#' @usage data(Ideogram_GRCm38, package="crisprViz")
"Ideogram_GRCm38"



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
