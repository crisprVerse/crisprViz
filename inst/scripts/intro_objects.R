## generate objects used for intro vignette
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(crisprDesign)
## crisprDesignData can be obtained from github: devtools::install_github("Jfortin1/crisprDesignData")
library(crisprDesignData)
data("txdb_human", package="crisprDesignData")
data("tss_human", package="crisprDesignData")
library(usethis)


## targeting KRAS =============================================================
krasCDS <- queryTxObject(txdb_human,
                         featureType="cds",
                         queryColumn="gene_symbol",
                         queryValue="KRAS")
krasGuideSet <- findSpacers(krasCDS,
                            crisprNuclease=SpCas9,
                            bsgenome=BSgenome.Hsapiens.UCSC.hg38)
krasGuideSet <- unique(krasGuideSet)
krasGuideSet <- addOnTargetScores(krasGuideSet,
                                  methods="deephf")
use_data(krasGuideSet, compress="xz", overwrite=TRUE)

krasGeneModel <- lapply(txdb_human, function(x){
    x[x$gene_symbol == "KRAS"]
})
krasGeneModel <- as(krasGeneModel, "CompressedGRangesList")
use_data(krasGeneModel, compress="xz", overwrite=TRUE)


## targeting GPR21 ============================================================
gpr21 <- queryTxObject(txdb_human,
                       featureType="cds",
                       queryColumn="gene_symbol",
                       queryValue="GPR21")
gpr21GuideSet <- findSpacers(gpr21,
                             crisprNuclease=BE4max,
                             bsgenome=BSgenome.Hsapiens.UCSC.hg38)
gpr21GuideSet <- gpr21GuideSet[1:4]
use_data(gpr21GuideSet, compress="xz", overwrite=TRUE)

gpr21GeneModel <- lapply(txdb_human, function(x){
    x[x$gene_symbol == "GPR21"]
})
gpr21GeneModel <- as(gpr21GeneModel, "CompressedGRangesList")
use_data(gpr21GeneModel, compress="xz", overwrite=TRUE)


## targeting MMP7 =============================================================
mmp7 <- queryTss(tss_human,
                 queryColumn="gene_symbol",
                 queryValue="MMP7",
                 tss_window=c(-2000, 0))
mmp7GuideSet <- findSpacers(mmp7,
                            crisprNuclease=SpCas9,
                            bsgenome=BSgenome.Hsapiens.UCSC.hg38)
## VCF of common SNPs obtained from dbSNP website: https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
## and downloaded to local machine
vcf <- "/Users/hoberecl/crisprIndices/snps/dbsnp151.grch38/00-common_all_snps_only.vcf.gz"
mmp7GuideSet <- addSNPAnnotation(mmp7GuideSet,
                                 vcf=vcf)
use_data(mmp7GuideSet, compress="xz", overwrite=TRUE)

mmp7GeneModel <- lapply(txdb_human, function(x){
    x[x$gene_symbol == "MMP7"]
})
mmp7GeneModel <- as(mmp7GeneModel, "CompressedGRangesList")
use_data(mmp7GeneModel, compress="xz", overwrite=TRUE)


## repeats / CAGE / DNase ==============================================
library(AnnotationHub)
ah <- AnnotationHub()
chain <- ah[["AH14150"]]
targetRange <- mmp7
start(targetRange) <- start(targetRange) - 5000
end(targetRange) <- end(targetRange) + 5000

data("gr.repeats.hg38", package="crisprDesignData")
repeats <- subsetByOverlaps(gr.repeats.hg38, targetRange, ignore.strand=TRUE)
use_data(repeats, compress="xz", overwrite=TRUE)

cage <- ah[["AH5084"]]
cage <- rtracklayer::liftOver(cage, chain)
cage <- unlist(cage)
cage <- subsetByOverlaps(cage, targetRange, ignore.strand=TRUE)
use_data(cage, compress="xz", overwrite=TRUE)

dnase <- ah[['AH30743']]
dnase <- rtracklayer::liftOver(dnase, chain)
dnase <- unlist(dnase)
dnase <- subsetByOverlaps(dnase, targetRange, ignore.strand=TRUE)
use_data(dnase, compress="xz", overwrite=TRUE)


## Cas9 / Cas12a GuideSets ====================================================
target <- queryTxObject(txdb_human,
                        featureType="cds",
                        queryColumn="exon_id",
                        queryValue="ENSE00001705186")
cas9GuideSet <- findSpacers(target,
                            crisprNuclease=SpCas9,
                            bsgenome=BSgenome.Hsapiens.UCSC.hg38)
cas9GuideSet <- unique(cas9GuideSet)
use_data(cas9GuideSet, compress="xz", overwrite=TRUE)

cas12aGuideSet <- findSpacers(target,
                         crisprNuclease=AsCas12a,
                         bsgenome=BSgenome.Hsapiens.UCSC.hg38)
cas12aGuideSet <- unique(cas12aGuideSet)
use_data(cas12aGuideSet, compress="xz", overwrite=TRUE)

ltn1GeneModel <- lapply(txdb_human, function(x){
    x[x$gene_symbol == "LTN1"]
})
ltn1GeneModel <- as(ltn1GeneModel, "CompressedGRangesList")
use_data(ltn1GeneModel, compress="xz", overwrite=TRUE)

