context('plotGuideSet')
library(crisprViz)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
data(krasGuideSet, package="crisprViz")
data(krasGeneModel, package="crisprViz")



myPlotGuideSet <- function(file,...){
    pdf(file)
    plotGuideSet(...)
    dev.off()
}
plotOutput <- c(`null device`=1L)


test_that('Produces a plot as expected with one gRNA', {
    file <- tempfile(fileext=".pdf")
    out <- myPlotGuideSet(file=file,
                          x=krasGuideSet[1],
                          geneModel=krasGeneModel,
                          targetGene="KRAS")
    expect_true(file.exists(file))
    expect_equal(out, plotOutput)
})


test_that('Produces a plot as expected with multiple gRNAs', {
    file <- tempfile(fileext=".pdf")
    out <- myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          targetGene="KRAS")
    expect_true(file.exists(file))
    expect_equal(out, plotOutput)
})


test_that('Warning when gene misspecified', {
    file <- tempfile(fileext=".pdf")
    expect_warning(myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          targetGene="A"))
})


test_that('Testing pamSiteOnly', {
    file <- tempfile(fileext=".pdf")
    out <- myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          pamSiteOnly=TRUE)
    expect_true(file.exists(file))
})

test_that('Testing onTargetScore option', {
    file <- tempfile(fileext=".pdf")
    out <- myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          onTargetScore="score_deephf")
    expect_true(file.exists(file))
})


test_that('Testing onTargetScore misspecification', {
    file <- tempfile(fileext=".pdf")
    expect_error(myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          onTargetScore="score_not_existing"))
})


test_that('Testing gcWindow misspecification', {
    file <- tempfile(fileext=".pdf")
    expect_error(myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          gcWindow=-100))
})


test_that('Testing required bsgenome for gcWindow specification', {
    file <- tempfile(fileext=".pdf")
    expect_error(myPlotGuideSet(file=file,
                          x=krasGuideSet[1:4],
                          geneModel=krasGeneModel,
                          gcWindow=100))
})


test_that('Testing gcWindow specification', {
    file <- tempfile(fileext=".pdf")
    out <- myPlotGuideSet(file=file,
                          x=krasGuideSet[1],
                          geneModel=krasGeneModel,
                          bsgenome=bsgenome,
                          gcWindow=2)
    expect_true(file.exists(file))
})



test_that('Testing gRNA labels turned off', {
    file <- tempfile(fileext=".pdf")
    out <- myPlotGuideSet(file=file,
                          x=krasGuideSet[1],
                          geneModel=krasGeneModel,
                          showGuideLabels=FALSE)
    expect_true(file.exists(file))
})






