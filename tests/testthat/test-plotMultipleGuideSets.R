context('plotMultipeGuideSets')
library(crisprViz)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
data(cas9GuideSet, package="crisprViz")
data(cas12aGuideSet, package="crisprViz")
data(ltn1GeneModel, package="crisprViz")

myPlotMultipleGuideSets <- function(file, ...){
    pdf(file)
    plotMultipleGuideSets(...)
    dev.off()
}
plotOutput <- c(`null device`=1L)


test_that('Produces a plot as expected with 2 gRNA tracks ', {
    file <- tempfile(fileext=".pdf")
    input <- list(SpCas9=cas9GuideSet,
                  AsCas12a=cas12aGuideSet)
    out <- myPlotMultipleGuideSets(file=file,
                                   x=input,
                                   geneModel=ltn1GeneModel,
                                   targetGene="LTN1")
    expect_true(file.exists(file))
    expect_equal(out, plotOutput)
})



test_that('Produces a plot as expected with 1 gRNA track', {
    file <- tempfile(fileext=".pdf")
    input <- list(SpCas9=cas9GuideSet)
    out <- myPlotMultipleGuideSets(file=file,
                                   x=input,
                                   geneModel=ltn1GeneModel,
                                   targetGene="LTN1")
    expect_true(file.exists(file))
    expect_equal(out, plotOutput)
})


test_that('Error when gRNA track is not a list', {
    file <- tempfile(fileext=".pdf")
    input <- cas9GuideSet
    expect_error(myPlotMultipleGuideSets(file=file,
                                   x=input,
                                   geneModel=ltn1GeneModel,
                                   targetGene="LTN1"))
})



test_that('Warning when gene misspecified', {
    file <- tempfile(fileext=".pdf")
    input <- list(SpCas9=cas9GuideSet,
                  AsCas12a=cas12aGuideSet)
    expect_warning(myPlotMultipleGuideSets(file=file,
                                   x=input,
                                   geneModel=ltn1GeneModel,
                                   targetGene="hello"))
})




