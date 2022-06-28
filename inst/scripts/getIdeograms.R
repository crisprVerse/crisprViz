library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
data(UCSC.Mouse.GRCm38.CytoBandIdeogram)

# Human
dat <- as.data.frame(UCSC.HG38.Human.CytoBandIdeogram,
                     stringsAsFactors=FALSE)
colnames(dat) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
Ideogram_GRCh38 <- dat
save(Ideogram_GRCh38,
     file="../../data/Ideogram_GRCh38.rda")

# Mouse
dat <- as.data.frame(UCSC.Mouse.GRCm38.CytoBandIdeogram,
                     stringsAsFactors=FALSE)
colnames(dat) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
Ideogram_GRCm38 <- dat
save(Ideogram_GRCm38,
     file="../../data/Ideogram_GRCm38.rda")
