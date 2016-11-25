library(heatmaps)
library(Biostrings)
context("test constructors")

pwm = matrix(runif(16), nrow=4)
pwm = log10(pwm/colSums(pwm))
rownames(pwm) = c("A", "C", "G", "T")

test_that("PWMScanHeatmap initialise default values", {
    seq = DNAStringSet("ATCGTACAGCATACAGACAT")
    hm = PWMScanHeatmap(seq, pwm)
    expect_equal(coords(hm), c(0, 20))
    expect_equal(label(hm), "")
})


