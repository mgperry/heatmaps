library(heatmaps)
context("smoothing heatmaps")

kernelSmooth = heatmaps:::kernelSmooth
gaussianFilter = heatmaps:::gaussianFilter
detect_algorithm = heatmaps:::detect_algorithm

rbern = function(n, p=0.5) {
    sample(c(0, 1), n, prob=c(1-p, p), replace=TRUE)
}

expect_dim = function(x, d) {
    expect_equal(dim(x), d)
}

test_that("smoothing method is selected correctly", {
    expect_equal(detect_algorithm(c(0, 1, 1, 0, 1)), "kernel")
    expect_equal(detect_algorithm(c(0, 1, 1, 0, -1)), "blur")
    expect_equal(detect_algorithm(c(0, 1, 1, 2, 0.5)), "blur")
})

test_that("kernelSmooth output size is correct", {
    mat = matrix(rbern(50), nrow=10)
    expect_dim(kernelSmooth(mat, c(3, 3), c(4, 2)), c(4, 2))
})

test_that("kernelSmooth performs some smoothing", {
    mat = matrix(rbern(50), nrow=10)
    expect_lt(max(kernelSmooth(mat, c(3, 3), c(4, 2))), 1)
})

test_that("gaussianFilter output size is correct", {
    mat = matrix(rbern(5000), nrow=100)
    expect_dim(gaussianFilter(mat, c(3, 3), c(50, 10)), c(50, 10))
})

test_that("gaussianFilter performs some smoothing", {
    mat = matrix(rbern(5000), nrow=100)
    expect_lt(max(gaussianFilter(mat, c(3, 3), c(4, 2))), 1)
})

test_that("makeGaussian2D produces the correct size output", {
    expect_dim(makeGaussian2D(c(3,5)), c(19, 31))
})

test_that("getScale works as expected", {
    expect_equal(getScale(0, 10), c(0, 10))
    expect_equal(getScale(4, 10), c(0, 10))
    expect_equal(getScale(-2, 10), c(-10, 10))
    expect_equal(getScale(-10, -2), c(-10, 0))
    expect_equal(getScale(-10, 0), c(-10, 0))
})

test_that("smooth checks inputs correctly", {
    hm = Heatmap(matrix(rbern(50), nrow=10))
    expect_error(smoothHeatmap(hm, output.size=c(0, 5)))
    expect_error(smoothHeatmap(hm, output.size=c(2.5, 5)))
    expect_error(smoothHeatmap(hm, output.size=c(2.5, 5, 10)))
    expect_error(smoothHeatmap(hm, sigma=c(5, 5, 5)))
    expect_error(smoothHeatmap(hm, sigma=c(-2)))
    expect_error(smoothHeatmap(hm, algorithm="box_blur"))
})

test_that("smoothHeatmap returns a heatmap", {
    hm = Heatmap(matrix(rbern(500), nrow=25))
    expect_true(class(smoothHeatmap(hm, sigma=c(3,3), output.size=c(8, 10))) == "Heatmap")
})

