#' Smooth a heatmap
#'
#' @param heatmap A heatmap object
#' @param sigma Numeric, lengt2, (recycled if length 1)
#' @param output.size  Numeric, length 2
#' @param algorithm "kernel" or "blur"
#' @param ... additional arguments to S4 methods
#'
#' This function smooths a heatmap using either binned kernel density
#' (more efficient for binary heatmaps) or gaussian blur.
#'
#' Sigma controls the SD of the kernel in both cases, defined in terms of
#' pixels. This means that if you have very diffirent x and y dimensions (eg. a
#' 200bp heatmap around 10000 promoters) you will need to compensate by setting
#' sigma[2] higher to get the same visual effect in both dimensions
#'
#' "output.size" specifies the dimensions of the output matrix. This can be
#' useful to reduce plotting time significantly.
#'
#' Smoothing can use either a kernel density estimate or a blurring function.
#' The methods implemented are KernSmooth:bkde2D and EBImage::filter2 with a
#' gaussian filter. The kernel based method assumes we are smoothing individual
#' points so the value of these points are ignored. This is most useful for
#' smoothing PatternHeatmaps where each cell in the matrix is either 1 or 0.
#' For non-binary heatmaps, blur is most appropriate. Not setting this parameter
#' will choose the method automatically.
#'
#' Scaling the output heatmap is handled as in CoverageHeatmap.
#'
#' @return A heatmap
#'
#' @export
#' @examples
#' data(HeatmapExamples)
#' hm_smoothed = smoothHeatmap(hm, sigma=c(5,5), algorithm="blur")
setGeneric("smoothHeatmap", function(heatmap, ...) standardGeneric("smoothHeatmap"))

#' @describeIn smoothHeatmap Smooth a heatmap
#' @export
#' @importFrom EBImage resize filter2
#' @importFrom KernSmooth bkde2D
#' @importFrom methods as
#' @importClassesFrom Matrix sparseMatrix
#' @importMethodsFrom Matrix summary
setMethod("smoothHeatmap", signature(heatmap="Heatmap"),
    function(heatmap, sigma=c(3,3), output.size=dim(image(heatmap)), algorithm=NULL) {

    if (!all(output.size %% 1 == 0)) stop("output.size must have positive integer values")
    if (!all(output.size > 0)) stop("output.size must have positive integer values")
    if (length(output.size) != 2) stop("output.size must be length 2")
    if (any(sigma <= 0)) stop("sigma must be positive")
    if (length(sigma) == 1) {
        sigma = c(sigma, sigma)
    } else if (length(sigma) != 2) {
        stop("sigma must be length one or two")
    }

    default_choice = detect_algorithm(image(heatmap))

    if (is.null(algorithm)) {
        algorithm = default_choice
    }

    if (algorithm == "kernel") {
        smoother = kernelSmooth
        if (default_choice == "blur") warning("kernel method selected: values will be coerced to 1 or 0")
    } else if (algorithm == "blur") {
        smoother = gaussianFilter
    } else {
        stop("method must be one of kernel or blur")
    }

    new_matrix = smoother(image(heatmap), sigma, output.size)
    scale = getScale(min(new_matrix), max(new_matrix))

    image(heatmap) = new_matrix
    scale(heatmap) = scale

    return(heatmap)
})

detect_algorithm = function(mat) {
    if (all(mat %in% c(0, 1))) {
        alg = "kernel"
    } else {
        alg = "blur"
    }
    alg
}

kernelSmooth = function(mat, sigma, output.size) {
    message("\nCalculating kernel density...")
    sm = as(mat, "sparseMatrix")
    df = summary(sm)
    map = bkde2D(cbind(df$i, df$j), bandwidth=sigma, gridsize=output.size,
                range.x=list(c(0, dim(mat)[1]), c(0, dim(mat)[2])))
    mat.new = sum(mat)*map$fhat
    mat.new
}

gaussianFilter = function(mat, sigma, output.size) {
    message("\nApplying Gaussian blur...")
    brush = makeGaussian2D(sigma)
    mat.new = filter2(mat, brush)
    if (!all(output.size == dim(mat))) {
        mat.new = resize(mat.new, output.size[1], output.size[2])
    }
    mat.new
}

makeGaussian2D = function(sigma) {
    dim = sigma*3
    x = -dim[1]:dim[1]
    y = -dim[2]:dim[2]
    x_d = vapply(x, dnorm, numeric(1), sd=sigma[1])
    y_d = vapply(y, dnorm, numeric(1), sd=sigma[2])
    brush = outer(x_d, y_d)
    brush
}

#' Make an appropriate scale for a heatmap
#'
#' @param x Min/max values for the heatmap
#' @param y Min/max values for the heatmap
#'
#' This function takes min/max values for a heatmap and
#' generates a scale either starting, ending or centered on
#' zero.
#'
#' @return numeric, length 2, a new scale
#' @export
#' @examples
#' getScale(0.5, 5) # c(0, 5)
#' getScale(-6, -2) # c(-6, 6)
#' getScale(-6, 2) # c(-6, 6)
getScale = function(x, y) {
    eps = 10^(-10)
    if (x >= -eps && y >= -eps) {
        scale = c(0, max(x, y))
    } else if (x <= eps && y <= eps) {
        scale = c(min(x, y), 0)
    } else {
        max_abs = max(abs(c(x, y)))
        scale = c(-max_abs, max_abs)
    }
    scale
}

