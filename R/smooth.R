#' Smooth a heatmap
#'
#' @param heatmap A heatmap object
#' @param sigma Numeric, length 2
#' @param output.ratio  Numeric, length 2
#' @param output.size  Numeric, length 2
#' @param method One of "auto", "kernel" or "blur"
#'
#' This function smooths a heatmap using either binned kernel density
#' (more efficient for binary heatmaps) or gaussian blur.
#'
#' Sigma controls
#' the SD of the kernel in both cases, defined in terms of pixels. This means
#' that if you have very diffirent x and y dimensions (eg. a 200bp heatmap around
#' 10000 promoters) you will need to compensate by setting sigma[2] higher to get
#' the same visual effect in both dimensions
#'
#' The output size can be determined by output.ratio or output.size. "output.ratio"
#' specifies the size of the output as a negative scaling factor, so c(2, 4) would
#' result in an image half as wide and a quarter as high. "output.size" specifies
#' the dimensions of the output matrix explicitly. One or the other should be used:
#' if they are in conflict output.ratio will override output.size and there will be
#' a warning.
#'
#' Smoothing can use either a kernel density estimate or a blurring function.
#' The methods implemented are KernSmooth:bkde2D and spatstat::blur. The kernel
#' based method assumes we are smoothing individual points so the value of these points
#' are ignored. This is most useful for smoothing PatternHeatmaps where each cell in
#' the matrix is either 1 or 0. For non-binary heatmaps, blur is most appropriate. The
#' "auto" method will choose "kernel" for binary heatmaps and "blur" for any others.
#'
#' Scaling the output heatmap is handled as in CoverageHeatmap.
#'
#' @export
#' @examples
#' hm_smoothed = smooth(hm, sigma=c(5,5), output.ratio=c(2,2), method="blur")
setGeneric("smooth", function(heatmap, ...) {
    StandardGeneric("smooth")
})

#' @describeIn smooth Smooth a heatmap
#' @export
#' @importFrom spatstat blur
#' @importFrom EBImage resize
#' @importFrom KernSmooth bkde2D
setMethod("smooth", signature(heatmap="Heatmap"),
    function(heatmap, sigma=c(3,3), output.ratio=c(1,1), output.size=NULL,
             method=c("auto", "kernel", "blur")) {

    if(!is.null(output.size)) {
        if (!length(output.size) == 2) stop("output.size must have length 2")
        if (!all(output.size %% 1 == 0)) stop("output.size must have integer values")
    }

    method = match.arg(method)
    is_binary = all(image(heatmap) %in% c(0,1))

    if (method=="auto") {
        if (is_binary) {
            method = "kernel"
        } else {
            method = "blur"
        }
    }

    if (method == "kernel" && !is_binary) {
        warning("kernel method expects a binary matrix: non-binary values will be coerced to 1 or 0")
    }

    if (!all(output.ratio == c(1,1))) {
        if(length(output.ratio) != 2) stop("output.ratio must have length 2")
        if (!is.null(output.size)) warning("output.ratio is set, overiding output.size")
        output.size = rev(dim(image(heatmap))/output.ratio)
        resize_img = TRUE
    } else {
        if (is.null(output.size)) {
            resize_img = FALSE
            output.size = rev(dim(image(heatmap)))
        } else {
            resize_img = TRUE
        }
    }

    if (method=="kernel") {
        message("\nCalculating kernel density...")
        sm = as(image(heatmap), "sparseMatrix")
        df = summary(sm)
        map = bkde2D(cbind(df$i, df$j), bandwidth=sigma, gridsize=output.size,
                    range.x=list(range(ym(heatmap)), range(xm(heatmap))))
        mat.new = sum(image(heatmap))*map$fhat
        scale = get_scale(min(mat.new), max(mat.new))
    } else if (method=="blur") {
        message("\nApplying Gaussian blur...")
        mat.new = as.matrix(blur(im(image(heatmap)), sigma=sigma))
        if (resize_img == TRUE) {
            mat.new = resize(mat.new, output.size[1], output.size[2])
        }
        scale = get_scale(min(mat.new), max(mat.new))
    }

    image(heatmap) = mat.new
    scale(heatmap) = scale
    return(heatmap)
})

get_scale = function(x, y) {
    if (x >= 0 && y >=0) {
        scale = c(0, max(x, y))
    } else if (x <= 0 && y <= 0) {
        scale = c(min(x, y), 0)
    } else {
        max_abs = max(abs(c(x, y)))
        scale = c(-max_abs, max_abs)
    }
    scale
}

