#takes a sparse matrix
setGeneric("smooth", function(heatmap, ...) {
    StandardGeneric("smooth")
})

#! @importFrom spatstat blur
#! @importFrom EBImage resize
#! @importFrom KernSmooth bkde2D
setMethod("smooth", signature(heatmap="Heatmap"),
    function(heatmap, sigma=c(3,3), output.ratio=c(1,1), output.size=NULL,
             method=c("auto", "kernel", "blur")) {

    if(!is.null(output.size)) {
        if (!length(output.size) == 2) stop("output.size must have length 2")
        if (!all(output.size %% 1 == 0)) stop("output.size must have integer values")
    }

    method = match.arg(method)
    is_binary = all(heatmap@matrix %in% c(0,1))

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
        output.size = rev(dim(heatmap@matrix)/output.ratio)
        resize_img = TRUE
    } else {
        if (is.null(output.size)) {
            resize_img = FALSE
            output.size = rev(dim(heatmap@matrix))
        } else {
            resize_img = TRUE
        }
    }

    if (method=="kernel") {
        message("\nCalculating kernel density...")
        sm = as(heatmap@matrix, "sparseMatrix")
        df = summary(sm)
        map = bkde2D(cbind(df$i, df$j), bandwidth=sigma, gridsize=output.size,
                    range.x=list(range(ym(heatmap)), range(xm(heatmap))))
        mat.new = sum(heatmap@matrix)*map$fhat
        max_value = max(mat.new)
    } else if (method=="blur") {
        message("\nApplying Gaussian blur...")
        mat.new = as.matrix(blur(im(heatmap@matrix), sigma=sigma))
        if (resize_img == TRUE) {
            mat.new = resize(mat.new, output.size[1], output.size[2])
        }
        max_value = max(mat.new)
    }

    heatmap@matrix = mat.new
    heatmap@max_value = max_value
    return(heatmap)
})

