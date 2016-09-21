#' An S4 class to represent a heatmap
#'
#' @slot image A numeric Matrix
#' @slot scale A length-2 vector
#' @slot coords A length-2 integer vector
#' @slot nseq An integer
#' @slot label A character vector
#' @slot metadata A list containing arbitrary metadata
#'
#' A class used to represent a heatmap in a simple, self-contained way
#'
#' @export
#' @seealso CoverageHeatmap PatternHeatmap plotHeatmap plotHeatmapMeta
#'
#' @examples
#'
#' data(HeatmapExamples)
#'
#' hm = new("Hetamap",
#'          image=example_matrix,
#'          scale=c(0,max(example_matrix),
#'          coords=c(-500, 500),
#'          nseq=1000,
#'          label="Test",
#'          metadata=list()),
#'
#' # or use the constructor:
#' hm = Heatmap(example_matrix, coords=c(-500, 500), label="Test")
setClass("Heatmap",
         slots=c(
            image="matrix",
            scale="numeric",
            coords="integer",
            nseq="integer",
            label="character",
            metadata="list"),
        prototype=list(
            image=matrix(0, ncol=0, nrow=0),
            scale=c(0,0),
            coords=c(0L,0L),
            nseq=0L,
            label="",
            metadata=list()),
)

#' @describeIn Heatmap Return the number of sequences in a heatmap
#' @export
setMethod("length", signature="Heatmap", function(x) x@nseq)

#' @importMethodsFrom BiocGenerics width
NULL

#' @describeIn Heatmap Return the width of sequence represented in a heatmap
#' @importMethodsFrom BiocGenerics width
#' @export
setMethod("width", signature="Heatmap", function(x) x@coords[2] - x@coords[1])

setMethod("show", signature="Heatmap", function(object) {
    cat("coords:", object@coords, "\n")
    cat("nseq:", object@nseq, "\n")
    cat("obs:", length(object@matrix), "\n")
    cat("scale:", object@scale, "\n")
    cat("label:", object@label, "\n")
})

#' Reflect a heatmap in the y axis
#'
#' @param x A heatmap
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' mirror(hm)
setGeneric("mirror", def=function(x) standardGeneric("mirror"))

#' @describeIn mirror Heatmap method
#' @export
setMethod("mirror", signature="Heatmap", function(x) {
    image(x) = t(apply(image(x), 1, rev))
    x@coords = -rev(x@coords)
    x
})

#' Reflect a heatmap in the x axis
#'
#' @param x A heatmap
#' @export
setMethod("rev", signature="Heatmap", function(x) {
    image(x) = apply(image(x), 2, rev)
    x
})

#' Generate co-ordinates for each row of the image matrix of a Heatmap
#' @param x A Heatmap
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' xm(hm)
setGeneric("xm", def=function(x) standardGeneric("xm"))

#' @describeIn xm Generate co-ordinates for each frow of the image matrix of a Heatmap
#' @export
setMethod("xm", signature="Heatmap", function(x) {
    seq(1, width(x), length.out=ncol(image(x)))
})

#' Generate co-ordinates for each column of the image matrix of a Heatmap
#' @param x A Heatmap
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' ym(hm)
setGeneric("ym", def=function(x) standardGeneric("ym"))

#' @describeIn ym Generate co-ordinates for each column of the matrix
#' @export
setMethod("ym", signature="Heatmap", function(x) {
    seq(1, x@nseq, length.out=nrow(image(x)))
})

#' Return or set the scale in a Heatmap
#' @name scale
NULL

#' @rdname scale
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' scale(hm) = c(-1000, 1000)
setGeneric("scale", def=function(x) standardGeneric("scale"))

#' @rdname scale
#' @export
setMethod("scale", signature="Heatmap", function(x) {
    x@scale
})

#' @rdname scale
#' @export
setGeneric("scale<-", def=function(x, value, ...) standardGeneric("scale<-"))

#' @rdname scale
#' @export
setMethod("scale<-", signature="Heatmap", function(x, value) {
    x@scale = value
    x
})

#' Return or set the image in a Heatmap
#' @name image
NULL

#' @rdname image
#' @export
setMethod("image", signature="Heatmap", function(x) {
    x@matrix
})

#' @rdname image
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' image(hm = log(image(hm)
#' scale(hm) = c(0, max(image(hm)))
setGeneric("image<-", def=function(x, value, ...) standardGeneric("image<-"))

#' @rdname image
#' @export
setMethod("image<-", signature="Heatmap", function(x, value) {
    x@matrix = value
    x
})

#' Function to create a heatmap object
#'
#' @param image A numeric Matrix
#' @param scale A length-2 vector
#' @param coords A length-2 integer vector
#' @param nseq An integer
#' @param label A character vector
#' @param metadata A list containing arbitrary metadata
#'
#' Using this function avoids calling 'new' directly or manually setting coords and
#' nseq to integers. Other constructors exist for creating heatmaps from data, rather
#' than a raw matrix.
#'
#' @seealso PatternHeatmap CoverageHeatmap PWMScanHeatmap
#' @export
#'
#' @examples
#'
#' library(HeatmapExamples)
#' hm = Heatmap(example_matrix, coords=c(-500, 500), label="Test")
Heatmap = function(mat, coords=NULL, label="", nseq=NULL, scale=max(mat), metadata=list()) {
    if (is.null(coords)) {
        coords = c(0L, ncol(mat))
    } else {
        coords = as.integer(coords)
    }
    if (is.null(nseq)) {
        nseq = nrow(mat)
    } else {
        nseq = as.integer(nseq)
    }
    hm = new(
        "Heatmap",
        matrix=mat,
        scale=scale,
        coords=coords,
        nseq=nseq,
        label=label,
        metadata=metadata)
    hm
}


