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
#' Slots can be accessed and set using getters and setters with the same name.
#'
#' @export
#' @seealso CoverageHeatmap PatternHeatmap plotHeatmap plotHeatmapMeta
#'
#' @examples
#'
#' data(HeatmapExamples)
#'
#' hm = new("Heatmap",
#'          image=mat,
#'          scale=c(0,max(mat)),
#'          coords=c(-100L, 100L),
#'          nseq=1000L,
#'          label="Test",
#'          metadata=list())
#'
#' # or use the constructor:
#' hm = Heatmap(mat, coords=c(-100, 100), label="Test")
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
        validity=function(object) {
            errors = character()
            if (length(scale(object)) != 2) {
                errors = c(errors, "Scale should be length 2")
            }
            if (length(coords(object)) != 2) {
                errors = c(errors, "coords should be length 2")
            }
            if (length(nseq(object)) != 1) {
                errors = c(errors, "nseq should be length 1")
            }
            if (length(label(object)) != 1) {
                errors = c(errors, "label should be length 1")
            }
            if (length(errors) == 0) TRUE else errors }
)

#' Return the number of sequences in a heatmap
#' @param x A heatmap
#' @return integer, value of x@nseq
#' @export
setMethod("length", signature="Heatmap", function(x) x@nseq)

#' Return the width of sequence represented in a heatmap
#' @param x A heatmap
#' @return integer
#' @importMethodsFrom BiocGenerics width
#' @export
setMethod("width", signature="Heatmap", function(x) x@coords[2] - x@coords[1])

setMethod("show", signature="Heatmap", function(object) {
cat("coords:", coords(object), "\n")
cat("nseq:", nseq(object), "\n")
cat("obs:", length(image(object)), "\n")
cat("scale:", scale(object), "\n")
cat("label:", label(object), "\n")
})

#' Reflect a heatmap in the x axis
#'
#' @param x A heatmap
#' @return A heatmap
#' @export
setMethod("rev", signature="Heatmap", function(x) {
image(x) = apply(image(x), 2, rev)
x
})

#' Generate co-ordinates for each row of the image matrix of a Heatmap
#' @param x A Heatmap
#' @return numeric, a list of co-ordinates for plotting values in hm@image
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
#' @return numeric, a list of co-ordinates for plotting values in hm@image
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
#' @param x A heatmap
#' @param value Replacement value
#' @return numeric, length 2, the value of hm@scale
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
setGeneric("scale<-", def=function(x, value) standardGeneric("scale<-"))

#' @rdname scale
#' @export
setMethod("scale<-", signature="Heatmap", function(x, value) {
    x@scale = value
    x
})

#' Return or set the image in a Heatmap
#' @param x A heatmap
#' @param value Replacement value
#' @return matrix, from hm@image
#' @name image
NULL

#' @rdname image
#' @export
setMethod("image", signature="Heatmap", function(x) {
    x@image
})

#' @rdname image
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' image(hm) = log(image(hm))
#' scale(hm) = c(0, max(image(hm)))
setGeneric("image<-", def=function(x, value) standardGeneric("image<-"))

#' @rdname image
#' @export
setMethod("image<-", signature="Heatmap", function(x, value) {
    x@image = value
    x
})

#' Return or set the coords in a Heatmap
#' @param x A heatmap
#' @param value Replacement value
#' @return integer, length 2, value of x@coords
#' @name coords
NULL

#' @rdname coords
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' coords(hm) = c(-100, 100)
setGeneric("coords", def=function(x) standardGeneric("coords"))

#' @rdname coords
#' @export
setMethod("coords", signature="Heatmap", function(x) {
    x@coords
})

#' @rdname coords
#' @export
setGeneric("coords<-", def=function(x, value) standardGeneric("coords<-"))

#' @rdname coords
#' @export
setMethod("coords<-", signature="Heatmap", function(x, value) {
    x@coords = as.integer(value)
    x
})

#' Return or set the label in a Heatmap
#' @param x A heatmap
#' @param value Replacement value
#' @return character, value of hm@label
#' @name label
NULL

#' @rdname label
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' label(hm) = "NewLabel"
#' label(hm) # "NewLabel"
setGeneric("label", def=function(x) standardGeneric("label"))

#' @rdname label
#' @export
setMethod("label", signature="Heatmap", function(x) {
    x@label
})

#' @rdname label
#' @export
setGeneric("label<-", def=function(x, value) standardGeneric("label<-"))

#' @rdname label
#' @export
setMethod("label<-", signature="Heatmap", function(x, value) {
    x@label = value
    x
})

#' Return or set nseq in a Heatmap
#' @param x A heatmap
#' @param value Replacement value
#' @return integer, value of hm@nseq
#' @name nseq
NULL

#' @rdname nseq
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' nseq(hm) = 1000
setGeneric("nseq", def=function(x) standardGeneric("nseq"))

#' @rdname nseq
#' @export
setMethod("nseq", signature="Heatmap", function(x) {
    x@nseq
})

#' @rdname nseq
#' @export
setGeneric("nseq<-", def=function(x, value) standardGeneric("nseq<-"))

#' @rdname nseq
#' @export
setMethod("nseq<-", signature="Heatmap", function(x, value) {
    x@nseq = as.integer(value)
    x
})

#' Return or set the metadata in a Heatmap
#'
#' Store arbitrary metadata in a list, if desired.
#'
#' @param x A heatmap
#' @param value Replacement value
#' @return list, value of hm@metadata
#' @name metadata
NULL

#' @rdname metadata
#' @export
#' @examples
#'
#' data(HeatmapExamples)
#' metadata(hm) = list(replicate=1, cell_line="ESC")
#' metadata(hm)$replicate == 1
setGeneric("metadata", def=function(x) standardGeneric("metadata"))

#' @rdname metadata
#' @export
setMethod("metadata", signature="Heatmap", function(x) {
    x@metadata
})

#' @rdname metadata
#' @export
setGeneric("metadata<-", def=function(x, value) standardGeneric("metadata<-"))

#' @rdname metadata
#' @export
setMethod("metadata<-", signature="Heatmap", function(x, value) {
    x@metadata = value
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
#' @return A Heatmap object
#'
#' @seealso PatternHeatmap CoverageHeatmap PWMScanHeatmap
#' @importFrom methods new
#' @export
#'
#' @examples
#'
#' data(HeatmapExamples)
#' hm = Heatmap(mat, coords=c(-100, 100), label="Test")
Heatmap = function(image, coords=NULL, label="", nseq=NULL, scale=NULL, metadata=list()) {
    if (is.null(coords)) {
        coords = c(0L, ncol(image))
    } else {
        coords = as.integer(coords)
    }
    if (is.null(nseq)) {
        nseq = nrow(image)
    } else {
        nseq = as.integer(nseq)
    }
    if (is.null(scale)) {
        scale = getScale(min(image), max(image))
    }
    hm = new(
        "Heatmap",
        image=image,
        scale=scale,
        coords=coords,
        nseq=nseq,
        label=label,
        metadata=metadata)
    hm
}


