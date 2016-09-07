setClass("Heatmap",
         slots=c(
            matrix="matrix",
            scale="numeric",
            coords="integer",
            nseq="integer",
            label="character",
            metadata="list"),
        prototype=list(
            matrix=matrix(0, ncol=0, nrow=0),
            scale=c(0,0),
            coords=c(0L,0L),
            nseq=0L,
            label="",
            metadata=list()),
)

setMethod("length", signature="Heatmap", function(x) x@nseq)

setMethod("width", signature="Heatmap", function(x) x@coords[2] - x@coords[1])

setMethod("show", signature="Heatmap", function(object) {
    cat("coords:", object@coords, "\n")
    cat("nseq:", object@nseq, "\n")
    cat("obs:", length(object@matrix), "\n")
    cat("scale:", object@scale, "\n")
    cat("label:", object@label, "\n")
})

setGeneric("mirror", def=function(x) standardGeneric("mirror"))

setMethod("mirror", signature="Heatmap", function(x) {
    image(x) = t(apply(image(x), 1, rev))
    x@coords = -rev(x@coords)
    x
})

setMethod("rev", signature="Heatmap", function(x) {
    image(x) = apply(image(x), 2, rev)
    x
})

setGeneric("xm", def=function(x) standardGeneric("xm"))

setMethod("xm", signature="Heatmap", function(x) {
    seq(1, width(x), length.out=ncol(image(x)))
})

setGeneric("ym", def=function(x) standardGeneric("ym"))

setMethod("ym", signature="Heatmap", function(x) {
    seq(1, x@nseq, length.out=nrow(image(x)))
})

setGeneric("scale", def=function(x) standardGeneric("scale"))

setMethod("scale", signature="Heatmap", function(x) {
    x@scale
})

setGeneric("scale<-", def=function(x, value, ...) standardGeneric("scale<-"))

setMethod("scale<-", signature="Heatmap", function(x, value) {
    x@scale = value
    x
})

setMethod("image", signature="Heatmap", function(x) {
    x@matrix
})

setGeneric("image<-", def=function(x, value, ...) standardGeneric("image<-"))

setMethod("image<-", signature="Heatmap", function(x, value) {
    x@matrix = value
    x
})

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


