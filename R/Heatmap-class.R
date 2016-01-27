setClass("Heatmap",
         slots=c(
            xm="numeric",
            ym="numeric",
            matrix="matrix",
            max_value="numeric",
            coords="integer",
            nseq="integer",
            label="character",
            metadata="list"),
        prototype=list(
            xm=numeric(0),
            ym=numeric(0),
            matrix=matrix(0, ncol=0, nrow=0),
            max_value=0,
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
    cat("max_value:", object@max_value, "\n")
    cat("label:", object@label, "\n")
})

setGeneric("mirror", def=function(x) standardGeneric("mirror"))

setMethod("mirror", signature="Heatmap", function(x) {
    x@matrix = t(apply(x@matrix, 1, rev))
    x@coords = -rev(x@coords)
    x
})

setMethod("rev", signature="Heatmap", function(x) {
    x@matrix = apply(x@matrix, 2, rev)
    x
})

Heatmap = function(mat, coords=NULL, label="", nseq=NULL, max_value=max(mat), metadata=list()) {
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
    xm = seq(1, coords[2] - coords[1], length.out=ncol(mat))
    ym = 1:nrow(mat)
    hm = new(
        "Heatmap",
        xm=xm,
        ym=ym,
        matrix=mat,
        max_value=max_value,
        coords=coords,
        nseq=nseq,
        label=label,
        metadata=metadata)
    hm
}


