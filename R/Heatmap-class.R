setClass("Heatmap",
         slots=c(
            xm="numeric",
            ym="numeric",
            value="matrix",
            transform="function",
            max_value="numeric",
            coords="integer",
            nseq="integer",
            label="character",
            metadata="list"),
        prototype=list(
            xm=numeric(0),
            ym=numeric(0),
            value=matrix(0, ncol=0, nrow=0),
            transform=NULL,
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
    cat("obs:", length(object@xm), "\n")
    cat("max_value:", object@max_value, "\n")
    cat("label:", object@label, "\n")
    cat("has_transform:", ifelse(is.null(object@transform), "yes", "no"), "\n")
})

