\name{zebrafishPromoters24h}
\docType{data}
\alias{zebrafishPromoters24h}
\title{Zebrafish promoters genomic coordinates}
\description{
This is a \code{\link{data.frame}} object that contains genomic coordinates of 
12078 promoters active in zebrafish \emph{Danio rerio} embryos at 24 hours past 
fertilisation. For each promoter additional information is provided, including 
position of the dominant (most frequently used) TSS position, number of CAGE 
tags per million supporting that promoter and the interquantile width of the 
promoter (width of the central region containing >= 80\% of the CAGE tags). The 
data is taken from Nepal \emph{et al.} Genome Research 2013, and it is intended 
to be used for running examples from \code{\link{seqPattern}} package vignette.
}
\usage{data(zebrafishPromoters24h)}
\value{A \code{\link{data.frame}} object}
\format{A \code{\link{data.frame}} object}
\keyword{datasets}
\references{
Nepal \emph{et al.} (2013) Dynamic regulation of the transcription initiation 
landscape at single nucleotide resolution during vertebrate embryogenesis, 
\emph{Genome Research} \bold{23}(11):1938-1950.
}