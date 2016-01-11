\name{motifScanScores}
\docType{methods}
\alias{motifScanScores}
\alias{motifScanScores,DNAStringSet,matrix-method}
\title{Motif scanning scores for a set of ordered sequences}
\description{
Provides motif scanning scores along the full length of a sequence for a list
of sequences of the same length ordered by a provided index. Motif is specified
by a position weight matrix (PWM) that contains estimated probability of base b
at position i and is usually constructed via call to \code{\link{PWM}} function.
Scanning scores are returned in the form of a two-dimensional matrix, where the
rows are sequences ordered by the specified index and the columns are relative
positions within the sequence. Each cell in the matrix contains the score of
the specified motif in the given sequence starting at the given position. The
resulting matrix can be used to visualise motif occurrences and their strength
in an ordered set of sequences centered at a common reference point.
}
\usage{
motifScanScores(regionsSeq, motifPWM, seqOrder = c(1:length(regionsSeq)),
        asPercentage = TRUE)
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        to be scanned with the motif.
    }
    \item{motifPWM}{
        A numeric matrix representing the Position Weight Matrix (PWM), such as
        returned by \code{\link{PWM}} function. Can contain either probabilities
        or log2 probability ratio of base b at position i.
    }
    \item{seqOrder}{
        Integer vector specifying the order of the provided input sequences.
        Must have the same length as the number of sequences in the
        \code{regionSeq}. The default value will order the sequences as they are
        ordered in the input \code{regionSeq} object.
    }
    \item{asPercentage}{
        Logical, should the scores represent percentage of the maximal motif PWM
        score (\code{TRUE}) or raw scores (\code{FALSE}).
    }
}
\details{
    This function uses the \code{\link{PWMscoreStartingAt}} function to get
    scores for a given motif starting at each position (nucleotide) in a set of
    input sequences. Input sequences must all be of the same length and are
    ordered according to the index provided in the \code{seqOrder} argument,
    creating an \code{n * m} matrix, where \code{n} is the number of sequences
    and \code{m} is the length of the sequences. Each cell in the matrix
    contains the score of the specified motif in the given sequence starting at
    the given position.
}
\value{
    The function returns a \code{matrix} with motif scanning scores for each
    position in the set of input sequences.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{plotMotifScanScores}}
    \cr
    \code{\link{motifScanHits}}
}
\examples{
library(GenomicRanges)
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))
promoterWidth <- elementMetadata(zebrafishPromoters)$interquantileWidth

load(system.file("data", "TBPpwm.RData", package="seqPattern"))

motifScores <- motifScanScores(regionsSeq = zebrafishPromoters,
                            motifPWM = TBPpwm, seqOrder = order(promoterWidth),
                            asPercentage = TRUE)
dim(motifScores)
motifScores[1:10,1:10]
}
