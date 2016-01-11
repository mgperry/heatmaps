\name{motifScanHits}
\docType{methods}
\alias{motifScanHits}
\alias{motifScanHits,DNAStringSet,matrix-method}
\title{Occurrence of motifs in a set of ordered sequences}
\description{
Finds positions of sequence motif hits above a specified threshold in a list of
sequences of the same length ordered by a provided index. Motif is specified by
a position weight matrix (PWM) that contains estimated probability of base b at
position i and is usually constructed via call to \code{\link{PWM}} function.
Position of each motif hit is specified in two-dimensional matrix, \emph{i.e.}
the first coordinate provides the ordinal number of the sequence and the second
coordinate gives the position within the sequence where the motif occurs.
}
\usage{
motifScanHits(regionsSeq, motifPWM, minScore = "80\%",
    seqOrder = c(1:length(regionsSeq)))
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        in which to search for the motif hits.
    }
    \item{motifPWM}{
        A numeric matrix representing the Position Weight Matrix (PWM), such as
        returned by \code{\link{PWM}} function. Can contain either probabilities
        or log2 probability ratio of base b at position i.
    }
    \item{minScore}{
        The minimum score for counting a motif hit. Can be given as a character
        string containing a percentage (\emph{e.g.} \code{"85\%"}) of the
        PWM score or a single number specifying score threshold. If a percentage
        is given, it is converted to a score value taking into account both
        minimal and maximal possible PWM scores as follows:
        \code{minPWMscore + percThreshold/100 * (maxPWMscore - minPWMscore)}
        This differs from the formula in the \code{\link{matchPWM}} function
        from the \code{Biostrings} package which takes into account only the 
        maximal possible PWM score and considers the given percentage as the 
        percentage of that maximal score:
        \code{percThreshold/100 * maxPWMscore}
    }
    \item{seqOrder}{
        Integer vector specifying the order of the provided input sequences.
        Must have the same length as the number of sequences in the
        \code{regionSeq}. The default value will order the sequences as they are
        ordered in the input \code{regionSeq} object.
    }
}
\details{
    This function uses the \code{\link{matchPWM}} function to find matches to
    given motif in a set of input sequences. Only matches above specified
    \code{minScore} are considered as hits. Input sequences must all be of the
    same length and are ordered according to the index provided in the
    \code{seqOrder} argument, creating a \code{n * m} matrix, where \code{n} is
    the number of sequences and \code{m} is the length of the sequences.
    Positions of motif hits in the resulting matrix are returned as
    two-dimensional coordinates.
}
\value{
    The function returns a \code{data.frame} with positions of the motif hits in
    the set of input sequences. The input sequences of the same length are
    sorted according to the index in \code{seqOrder} argument and the positions
    of motif hits in the resulting \code{n * m} matrix (where \code{n} is the
    number of sequences and \code{m} is the length of the sequence) are
    provided. The \code{sequence} column in the data.frame provides the ordinal
    number of the sequence in the ordered list of sequences and the
    \code{position} column provides the start position of the motif hit within
    that sequence.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{plotMotifDensityMap}}
    \cr
    \code{\link{getPatternOccurrenceList}}
}
\examples{
library(GenomicRanges)
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))
promoterWidth <- elementMetadata(zebrafishPromoters)$interquantileWidth

load(system.file("data", "TBPpwm.RData", package="seqPattern"))

motifOccurrence <- motifScanHits(regionsSeq = zebrafishPromoters,
                                motifPWM = TBPpwm, minScore = "85\%",
                                seqOrder = order(promoterWidth))
head(motifOccurrence)
}
