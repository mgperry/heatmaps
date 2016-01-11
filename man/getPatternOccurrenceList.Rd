\name{getPatternOccurrenceList}
\docType{methods}
\alias{getPatternOccurrenceList}
\alias{getPatternOccurrenceList,DNAStringSet-method}
\title{Occurrence of sequence patterns in a set of ordered sequences}
\description{
Finds positions of specified sequence patterns in a list of sequences of the
same length ordered by a provided index. Sequence patterns can be consensus
sequences of variable length and can contain IUPAC ambiguity code. Position of
each pattern occurrence is specified in two-dimensional matrix, \emph{i.e.} the
first coordinate provides the ordinal number of the sequence and the second
coordinate gives the position within the sequence where the pattern occurs.
}
\usage{
getPatternOccurrenceList(regionsSeq, patterns, seqOrder =
        c(1:length(regionsSeq)), useMulticore = FALSE, nrCores = NULL)
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        in which to search for the patterns.
    }
    \item{patterns}{
        Character vector specifying one or more DNA sequence patterns
        (oligonucleotides). IUPAC ambiguity codes can be used and will match any
        letter in the subject that is associated with the code.
    }
    \item{seqOrder}{
        Integer vector specifying the order of the provided input sequences.
        Must have the same length as the number of sequences in the
        \code{regionSeq}. The default value will order the sequences as they are
        ordered in the \code{regionSeq} object.
    }
    \item{useMulticore}{
        Logical, should multicore be used. \code{useMulticore = TRUE} is
        supported only on Unix-like platforms.
    }
    \item{nrCores}{
        Number of cores to use when \code{useMulticore = TRUE}. Default value
        \code{NULL} uses all detected cores.
    }
}
\details{
    This function uses the \code{\link{matchPattern}} function to find
    occurrences of given sequence patterns in a set of input sequences.
    Input sequences must all be of the same length and are ordered according to
    the index provided in the \code{seqOrder} argument, creating a \code{n * m}
    matrix, where \code{n} is the number of sequences and \code{m} is the length
    of the sequences. Positions of pattern matches in the resulting matrix are
    returned as two-dimensional coordinates.
}
\value{
    The function returns a named list with one element for each sequence pattern
    specified in the \code{patterns} argument. Each element of the list is a
    \code{data.frame} with positions of the corresponding pattern in the set of
    input sequences. The input sequences of the same length are sorted according
    to the index in \code{seqOrder} argument and the positions of pattern
    matches in the resulting \code{n * m} matrix (where \code{n} is the number
    of sequences and \code{m} is the length of the sequence) are provided. The
    \code{sequence} column in the data.frame provides the ordinal number of the
    sequence in the ordered list of sequences and the \code{position} column
    provides the start position of the pattern match within that sequence.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{plotPatternDensityMap}}
    \cr
    \code{\link{motifScanHits}}
}
\examples{
library(GenomicRanges)
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))
promoterWidth <- elementMetadata(zebrafishPromoters)$interquantileWidth

# dinucleotide patterns
patternsOccurrence <- getPatternOccurrenceList(regionsSeq = zebrafishPromoters,
                    patterns = c("TA", "GC"), seqOrder = order(promoterWidth))
names(patternsOccurrence)
head(patternsOccurrence[["GC"]])

# motif consensus sequence
patternsOccurrence <- getPatternOccurrenceList(regionsSeq = zebrafishPromoters,
                    patterns = "TATAWAWR", seqOrder = order(promoterWidth))
names(patternsOccurrence)
head(patternsOccurrence[["TATAWAWR"]])
}
