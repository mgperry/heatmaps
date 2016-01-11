\name{plotPatternOccurrenceAverage}
\docType{methods}
\alias{plotPatternOccurrenceAverage}
\alias{plotPatternOccurrenceAverage,DNAStringSet-method}
\title{Plotting average profile of sequence pattern occurrence}
\description{
Plots average profile of sequence pattern occurrence for a set of input sequences of the same length. Multiple sequence patterns can be processed at once and visualised in the same plot, allowing comparison across different patterns.
}
\usage{
plotPatternOccurrenceAverage(regionsSeq, patterns, flankUp = NULL, 
    flankDown = NULL, smoothingWindow = 1, color = rainbow(length(patterns)),
    xLabel = "Distance to reference point (bp)", yLabel = "Relative frequency",
    cexAxis = 1, addReferenceLine = TRUE, plotLegend = TRUE, cexLegend = 1,
    useMulticore = FALSE, nrCores = NULL, add = FALSE, ...)
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        for which the patterns occurrence profile should be visualised.
    }
    \item{patterns}{
        Character vector specifying one or more DNA sequence patterns
        (oligonucleotides). IUPAC ambiguity codes can be used and will match any
        letter in the subject that is associated with the code.
    }
    \item{flankUp, flankDown}{
        The number of base-pairs upstream and downstream of the reference
        position in the provided sequences, respectively.
        \code{flankUp + flankDown} must sum up to the length of the sequences. 
        If no values are provided both \code{flankUp} and \code{flankDown} are 
        set to be half of the length of the input sequences, \emph{i.e.} the 
        reference position is assumed to be in the middle of the sequences.
    }
    \item{smoothingWindow}{
        Integer specifying the size of a window (in base-pairs) to be used for 
        smoothing the signal. Default value of 1 corresponds to no smoothing.
    }
    \item{color}{
        A vector of values specifying the colors for plotting. Number of colors
        must match the number of patterns that should be plotted.
    }
    \item{xLabel, yLabel}{
        Character strings for x and y axis labels, respectively.
    }
    \item{cexAxis}{
        The magnification to be used for axis annotation relative to the current 
        setting of \code{cex}.
    }
    \item{addReferenceLine}{
        Logical, should the vertical dashed line be drawn at the reference
        point.
    }
    \item{plotLegend}{
        Logical, should the legend be plotted at the top.
    }
    \item{cexLegend}{
        The magnification to be used for legend relative to the current setting 
        of \code{cex}.
    }
    \item{useMulticore}{
        Logical, should multicore be used. \code{useMulticore = TRUE} is
        supported only on Unix-like platforms.
    }
    \item{nrCores}{
        Number of cores to use when \code{useMulticore = TRUE}.  Default value
        \code{NULL} uses all detected cores.
    }
    \item{add}{
        Logical, should the pattern occurrence profiles be added to the existing 
        plot.
    }
    \item{...}{
        Further arguments to be passed to \code{plot()} or \code{lines()} 
        methods, such as \code{lty}, \emph{etc}.
    }
}
\value{
    The function finds all hits matching the specified patterns in the set of 
    input sequences and plots one average profile per pattern reflecting the 
    occurrence of patterns across sequences.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{getPatternOccurrenceList}}
    \cr
    \code{\link{plotPatternDensityMap}}
}
\examples{
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))

# dinucleotide patterns
plotPatternOccurrenceAverage(regionsSeq = zebrafishPromoters, patterns = c("AT", 
            "TA", "CG", "GC"), flankUp = 400, flankDown = 600, smoothingWindow = 
            3, color = c("gold", "darkred", "forestgreen", "navy"))

# motif consensus sequence
plotPatternOccurrenceAverage(regionsSeq = zebrafishPromoters, patterns = 
            "TATAWAWR", flankUp = 400, flankDown = 600, smoothingWindow = 3, 
            color = "gray")
}
