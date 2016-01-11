\name{plotMotifOccurrenceAverage}
\docType{methods}
\alias{plotMotifOccurrenceAverage}
\alias{plotMotifOccurrenceAverage,DNAStringSet,matrix-method}
\title{Plotting average profile of motif occurrence}
\description{
Plots average profile of motif occurrence for a set of input sequences of the
same length. Motif is specified by a position weight matrix (PWM) that contains
estimated probability of base b at position i, and only motif hits above
specified threshold are taken into account.
}
\usage{
plotMotifOccurrenceAverage(regionsSeq, motifPWM, minScore = "80\%", flankUp = 
    NULL, flankDown = NULL, smoothingWindow = 1, color = "black", xLabel = 
    "Distance to reference point (bp)", yLabel = "Relative frequency", cexAxis = 
    1, addReferenceLine = TRUE, plotLegend = FALSE, cexLegend = 1, add = FALSE, 
    ...)
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        for which the patterns occurrence profile should be visualised.
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
        Value specifying the color for plotting.
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
    The function finds all hits matching the motif above the specified score    
    threshold in the set of input sequences and plots an average profile 
    reflecting the occurrence of the motif across input sequences.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{motifScanHits}}
    \cr
    \code{\link{plotMotifDensityMap}}
}
\examples{
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))
load(system.file("data", "TBPpwm.RData", package="seqPattern"))

plotMotifOccurrenceAverage(regionsSeq = zebrafishPromoters, motifPWM = TBPpwm,
                    minScore = "85\%", flankUp = 400, flankDown = 600, 
                    smoothingWindow = 3)
}
