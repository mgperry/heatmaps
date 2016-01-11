\name{plotMotifScanScores}
\docType{methods}
\alias{plotMotifScanScores}
\alias{plotMotifScanScores,DNAStringSet,matrix-method}
\title{Plotting heatmap of motif scanning scores}
\description{
Plots heatmap of motif scanning scores for a set of sorted sequences of the
same length in the form of a two dimensional map centered at a common reference
position. Motif is specified by a position weight matrix (PWM) that contains
estimated probability of base b at position i, and the percentage of the
maximal PWM score is plotted for each position (nucleotide) in each sequence.
}
\usage{
plotMotifScanScores(regionsSeq, motifPWM, seqOrder = c(1:length(regionsSeq)),
    flankUp = NULL, flankDown = NULL, xTicks = NULL, xTicksAt = NULL,
    xLabel = "", yTicks = NULL, yTicksAt = NULL, yLabel = "", cexAxis = 8,
    plotScale = TRUE, scaleLength = NULL, scaleWidth = 15,
    addReferenceLine = TRUE, plotColorLegend = TRUE, outFile =
    "MotifScanningScores.png", plotWidth = 2000, plotHeight = 2000)
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        for which the motif occurrence density should be visualised.
    }
    \item{motifPWM}{
        A numeric matrix representing the Position Weight Matrix (PWM), such as
        returned by \code{\link{PWM}} function. Can contain either probabilities
        or log2 probability ratio of base b at position i.
    }
    \item{seqOrder}{
        Integer vector specifying the order of the provided input sequences.
        Must have the same length as the number of sequences in the
        \code{regionSeq}. Input sequences will be sorted according to this index
        in an ascending order form top to the bottom of the plot, \emph{i.e.}
        the sequence labeled with the lowest number will appear at the top of
        the plot. The default value will order the sequences as they are ordered
        in the input \code{regionSeq} object.
    }
    \item{flankUp, flankDown}{
        The number of base-pairs upstream and downstream of the reference
        position in the provided sequences, respectively.
        \code{flankUp + flankDown} must sum up to the length of the sequences.
    }
    \item{xTicks}{
        Character vector of labels to be placed at the tick-marks on x-axis.
        The default \code{NULL} value produces five tick-marks: one at the
        reference point and two equally spaced tick-marks both upstream and
        downstream of the reference point.
    }
    \item{xTicksAt}{
        Numeric vector of positions of the tick-marks on the x-axis. The values
        can range from 1 (the position of the first base-pair in the sequence)
        to input sequence length. The default \code{NULL} value produces five
        tick-marks: one at the reference point and two equally spaced tick-marks
        both upstream and downstream of the reference point.
    }
    \item{xLabel}{
        The label for the x-axis. The default is no label, \emph{i.e.} empty
        string.
    }
    \item{yTicks}{
        Character vector of labels to be placed at the tick-marks on y-axis.
        The default \code{NULL} value produces no tick-marks and labels.
    }
    \item{yTicksAt}{
        Numeric vector of positions of the tick-marks on the y-axis. The values
        can range from 1 (the position of the last sequence on the bottom of the
        plot) to input sequence length (the position of the first sequence on
        the top of the plot). The default \code{NULL} value produces no
        tick-marks.
    }
    \item{yLabel}{
        The label for the y-axis. The default is no label, \emph{i.e.} empty
        string.
    }
    \item{cexAxis}{
        The magnification to be used for axis annotation.
    }
    \item{plotScale}{
        Logical, should the scale bar be plotted in the lower left corner of
        the plot.
    }
    \item{scaleLength}{
        The length of the scale bar to be plotted. Used only when
        \code{plotScale = TRUE}. If no value is provided, it defaults to one
        fifth of the input sequence length.
    }
    \item{scaleWidth}{
        The width of the line for the scale bar. Used only when
        \code{plotScale = TRUE}.
    }
    \item{addReferenceLine}{
        Logical, should the vertical dashed line be drawn at the reference
        point.
    }
    \item{plotColorLegend}{
        Logical, should the color legend for the scanning score be plotted on
        the right side of the plot.
    }
    \item{outFile}{
        Character vector specifying the base name of the output plot file. The
        final name of the plot file for each pattern will be
        \code{outFile}."pattern.jpg".
    }
    \item{plotWidth, plotHeight}{
        Width and height of the density plot in pixels.
    }
}

\value{
    The function produces a PNG file in the working directory, visualising motif
    scanning scores in the set of ordered input sequences.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{motifScanScores}}
    \cr
    \code{\link{plotPatternDensityMap}}
}
\examples{
library(GenomicRanges)
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))
promoterWidth <- elementMetadata(zebrafishPromoters)$interquantileWidth

load(system.file("data", "TBPpwm.RData", package="seqPattern"))

plotMotifScanScores(regionsSeq=zebrafishPromoters, motifPWM = TBPpwm,
        seqOrder=order(promoterWidth), flankUp = 400, flankDown = 600)
}
