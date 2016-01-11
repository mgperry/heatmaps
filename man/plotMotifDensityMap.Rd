\name{plotMotifDensityMap}
\docType{methods}
\alias{plotMotifDensityMap}
\alias{plotMotifDensityMap,DNAStringSet,matrix-method}
\title{Plotting density maps of motif occurrence}
\description{
Plots density of motif occurrences in an ordered set of sequences of the same
length in the form of a two dimensional map centered at a common reference
position. Motif is specified by a position weight matrix (PWM) that contains
estimated probability of base b at position i, and only motif hits above
specified threshold are taken into account and plotted.
}
\usage{
plotMotifDensityMap(regionsSeq, motifPWM, minScore = "80\%",
    seqOrder = c(1:length(regionsSeq)), flankUp = NULL, flankDown = NULL,
    nBin = NULL, bandWidth = NULL, color = "blue", transf = NULL, xTicks = NULL,
    xTicksAt = NULL, xLabel = "", yTicks = NULL, yTicksAt = NULL, yLabel = "",
    cexAxis = 8, plotScale = TRUE, scaleLength = NULL, scaleWidth = 15,
    addReferenceLine = TRUE, plotColorLegend = TRUE, outFile = "DensityMap",
    plotWidth = 2000, plotHeight = 2000)
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
        If no values are provided both \code{flankUp} and \code{flankDown} are 
        set to be half of the length of the input sequences, \emph{i.e.} the 
        reference position is assumed to be in the middle of the sequences.
    }
    \item{nBin}{
        Numeric vector with two values containing the number of equally spaced
        points in each direction over which the density is to be estimated. The
        first value specifies number of bins along x-axis, \emph{i.e.} along the
        nucleotides in the sequence, and the second value specifies the number
        of bins along y-axis, \emph{i.e.} across ordered input sequences. The
        values are passed on to the \code{gridsize} argument of the
        \code{\link{bkde2D}} function to compute a 2D binned kernel density
        estimate. If \code{nBin} is not specified it will default to
        \code{c(n, m)}, where \code{n} is the number of input sequences and
        \code{m} is the length of sequences.
    }
    \item{bandWidth}{
        Numeric vector of length 2, containing the bandwidth to be used in each
        coordinate direction. The first value specifies the bandwidth along the
        x-axis, \emph{i.e.} along the nucleotides in the sequence, and the
        second value specifies the bandwidth along y-axis, \emph{i.e.} across
        ordered input sequences. The values are passed on to the
        \code{bandwidth} argument of the \code{\link{bkde2D}} function to
        compute a 2D binned kernel density estimate and are used as standard
        deviation of the bivariate Gaussian kernel. If \code{bandWidth} is not
        specified it will default to \code{c(3,3)}.
    }
    \item{color}{
        Character specifying the color palette for the density plot. One of the 
        following color palettes can be specified: \code{"blue", "brown", 
        "cyan", "gold", "gray", "green", "pink", "purple", "red"}. Please refer 
        to the vignette for the appearance of these palettes.
    }
    \item{transf}{
        The function mapping the density scale to the color scale. See Details.
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
        Logical, should the scale bar be plotted in the lower left corner of the
        plot.
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
        Logical, should the color legend for the pattern density be plotted. If
        \code{TRUE} a separate .png file named \code{outFile}."ColorLegend.png"
        will be created, showing mapping of pattern density values to colours.
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
    The function produces a PNG file in the working directory, visualising
    density of the motif occurrence above specified threshold in the set of
    ordered input sequences.
}
\references{
    Haberle \emph{et al.} (2014) Two independent transcription initiation codes
    overlap on vertebrate core promoters, \emph{Nature} \bold{507}:381-385.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{motifScanHits}}
    \cr
    \code{\link{plotPatternDensityMap}}
}
\examples{
library(GenomicRanges)
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))
promoterWidth <- elementMetadata(zebrafishPromoters)$interquantileWidth

load(system.file("data", "TBPpwm.RData", package="seqPattern"))

plotMotifDensityMap(regionsSeq = zebrafishPromoters, motifPWM = TBPpwm,
                    minScore = "85\%", seqOrder = order(promoterWidth),
                    flankUp = 400, flankDown = 600, color = "red")
}
