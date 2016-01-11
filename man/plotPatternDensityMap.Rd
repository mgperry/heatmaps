\name{plotPatternDensityMap}
\docType{methods}
\alias{plotPatternDensityMap}
\alias{plotPatternDensityMap,DNAStringSet-method}
\title{Plotting density maps of sequence pattern occurrence}
\description{
Plots density of sequence pattern occurrences in an ordered set of sequences of
the same length in the form of a two dimensional map centered at a common
reference position. Multiple sequence patterns can be processed at once and one
plot per pattern will be created with the same color scale across all plots,
allowing visual density comparison across different patterns.
}
\usage{
plotPatternDensityMap(regionsSeq, patterns, seqOrder = c(1:length(regionsSeq)),
    flankUp = NULL, flankDown = NULL, nBin = NULL, bandWidth = NULL, 
    color = "blue", transf = NULL, xTicks = NULL, xTicksAt = NULL, xLabel = "", 
    yTicks = NULL, yTicksAt = NULL, yLabel = "", cexAxis = 8, plotScale = TRUE,
    scaleLength = NULL, scaleWidth = 15, addPatternLabel = TRUE, cexLabel = 8,
    labelCol = "black", addReferenceLine = TRUE, plotColorLegend = TRUE,
    outFile = "PatternDensityMap", plotWidth = 2000, plotHeight = 2000,
    useMulticore = FALSE, nrCores = NULL)
}
\arguments{
    \item{regionsSeq}{
        A \code{\link{DNAStringSet}} object. Set of sequences of the same length
        for which the patterns occurrence density should be visualised.
    }
    \item{patterns}{
        Character vector specifying one or more DNA sequence patterns
        (oligonucleotides). IUPAC ambiguity codes can be used and will match any
        letter in the subject that is associated with the code.
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
        first value specifies number of bins along x-axis, i.e. along the
        nucleotides in the sequence, and the second value specifies the number
        of bins along y-axis, i.e. across ordered input sequences. The values
        are passed on to the \code{gridsize} argument of the
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
    \item{addPatternLabel}{
        Logical, should the pattern label be written in the upper left corner
        of the plot.
    }
    \item{cexLabel}{
        The magnification to be used for pattern label.
    }
    \item{labelCol}{
        The color to be used for pattern label and scale bar.
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
        \code{outFile}."pattern.png".
    }
    \item{plotWidth, plotHeight}{
        Width and height of the density plot(s) in pixels.
    }
    \item{useMulticore}{
        Logical, should multicore be used. \code{useMulticore = TRUE} is
        supported only on Unix-like platforms.
    }
    \item{nrCores}{
        Number of cores to use when \code{useMulticore = TRUE}.  Default value
        \code{NULL} uses all detected cores.
    }
}
\value{
    The function produces PNG files in the working directory, visualising
    density of patterns occurrence in the set of ordered input sequences. One
    file/plot per specified pattern is created.
}
\references{
    Haberle \emph{et al.} (2014) Two independent transcription initiation codes
    overlap on vertebrate core promoters, \emph{Nature} \bold{507}:381-385.
}
\author{Vanja Haberle}

\seealso{
    \code{\link{getPatternOccurrenceList}}
    \cr
    \code{\link{plotMotifDensityMap}}
}
\examples{
library(GenomicRanges)
load(system.file("data", "zebrafishPromoters.RData", package="seqPattern"))

promoterWidth <- elementMetadata(zebrafishPromoters)$interquantileWidth

# dinucleotide patterns
plotPatternDensityMap(regionsSeq = zebrafishPromoters, patterns = c("TA", "GC"),
            seqOrder = order(promoterWidth), flankUp = 400, flankDown = 600, 
            color = "blue")

# motif consensus sequence
plotPatternDensityMap(regionsSeq = zebrafishPromoters, patterns = "TATAWAWR",
            seqOrder = order(promoterWidth), flankUp = 400, flankDown = 600,
            color = "cyan")
}
