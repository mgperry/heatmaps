% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotHeatmapList.R
\name{plotHeatmapList}
\alias{plotHeatmapList}
\title{Plot a list of heatmaps}
\usage{
plotHeatmapList(heatmap_list, groups = 1:length(heatmap_list),
  options = heatmapOptions(), ...)
}
\arguments{
\item{heatmap_list}{A list of Heatmaps}

\item{groups}{Optionally group heatmaps together}

\item{options}{Heatmap options}

\item{...}{Additional options

This function takes a list of one or more heatmaps and plots
them to a single image tiled horizontally.

The "groups" argument specifies heatmaps to be grouped together and
plotted using the same display parameters and a unified scale. plotHeatmapList
will try to guess the best scale, either starting or finishing at zero, or symetrical
aronud zero - if this is not the desired behaviour, make sure the scales are
identical before the heatmaps are passed to the function.

Options are specified as for plotHeatmap, but can be specified per group by passing
a list of options instead of a single vector. Note the difference between a length-2
character vector, c("Reds", "Blues"), and a list contatining two length-1 character
vectors: list("Reds", "Blues").


These are generally large, complex plots, so it can better to plot
straight to a file. PNG is preferred since pdf files generated can be
if the images are not downsized. The default settings are designed for plots
of about 10cm x 20cm per heatmap, but all of the relevant settigns can be tweaked
using the options. For display-quality images, it helps to increase the resolution
at to at least 150ppi, double the default of 72ppi on most systems.}
}
\value{
invisible(0)
}
\description{
Plot a list of heatmaps
}
\examples{
data(HeatmapExamples)
plotHeatmapList(list(hm, hm2), groups=c(1,2), color=list("Reds", "Blues"))
}
\seealso{
plotHeatmap heatmapOptions plot_legend
}
