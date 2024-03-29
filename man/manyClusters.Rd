% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manyClusters.R
\name{manyClusters}
\alias{manyClusters}
\title{Creates a pdf with different dendrograms (hierarchical clusters) combining a distance and a linkage method}
\usage{
manyClusters(
  estimates,
  resultsDir = NULL,
  fileName = "manyClusters",
  fmtPlot = "pdf",
  title = NULL,
  conditions = NULL,
  colors = NULL
)
}
\arguments{
\item{estimates}{Numeric matrix with the intensity values to plot. The matrix should have samples as columns and genes/transcripts/exons as rows}

\item{resultsDir}{Output directory. Default = NULL}

\item{fileName}{name of the output file, without extension. Default = "manyClusters"}

\item{fmtPlot}{Format for the image file. If none specified images will be generated in R session. Default = "pdf"}

\item{title}{Title for each plot. Default = NULL}

\item{conditions}{Vector with the different conditions. Default = NULL}

\item{colors}{Vector with the colors assigned to each condition. Default = NULL}
}
\value{
a pdf with the hierarchical clusters
}
\description{
Creates a pdf with different dendrograms (hierarchical clusters) combining a distance and a linkage method
}
