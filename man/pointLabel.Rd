% Generated by roxygen2 (4.0.2): do not edit by hand
\name{pointLabel}
\alias{pointLabel}
\title{Label placement for points to avoid overlaps}
\usage{
pointLabel(x, y = NULL, labels = seq(along = x), cex = 1,
  method = c("SANN", "GA"), allowSmallOverlap = FALSE, trace = FALSE,
  doPlot = TRUE, ...)
}
\arguments{
\item{x,y,labels,cex,method,allowSmallOverlap,trace,doPlot}{see maptools::pointLabel}

\item{...}{see maptools::pointLabel}
}
\description{
pointLabel was taken from package maptools. maptools is not imported
or mentioned in DESCRIPTION to reduce dependencies as maptools
requires sp and gpclib. Below is the exact maptools::pointLabel code.
See pointLabel in maptools for a description of the arguments
}
\keyword{internal}

