% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mesh2projector.R
\name{mesh2projector}
\alias{mesh2projector}
\title{Illustrative code to build the projector matrix for SPDE models.}
\usage{
mesh2projector(
  mesh,
  loc = NULL,
  lattice = NULL,
  xlim = NULL,
  ylim = NULL,
  dims = c(100, 100)
)
}
\arguments{
\item{mesh}{a 2d mesh object.}

\item{loc}{a two columns matrix with the locations to project for.}

\item{lattice}{Unused; feature not supported by this illustration.}

\item{xlim, ylim}{vector with the boundary limits.}

\item{dims}{the number of subdivisions over each boundary limits.}
}
\value{
the projector matrix as a list with sparse matrix object at \code{x$proj$A}..
}
\description{
Creates a projector matrix object.
}
\section{Warning}{

This is just for illustration purpose and one should consider the
efficient functions available in the fmesher package,
e.g. \code{\link[fmesher:fm_evaluate]{fmesher::fm_evaluator()}} and \code{\link[fmesher:fm_basis]{fmesher::fm_basis()}}.
}

