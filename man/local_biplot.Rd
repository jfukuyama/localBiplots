% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdessence-functions.R
\name{local_biplot}
\alias{local_biplot}
\title{Create local biplot axes}
\usage{
local_biplot(
  X,
  dist,
  dist_deriv = NULL,
  k = 2,
  samples = 1:nrow(X),
  new_points = list()
)
}
\arguments{
\item{X}{A data matrix, samples as rows.}

\item{dist}{Either a string describing one of the supported
distances or a function that takes a matrix and returns the
distances between the rows of the matrix.}

\item{dist_deriv}{Either NULL (if dist is a string describing one
of the supported distances) or a function that takes two
vectors and computes \eqn{\frac{\partial}{\partial y_j}d(x,y)}.}

\item{k}{The number of embedding dimensions.}

\item{samples}{The samples to compute local biplot axes
at. Defaults to all of the original samples.}

\item{new_points}{New points (not one of the original samples) to
compute local biplot axes at.}
}
\value{
A data frame. Each row describes one LB axis for one
    sample. Columns labeled 'Embedding' give the embedding of the
    sample in MDS space, columns labeled 'Axis' give the LB axis
    for a given sample and variable. Columns labeled 'variable' and
    'sample' give the variable and sample for the LB axis values.
}
\description{
Create local biplot axes
}
