% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fqpca.R
\name{compute_loadings}
\alias{compute_loadings}
\title{Compute loadings (aka principal components)}
\usage{
compute_loadings(
  spline.basis,
  intercept.spline.basis,
  spline.coefficients,
  method
)
}
\arguments{
\item{spline.basis}{The spline basis matrix.}

\item{intercept.spline.basis}{The spline basis matrix for the conquer intercept.}

\item{spline.coefficients}{the matrix of spline coefficients.}

\item{method}{Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')} along with any available solver in \code{CVXR} package.}
}
\value{
The matrix of loadings
}
\description{
Inner function to compute the loadings of the fqpca methodology.
}
