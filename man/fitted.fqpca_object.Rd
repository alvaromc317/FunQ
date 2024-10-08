% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fqpca.R
\name{fitted.fqpca_object}
\alias{fitted.fqpca_object}
\title{Fit Yhat}
\usage{
\method{fitted}{fqpca_object}(object, pve = 0.95, ...)
}
\arguments{
\item{object}{An object output of the fqpca function.}

\item{pve}{If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The normalized matrix of scores.
}
\description{
S3 method for class 'fqpca_object'. Given an fqpca_object model, estimates Yhat for different pve values.
}
\examples{
# Generate fake dataset with 150 observations and 144 time points

Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)

# Add missing observations
Y[sample(150*144, as.integer(0.2*150*144))] <- NA

results <- fqpca(data = Y[1:100,], npc = 2, quantile.value = 0.5)

Yhat <- fitted(object = results, pve=0.99)
}
