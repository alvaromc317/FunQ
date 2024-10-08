% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fqpca-utils.R
\name{train_test_split}
\alias{train_test_split}
\title{Split a given matrix Y into train / test}
\usage{
train_test_split(
  Y,
  criteria = "points",
  train.pct = NULL,
  train.size = NULL,
  seed = NULL
)
}
\arguments{
\item{Y}{\eqn{(N \times P)} matrix of predictive variables.}

\item{criteria}{Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.}

\item{train.pct}{Float number indicating the \% of rows used for training. This takes precedence over \code{train.size}.}

\item{train.size}{Integer number indicating number of rows used for training. \code{train.size} is superseded by \code{train.pct}.}

\item{seed}{Seed for the random generator number.}
}
\value{
A list containing a matrix Y.train and a matrix Y.test
}
\description{
Splits a given matrix into a train / test split based on two possible criteria: either based on the total number of rows or the total number of data points.
}
\examples{
# Generate a small matrix
Y <- matrix(rnorm(50), nrow = 10)

# Divide based on full rows
train_test_rows <- train_test_split(Y, criteria = 'rows', train.pct = 0.5)
Y.train <- train_test_rows$Y.train
Y.test <- train_test_rows$Y.test

train_test_points <- train_test_split(Y, criteria = 'points', train.pct = 0.5)
Y.train <- train_test_points$Y.train
Y.test <- train_test_points$Y.test
}
