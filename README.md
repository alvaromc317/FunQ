
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `FunQuantPCA` package

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Package
Version](https://img.shields.io/badge/version-0.1.1-blue.svg)](https://cran.r-project.org/package=yourPackageName)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-brightgreen.svg)](https://www.tidyverse.org/lifecycle/)
<!-- badges: end -->

`FunQuantPCA` is an R package that solves the functional quantile
principal component analysis (fqpca) methodology. FQPCA extends the
concept of functional principal component analysis (FPCA) to the
quantile regression framework. The goal of many methods in FDA is to
recover the curve-specific mean by leveraging information across time
points, subjects, or both. Our goal is broader: we seek to describe the
full curve and time-specific probability distribution that underlies
individual measurements. Although only one draw from the curve-and
time-specific distribution of each individual is available, we will
leverage information across time points and curves to estimate smooth,
curve-specific quantile functions. This approach allows a richer
understanding of functional data than considering only the expected
value, and may be particularly useful when distributions are skewed,
vary across subjects or present outliers.

## Installation

You can install the development version of FunQuantPCA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alvaromc317/FunQuantPCA")
```

## Main characteristics

The key methodology in the `FunQuantPCA` package is the `fqpca`
function, that implements the Functional Quantile Principal COmponent
Analysis methodology.

- This function can receive the functional data as a `matrix`,
  `dataframe`, `tibble` or `tf` (from the `tidyfun` package) objects.
- It can deal with irregular time grids, which means that can deal with
  missing data.
- Can control the level of smoothness of the results based on two
  approaches:
  - Based on the degrees of freedom of the spline basis reconstruction,
    using the parameter `splines.df` This is our preferred approach.
  - Based on the inclusion of a second derivative penalty on the splines
    coefficients, changing the `method` parameter to use a valid `CVXR`
    solver (for example, setting `method='SCS'`) and then selecting the
    desired hyper-parameter value (for example `alpha.ridge=1e-7`). This
    approach is experimental and is prone to show computational issues
    for large values of the hyper-parameter.

The `FunQuantPCA` also implements functions to perform cross validation
on either the `splines.df` parameter (`cross_validation_df`) or the
`alpha.ridge` parameter (`cross_validation_alpha`). These cross
validation functions consider the quantile error as the reference
prediction error. This error is available using the function
`quantile_error`.

## Example 1: setting the basics

`fqpca` is mainly designed to deal with functional data. The following
example generates a fake dataset with 200 observations taken every 10
minutes during one day. This defines a data matrix with 200 rows and 144
columns following the formula:

$$x_i = \lambda_1(sin(t)+sin(0.5t))+\varepsilon_i$$ where

- $\lambda_1\sim N(0,0.4)$
- $\varepsilon_i\sim\chi^2(3)$

``` r
set.seed(5)

n = 200
t = 144
time.points = seq(0, 2*pi, length.out=t)
Y = matrix(rep(sin(time.points) + sin(0.5*time.points), n), byrow=TRUE, nrow=n)
Y = Y + matrix(rnorm(n*t, 0, 0.4), nrow=n) + rchisq(n, 3)
matplot(t(Y[1:20,]), type="l")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

The above plot visualizes a subset of the data generated this way. Since
the fqpca methodology can deal with sparse and irregular time
measurements, we will include 50% of missing observations in the data
matrix.

``` r
Y[sample(n*t, as.integer(0.50*n*t))] = NA
```

Now, we apply the `fqpca` methodology on this dataset and obtain the
decomposition of the data in terms of the median (`quantile.value=0.5`),
which is a robust alternative to the mean based predictions of
traditional FPCA.

``` r
library(FunQuantPCA)

Y.train = Y[1:150,]
Y.test = Y[151:n,]

results = fqpca(Y=Y.train, npc=2, quantile.value=0.5)

loadings = results$loadings
scores = results$scores

# Recover x_train based on decomposition
Y.train.estimated = scores %*% t(loadings)
```

Finally, given a new set of observations, it is possible to decompose
the new observations using the loadings already computed.

``` r
scores.test = predict(results, newdata=Y.test)
Y.test.estimated = scores.test %*% t(loadings)
```

You can plot the computed loadings on a somewhat not very pretty plot,
but still useful plot

``` r
plot(results)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

And you can also compute the quantile error between the curve
reconstruction and the true data, which is the metric we recommend to
use as prediction error metric.

``` r
quantile_error(Y=Y.train, Y.pred=Y.train.estimated, quantile.value=0.5)
#> [1] 0.1575061
```

## Example 2: cross validating

The `FunQuantPCA` package implements functions that allow to perform
cross validation based on both the `splines.df` or the `alpha.ridge`
criterias. Let’s see an example using the well known weather dataset.

``` r
# We use the tidy structure of the tidyfun package to deal with the functional data
devtools::install_github("tidyfun/tidyfun")
```

``` r
data = t(fda::CanadianWeather$dailyAv[,,1])
tf_data = tf::tfd(data, arg = 1:365)

plot(tf_data, col='black')
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
splines.df.grid = c(5, 10, 15)
cv_result = cross_validation_df(tf_data, splines.df.grid=splines.df.grid, n.folds=2)
#> Degrees of freedom: 5 ---------------------
#> Fold: 1
#> Fold: 2
#> Degrees of freedom: 5 .Execution completed in: 0.55 seconds
#> Degrees of freedom: 10 ---------------------
#> Fold: 1
#> Fold: 2
#> Degrees of freedom: 10 .Execution completed in: 0.9 seconds
#> Degrees of freedom: 15 ---------------------
#> Fold: 1
#> Fold: 2
#> Degrees of freedom: 15 .Execution completed in: 1.15 seconds

cv_result$error.matrix
#>           [,1]      [,2]
#> [1,] 0.5597133 0.5693647
#> [2,] 0.5267699 0.5433952
#> [3,] 0.5244181 0.5390324
```

The dimensions of the error matrix are (length(splines.df.grid),
n.folds). We can find the optimal number of degrees of freedom by
taking, for example, the mean of each row and picking the minimum.

``` r
optimal_df = which.min(rowMeans(cv_result$error.matrix))

paste0('Optimal number of degrees of freedom: ', splines.df.grid[optimal_df])
#> [1] "Optimal number of degrees of freedom: 15"
```

Now we can build the final model using the optimal number of degrees of
freedom, and check the number of components based on the percentage of
explained variability.

``` r
results = fqpca(Y=tf_data, npc=10, quantile.value=0.5, splines.df=15, seed=5)

cumsum(results$pve)
#>  [1] 0.8881699 0.9733824 0.9922474 0.9974198 0.9984257 0.9992301 0.9996471
#>  [8] 0.9998388 0.9999582 1.0000000
```

This shows that with 2 components we are able to explain 97% of the
variability in the data.

The package also includes a function that allows to perform cross
validation on the hyper-parameter controlling the effect of a second
derivative penalty on the splines. Be aware that this smoothness
controlling process is experimental and may be subject to computation
issues.

``` r
cv_result = cross_validation_alpha(tf_data, alpha.grid=c(0, 1e-10, 1e-5), n.folds=2)
#> alpha=0 ---------------------
#> Fold: 1
#> Fold: 2
#> Alpha 0 execution completed in: 33.83 seconds
#> alpha=1e-10 ---------------------
#> Fold: 1
#> Fold: 2
#> Alpha 1e-10 execution completed in: 46.52 seconds
#> alpha=1e-05 ---------------------
#> Fold: 1
#> Fold: 2
#> Alpha 1e-05 execution completed in: 33.19 seconds

cv_result$error.matrix
#>           [,1]      [,2]
#> [1,] 0.5337405 0.5356558
#> [2,] 0.5399521 0.5355772
#> [3,] 0.6243240 0.6118602
```

The dimensions of the error matrix are (length(alpha.grid), n.folds)
