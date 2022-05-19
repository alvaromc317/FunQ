
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `fqpca` package

<!-- badges: start -->
<!-- badges: end -->

`fqpca` is an R package that solves the functional quantile principal
component analysis (fqpca) methodology. This is a dimensionality
reduction technique that extends the concept of functional principal
component analysis (FPCA) to the quantile regression framework. FPCA is
commonly used to identify and describe variation in samples of curves,
but is only able to provide mean based estimations which are known to be
affected by outliers and skewness. On the other hand, FQPCA is able to
capture shifts on the scale of the data affecting the quantiles, and is
also a robust methodology suitable for dealing with outliers,
heteroscedastic data or skewed data. The model is estimated using
penalized M splines, and can deal with sparse and irregular time
measurements.

## Installation

You can install the development version of fqpca from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alvaromc317/fqpca")
```

## Example

`fqpca` is mainly designed to deal with functional data. The following
example generates a fake dataset with 200 observations taken every 10
minutes during one day. This defines a data matrix with 200 rows and 144
columns following the formula:

![x_i = \\lambda_1(sin(t)+sin(0.5t))+\\varepsilon_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i%20%3D%20%5Clambda_1%28sin%28t%29%2Bsin%280.5t%29%29%2B%5Cvarepsilon_i "x_i = \lambda_1(sin(t)+sin(0.5t))+\varepsilon_i")

where \*
![\\lambda_1\\sim N(0,0.4)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_1%5Csim%20N%280%2C0.4%29 "\lambda_1\sim N(0,0.4)")
\*
![\\varepsilon_i\\sim\\chi^2(3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cvarepsilon_i%5Csim%5Cchi%5E2%283%29 "\varepsilon_i\sim\chi^2(3)")

``` r
n = 200
t = 144
time.points = seq(0, 2*pi, length.out=t)
x = matrix(rep(sin(time.points) + sin(0.5*time.points), n), byrow=TRUE, nrow=n)
x = x + matrix(rnorm(n*t, 0, 0.4), nrow=n) + rchisq(n, 3)
matplot(t(x[1:20,]), type="l")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

The above plot visualizes a subset of the data generated this way. Since
the main methodology in the package can deal with sparse and irregular
time measurements, we will include 50% of missing observations in the
data matrix.

``` r
x[sample(n*t, as.integer(0.5*n*t))] = NA
```

Now, we apply the `fqpca` methodology on this dataset and obtain the
decomposition of the data in terms of the median (`quantile.value=0.5`),
which is a robust alternative to the mean based predictions of
traditional FPCA.

``` r
library(fqpca)
x.train = x[1:150,]
x.test = x[151:n,]
results = fqpca(x=x.train, n.components=2, quantile.value=0.5)

loadings = results$loadings
scores = results$scores

# Recover x.train based on decomposition
x.train.estimated = scores %*% t(loadings)
```

Finally, given a new set of observations, it is possible to decompose
the new observations using the loadings already computed.

``` r
scores.test = predict_scores_fqpca(results, x=x.test)
x.test.estimated = scores.test %*% t(loadings)
```
