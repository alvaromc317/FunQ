
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
function, that implements the Functional Quantile Principal Component
Analysis methodology.

- This function can receive the functional data as an $(N\times T)$
  `matrix` (through parameter `Y`) or as a dataframe containing a column
  with a functional `tf` vector (through parameters `data` and
  `colname`).
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

Y = matrix(rep(sin(time.points) + sin(0.5*time.points), n), byrow=TRUE, nrow=n)
Y = Y + matrix(rnorm(n*t, 0, 0.4), nrow=n) + rchisq(n, 3)
matplot(t(Y[1:20,]), type="l")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

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
Y.train.estimated = fitted(results, pve = 0.95)
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
plot(results, pve=0.95)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

And you can also compute the quantile error between the curve
reconstruction and the true data, which is the metric we recommend to
use as prediction error metric.

``` r
quantile_error(Y=Y.train, Y.pred=Y.train.estimated, quantile.value=0.5)
#> [1] 0.1597887
```

## Example 2: cross validating the Canadian Weather dataset

The `FunQuantPCA` package implements functions that allow to perform
cross validation based on both the `splines.df` or the `alpha.ridge`
criterias. Let’s see an example using the well known weather dataset.

``` r
# We use the tidy structure of the tidyfun package to deal with the functional data
devtools::install_github("tidyfun/tidyfun")
```

``` r

library(fda)
#> Loading required package: splines
#> Loading required package: fds
#> Loading required package: rainbow
#> Loading required package: MASS
#> Loading required package: pcaPP
#> Loading required package: RCurl
#> Loading required package: deSolve
#> 
#> Attaching package: 'fda'
#> The following object is masked from 'package:graphics':
#> 
#>     matplot
library(tidyverse)
#> Warning: package 'ggplot2' was built under R version 4.3.3
#> Warning: package 'readr' was built under R version 4.3.3
#> Warning: package 'stringr' was built under R version 4.3.3
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.3     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
#> ✔ purrr     1.0.2
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ tidyr::complete() masks RCurl::complete()
#> ✖ dplyr::filter()   masks stats::filter()
#> ✖ dplyr::lag()      masks stats::lag()
#> ✖ dplyr::select()   masks MASS::select()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(tidyfun)
#> Loading required package: tf
#> 
#> Attaching package: 'tf'
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     sd, var
#> 
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2

data = tibble(temperature = tf::tfd(t(fda::CanadianWeather$dailyAv[,,1]), arg = 1:365),
              province = CanadianWeather$province)

head(data)
#> # A tibble: 6 × 2
#>                      temperature province     
#>                        <tfd_reg> <chr>        
#> 1 [1]: (1,-4);(2,-3);(3,-3); ... Newfoundland 
#> 2 [2]: (1,-4);(2,-4);(3,-5); ... Nova Scotia  
#> 3 [3]: (1,-4);(2,-4);(3,-5); ... Nova Scotia  
#> 4 [4]: (1,-1);(2,-2);(3,-2); ... Nova Scotia  
#> 5 [5]: (1,-6);(2,-6);(3,-7); ... Ontario      
#> 6 [6]: (1,-8);(2,-8);(3,-9); ... New Brunswick
```

``` r

data %>% 
  ggplot(aes(y=temperature, color=province)) + 
  geom_spaghetti() + 
  theme_light()
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

The `fqpca` function admits the data as a `tf` functional vector.

``` r
splines.df.grid = c(5, 10, 15, 20)
cv_result = cross_validation_df(data=data, colname='temperature', splines.df.grid=splines.df.grid, n.folds=3)
#> Degrees of freedom: 5 ---------------------
#> Fold: 1
#> Fold: 2
#> Fold: 3
#> Degrees of freedom: 5. Execution completed in: 0.91 seconds.
#> Degrees of freedom: 10 ---------------------
#> Fold: 1
#> Fold: 2
#> Fold: 3
#> Degrees of freedom: 10. Execution completed in: 0.71 seconds.
#> Degrees of freedom: 15 ---------------------
#> Fold: 1
#> Fold: 2
#> Fold: 3
#> Degrees of freedom: 15. Execution completed in: 1.4 seconds.
#> Degrees of freedom: 20 ---------------------
#> Fold: 1
#> Fold: 2
#> Fold: 3
#> Degrees of freedom: 20. Execution completed in: 2.74 seconds.
cv_result$error.matrix
#>           [,1]      [,2]      [,3]
#> [1,] 0.5600411 0.5580592 0.5632948
#> [2,] 0.5346743 0.5288065 0.5340524
#> [3,] 0.5315561 0.5240813 0.5298265
#> [4,] 0.5343425 0.5237205 0.5289621
```

The dimensions of the error matrix are
`(length(splines.df.grid), n.folds)`. We can find the optimal number of
degrees of freedom by taking the mean of each row and picking the
minimum.

``` r
optimal_df = which.min(rowMeans(cv_result$error.matrix))

paste0('Optimal number of degrees of freedom: ', splines.df.grid[optimal_df])
#> [1] "Optimal number of degrees of freedom: 15"
```

Now we can build the final model using the optimal number of degrees of
freedom, and check the number of components based on the percentage of
explained variability.

``` r
results = fqpca(Y=data$temperature, npc=10, quantile.value=0.5, splines.df=15, seed=5)

cumsum(results$pve)
#>  [1] 0.8881699 0.9733824 0.9922474 0.9974198 0.9984257 0.9992301 0.9996471
#>  [8] 0.9998388 0.9999582 1.0000000
```

This shows that with 2 components we are able to explain 97% of the
variability in the data. Lets see these components.

``` r

plot(results, pve = 0.95)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

The package also includes a function that allows to perform cross
validation on the hyper-parameter controlling the effect of a second
derivative penalty on the splines. Be aware that this smoothness
controlling process is experimental and may be subject to computation
issues.

``` r
cv_result = cross_validation_alpha(data$temperature, alpha.grid=c(0, 1e-10, 1e-5), n.folds=2)
#> alpha=0 ---------------------
#> Fold: 1
#> Fold: 2
#> alpha: 0. Execution completed in: 15.27 seconds.
#> alpha=1e-10 ---------------------
#> Fold: 1
#> Fold: 2
#> alpha: 1e-10. Execution completed in: 19.07 seconds.
#> alpha=1e-05 ---------------------
#> Fold: 1
#> Warning in fqpca(Y = Y.train, npc = npc, quantile.value = quantile.value, :
#> Algorithm reached maximum number of iterations without convergence: 20
#> iterations
#> Fold: 2
#> Warning in fqpca(Y = Y.train, npc = npc, quantile.value = quantile.value, :
#> Algorithm reached maximum number of iterations without convergence: 20
#> iterations
#> alpha: 1e-05. Execution completed in: 72.38 seconds.

cv_result$error.matrix
#>           [,1]      [,2]
#> [1,] 0.5340146 0.5355163
#> [2,] 0.5341246 0.5359910
#> [3,] 1.0199191 0.8868765
```

The dimensions of the error matrix are (length(alpha.grid), n.folds)
