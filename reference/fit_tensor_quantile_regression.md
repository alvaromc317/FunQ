# Fit the tensor-design quantile regression

Inner function to fit the spline coefficients of the alternating
quantile algorithms: one large quantile regression of the vectorised
data on the tensor-product design built from the current scores and the
spline basis.

## Usage

``` r
fit_tensor_quantile_regression(
  Y.vector,
  Y.mask,
  scores,
  spline.basis,
  quantile.value,
  method,
  intercept = TRUE,
  penalized = FALSE,
  lambda.ridge = 0,
  return.model = FALSE
)
```

## Arguments

- Y.vector:

  vectorised version of Y, the \\(N \times T)\\ matrix of observed time
  instants.

- Y.mask:

  Mask matrix of the same dimensions as Y indicating which observations
  in Y are known.

- scores:

  Matrix of scores (and / or regressors) scaling the spline blocks.

- spline.basis:

  The spline basis matrix.

- quantile.value:

  The quantile considered.

- method:

  Method used in the resolution of the quantile regression model. It
  currently accepts the methods `c('conquer', 'quantreg')`.

- intercept:

  Boolean. If TRUE (default) the model contains a functional intercept;
  if FALSE the fit has no intercept of any kind (used by the mfqpca
  within level).

- penalized:

  Boolean indicating if the smoothness should be controlled using a
  second derivative penalty (conquer only).

- lambda.ridge:

  Hyper parameter controlling the penalization on the second derivative
  of the splines. Only used when `penalized=TRUE`.

- return.model:

  Boolean. If TRUE, returns a list with the fitted model object, the
  spline coefficient matrix, the flat coefficient vector and the dense
  tensor design matrix (needed for the fosqr variance estimation).

## Value

The matrix of spline coefficients (with `npc + 1` columns when
`intercept=TRUE`, `npc` otherwise), or the extended list when
`return.model=TRUE`.
