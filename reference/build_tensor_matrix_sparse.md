# Build Sparse Tensor Design Matrix

Inner function to construct the tensor product design matrix as a sparse
[`SparseM::matrix.csr`](https://rdrr.io/pkg/SparseM/man/SparseM.ontology.html),
exploiting the local support of the B-spline basis (each row has only a
handful of non-zero entries). The column layout is identical to
`cbind(1, build_tensor_matrix(...))`, i.e. it includes the leading
global intercept column, so the result can be passed directly to
[`quantreg::rq.fit.sfn`](https://rdrr.io/pkg/quantreg/man/rq.fit.sfn.html).

## Usage

``` r
build_tensor_matrix_sparse(scores, Y.mask, spline.basis)
```

## Arguments

- scores:

  A matrix of functional Principal Component Analysis (fPCA) scores for
  each subject.

- Y.mask:

  A logical mask matrix of dimensions (N x T) indicating which time
  points are observed (`TRUE`) or missing (`FALSE`) for each subject.

- spline.basis:

  The evaluation matrix of the B-spline basis functions at the time
  points.

## Value

A `matrix.csr` sparse matrix containing the intercept column, the
intercept bases and the scaled spline bases for all valid observations.
