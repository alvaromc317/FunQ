# Package index

## Core Methodologies

Functions implementing main quantile PCA algorithms

- [`fqpca()`](https://alvaromc317.github.io/FunQ/reference/fqpca.md) :
  FQPCA (Functional Quantile Principal Component Analysis)
- [`mfqpca()`](https://alvaromc317.github.io/FunQ/reference/mfqpca.md) :
  MFQPCA (Multilevel Functional Quantile Principal Component Analysis)

## Parameter Cross-Validation

Functions to perform k-fold cross-validation on hyperparameters

- [`fqpca_cv_df()`](https://alvaromc317.github.io/FunQ/reference/fqpca_cv_df.md)
  : Cross validation of degrees of freedom
- [`fqpca_cv_lambda()`](https://alvaromc317.github.io/FunQ/reference/fqpca_cv_lambda.md)
  : Cross validation of lambda.ridge
- [`mfqpca_cv_df()`](https://alvaromc317.github.io/FunQ/reference/mfqpca_cv_df.md)
  : Cross validation of degrees of freedom

## S3 Methods & Plotting

Methods to plot, fit, and predict from fitted model objects

- [`fitted(`*`<fqpca_object>`*`)`](https://alvaromc317.github.io/FunQ/reference/fitted.fqpca_object.md)
  : Fit Yhat
- [`fitted(`*`<mfqpca_object>`*`)`](https://alvaromc317.github.io/FunQ/reference/fitted.mfqpca_object.md)
  : Fit Yhat
- [`plot(`*`<fqpca_object>`*`)`](https://alvaromc317.github.io/FunQ/reference/plot.fqpca_object.md)
  : Plot fqpca loading functions with ggplot2
- [`plot(`*`<mfqpca_object>`*`)`](https://alvaromc317.github.io/FunQ/reference/plot.mfqpca_object.md)
  : Plot fqpca loading functions with ggplot2
- [`predict(`*`<fqpca_object>`*`)`](https://alvaromc317.github.io/FunQ/reference/predict.fqpca_object.md)
  : Predict fqpca scores
- [`predict(`*`<mfqpca_object>`*`)`](https://alvaromc317.github.io/FunQ/reference/predict.mfqpca_object.md)
  : Predict mfqpca scores

## Utilities & Helper Metrics

Internal calculations and error metrics

- [`create_folds()`](https://alvaromc317.github.io/FunQ/reference/create_folds.md)
  : Split a given matrix Y into k-folds
- [`train_test_split()`](https://alvaromc317.github.io/FunQ/reference/train_test_split.md)
  : Split a given matrix Y into train / test
- [`quantile_error()`](https://alvaromc317.github.io/FunQ/reference/quantile_error.md)
  : Quantile error computation
- [`proportion_under_quantile()`](https://alvaromc317.github.io/FunQ/reference/proportion_under_quantile.md)
  : Proportion of points under quantile

## Internal Utilities

Internal helper functions not exported by the package

- [`build_tensor_matrix()`](https://alvaromc317.github.io/FunQ/reference/build_tensor_matrix.md)
  : Build Tensor Design Matrix
- [`build_tensor_matrix_sparse()`](https://alvaromc317.github.io/FunQ/reference/build_tensor_matrix_sparse.md)
  : Build Sparse Tensor Design Matrix
- [`check_flag()`](https://alvaromc317.github.io/FunQ/reference/check_flag.md)
  : Check that a parameter is a single Boolean
- [`check_fqpca_params()`](https://alvaromc317.github.io/FunQ/reference/check_fqpca_params.md)
  : Check fqpca input parameters
- [`check_iters()`](https://alvaromc317.github.io/FunQ/reference/check_iters.md)
  : Check that a maximum-iterations parameter is an integer \>= 2
- [`check_loss()`](https://alvaromc317.github.io/FunQ/reference/check_loss.md)
  : Quantile check (pinball) loss
- [`check_nonnegative_number()`](https://alvaromc317.github.io/FunQ/reference/check_nonnegative_number.md)
  : Check that a parameter is a non-negative number
- [`check_positive_float()`](https://alvaromc317.github.io/FunQ/reference/check_positive_float.md)
  : Check that a parameter is a positive float
- [`check_positive_integer()`](https://alvaromc317.github.io/FunQ/reference/check_positive_integer.md)
  : Check that a parameter is a positive integer
- [`check_quantile_value()`](https://alvaromc317.github.io/FunQ/reference/check_quantile_value.md)
  : Check that the quantile value lies in (0, 1)
- [`check_regressors()`](https://alvaromc317.github.io/FunQ/reference/check_regressors.md)
  : Checks input regressors
- [`check_seed()`](https://alvaromc317.github.io/FunQ/reference/check_seed.md)
  : Check that the seed is a positive integer or NULL
- [`check_splines_method()`](https://alvaromc317.github.io/FunQ/reference/check_splines_method.md)
  : Check that the splines method is supported
- [`coerce_functional_input()`](https://alvaromc317.github.io/FunQ/reference/coerce_functional_input.md)
  : Coerce functional input data
- [`compute_loadings()`](https://alvaromc317.github.io/FunQ/reference/compute_loadings.md)
  : Compute loadings (aka principal components)
- [`compute_objective_value()`](https://alvaromc317.github.io/FunQ/reference/compute_objective_value.md)
  : Compute objective value function.
- [`compute_scores()`](https://alvaromc317.github.io/FunQ/reference/compute_scores.md)
  : Compute scores
- [`cv_framework()`](https://alvaromc317.github.io/FunQ/reference/cv_framework.md)
  : Generic cross-validation framework
- [`explained_variance_ratio()`](https://alvaromc317.github.io/FunQ/reference/explained_variance_ratio.md)
  : Explained variance ratio of the scores
- [`extract_coefficients()`](https://alvaromc317.github.io/FunQ/reference/extract_coefficients.md)
  : Extract the coefficient vector from a quantile regression model
  object
- [`fit_tensor_quantile_regression()`](https://alvaromc317.github.io/FunQ/reference/fit_tensor_quantile_regression.md)
  : Fit the tensor-design quantile regression
- [`kfold_cv_points()`](https://alvaromc317.github.io/FunQ/reference/kfold_cv_points.md)
  : Split a given matrix Y into k-folds
- [`kfold_cv_rows()`](https://alvaromc317.github.io/FunQ/reference/kfold_cv_rows.md)
  : Split a given matrix Y into k-folds
- [`mfqpca_check_params()`](https://alvaromc317.github.io/FunQ/reference/mfqpca_check_params.md)
  : Check mfqpca input parameters
- [`mfqpca_compute_scores_between()`](https://alvaromc317.github.io/FunQ/reference/mfqpca_compute_scores_between.md)
  : Compute scores between
- [`reconstruct_scores_loadings()`](https://alvaromc317.github.io/FunQ/reference/reconstruct_scores_loadings.md)
  : Truncated reconstruction from scores and loadings
- [`rotate_factors()`](https://alvaromc317.github.io/FunQ/reference/rotate_factors.md)
  : Rotation of loadings and scores
- [`select_npc()`](https://alvaromc317.github.io/FunQ/reference/select_npc.md)
  : Select the number of components
- [`train_test_split_points()`](https://alvaromc317.github.io/FunQ/reference/train_test_split_points.md)
  : Split points of a given matrix Y into train / test
- [`train_test_split_rows()`](https://alvaromc317.github.io/FunQ/reference/train_test_split_rows.md)
  : Split rows of a given matrix Y into train / test
