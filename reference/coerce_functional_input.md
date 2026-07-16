# Coerce functional input data

Coerces the supported functional-data containers (matrix, tf object, or
data.frame with a tf column) into an unnamed matrix, optionally
resolving a grouping structure.

## Usage

``` r
coerce_functional_input(data, colname = NULL, group = NULL)
```

## Arguments

- data:

  An \\(N \times T)\\ matrix, a tf object from the tidyfun package or a
  data.frame containing the functional data as a tf column.

- colname:

  The name of the column containing the functional data. Use only if
  data is a dataframe and colname is a column in the dataframe.

- group:

  Optional grouping structure of hierarchical data. If data is a
  dataframe this can be a name pointing to a column in the dataframe. If
  data is any other object it must be an array with length equal to the
  number of rows of data. Set to NULL when no grouping is used.

## Value

A list with elements `Y` (unnamed matrix of functional data) and `group`
(the processed group vector, or NULL).
