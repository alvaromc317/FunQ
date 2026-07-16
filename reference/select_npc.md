# Select the number of components

Computes the required number of components to achieve the desired pve
value

## Usage

``` r
select_npc(scores, pve, npc = NULL)
```

## Arguments

- scores:

  matrix of scores

- pve:

  Percentage of explained variability, a number between 0 and 1. Set to
  NULL to use all components.

- npc:

  Number of components. If provided (non NULL), it supersedes `pve`.
  Must be an integer between 0 and the number of components in `scores`.
