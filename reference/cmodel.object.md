# Developer documentation: censoring model object

A list giving information about censored states, their labels in the
data and what true states they represent.

## Value

- ncens:

  The number of distinct values used for censored observations in the
  `state` data supplied to
  [`msm`](https://chjackson.github.io/msm/reference/msm.md).

- censor:

  A vector of length `ncens`, giving the labels used for censored states
  in the data.

- states:

  A vector obtained by
  [`unlist()`](https://rdrr.io/r/base/unlist.html)ing a list with
  `ncens` elements, each giving the set of true states that an
  observation with this label could be.

- index:

  Index into `states` for the first state corresponding to each
  `censor`, plus an extra `length(states)+1`.

## See also

[`msm.object`](https://chjackson.github.io/msm/reference/msm.object.md).
