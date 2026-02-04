# Explore the likelihood surface

Plot the log-likelihood surface with respect to two parameters.

## Usage

``` r
surface.msm(
  x,
  params = c(1, 2),
  np = 10,
  type = c("contour", "filled.contour", "persp", "image"),
  point = NULL,
  xrange = NULL,
  yrange = NULL,
  ...
)

# S3 method for class 'msm'
contour(x, ...)

# S3 method for class 'msm'
persp(x, ...)

# S3 method for class 'msm'
image(x, ...)
```

## Arguments

- x:

  Output from [`msm`](https://chjackson.github.io/msm/reference/msm.md),
  representing a fitted msm model.

- params:

  Integer vector with two elements, giving the indices of the parameters
  to vary. All other parameters will be fixed. Defaults to `c(1,2)`,
  representing the first two log transition intensities. See the
  `fixedpars` argument to `msm` for a definition of these indices.

- np:

  Number of grid points to use in each direction, by default 10. An
  `np x np` grid will be used to evaluate the likelihood surface. If 100
  likelihood function evaluations is slow, then reduce this.

- type:

  Character string specifying the type of plot to produce.

  |                    |                                                                                                                    |
  |--------------------|--------------------------------------------------------------------------------------------------------------------|
  | `"contour"`        | Contour plot, using the R function [`contour`](https://rdrr.io/r/graphics/contour.html).                           |
  | `"filled.contour"` | Solid-color contour plot, using the R function [`filled.contour`](https://rdrr.io/r/graphics/filled.contour.html). |
  | `"persp"`          | Perspective plot, using the R function [`persp`](https://rdrr.io/r/graphics/persp.html).                           |
  | `"image"`          | Grid color plot, using the R function [`image`](https://rdrr.io/r/graphics/image.html).                            |

- point:

  Vector of length `n`, where `n` is the number of parameters in the
  model, including the parameters that will be varied here. This
  specifies the point at which to fix the likelihood. By default, this
  is the maximum likelihood estimates stored in the fitted model `x`,
  `x$estimates`.

- xrange:

  Range to plot for the first varied parameter. Defaults to plus and
  minus two standard errors, obtained from the Hessian at the maximum
  likelihood estimate.

- yrange:

  Range to plot for the second varied parameter. Defaults to plus and
  minus two standard errors, obtained from the Hessian at the maximum
  likelihood estimate.

- ...:

  Further arguments to be passed to the plotting function.

## Details

Draws a contour or perspective plot. Useful for diagnosing
irregularities in the likelihood surface. If you want to use these plots
before running the maximum likelihood estimation, then just run `msm`
with all estimates fixed at their initial values.

`contour.msm` just calls surface.msm with `type = "contour"`.

`persp.msm` just calls surface.msm with `type = "persp"`.

`image.msm` just calls surface.msm with `type = "image"`.

As these three functions are methods of the generic functions `contour`,
`persp` and `image`, they can be invoked as `contour(x)`, `persp(x)` or
`image(x)`, where `x` is a fitted `msm` object.

## See also

[`msm`](https://chjackson.github.io/msm/reference/msm.md),
[`contour`](https://rdrr.io/r/graphics/contour.html),
[`filled.contour`](https://rdrr.io/r/graphics/filled.contour.html),
[`persp`](https://rdrr.io/r/graphics/persp.html),
[`image`](https://rdrr.io/r/graphics/image.html).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
