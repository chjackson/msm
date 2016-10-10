\name{surface.msm}
\alias{surface.msm}
\alias{persp.msm}
\alias{contour.msm}
\alias{image.msm}
\title{Explore the likelihood surface}
\description{
  Plot the log-likelihood surface with respect to two parameters. 
}
\usage{
surface.msm(x, params=c(1,2), np=10, type=c("contour","filled.contour","persp","image"),
                        point=NULL, xrange=NULL, yrange=NULL,...)
\method{contour}{msm}(x, ...)
\method{persp}{msm}(x, ...)
\method{image}{msm}(x, ...)
}
\arguments{
  \item{x}{Output from \code{\link{msm}}, representing a fitted msm model.}
  \item{params}{Integer vector with two elements, giving the indices of
    the parameters to vary. All other parameters
    will be fixed. Defaults to \code{c(1,2)}, representing the first two
  log transition intensities. See the \code{fixedpars} argument to
  \code{msm} for a definition of these indices.}
  \item{np}{Number of grid points to use in each direction, by default
    10.   An \code{np x np} grid will be used to evaluate the likelihood
    surface. If 100 likelihood function evaluations is slow, then reduce this. }
  \item{type}{Character string specifying the type of plot to produce.
    \tabular{ll}{
    \code{"contour"} \tab Contour plot, using the R function
    \code{\link{contour}}. \cr
    \code{"filled.contour"} \tab Solid-color contour plot, using the R function \code{\link{filled.contour}}. \cr
    \code{"persp"} \tab Perspective plot, using the R function \code{\link{persp}}. \cr
    \code{"image"} \tab Grid color plot, using the R function \code{\link{image}}. \cr
    }
  }
  \item{point}{Vector of length \code{n}, where \code{n} is the number
    of parameters in the model, including the parameters that will be
    varied here.  This specifies the point at which to fix the likelihood.  By default, this
    is the maximum likelihood estimates stored in the fitted model
    \code{x}, \code{x$estimates}.}
  \item{xrange}{Range to plot for the first varied parameter.  Defaults
    to plus and minus two standard errors, obtained from the Hessian at
    the maximum likelihood estimate.}
  \item{yrange}{Range to plot for the second varied parameter.  Defaults
    to plus and minus two standard errors, obtained from the Hessian at
    the maximum likelihood estimate.}
  \item{...}{Further arguments to be passed to the plotting function.}
}
\details{
  Draws a contour or perspective plot.  Useful for diagnosing
  irregularities in the likelihood surface.  If you want to use these
  plots before running the maximum likelihood estimation, then just run
  \code{msm} with all estimates fixed at their initial values. 

  \code{contour.msm} just calls surface.msm with \code{type = "contour"}.

  \code{persp.msm} just calls surface.msm with \code{type = "persp"}.

  \code{image.msm} just calls surface.msm with \code{type = "image"}.

  As these three functions are methods of the generic functions
  \code{contour}, \code{persp} and \code{image}, they
  can be invoked as \code{contour(x)},  \code{persp(x)} or
  \code{image(x)}, where \code{x} is a fitted \code{msm} object.  
}
\seealso{
  \code{\link{msm}}, \code{\link{contour}},
  \code{\link{filled.contour}}, \code{\link{persp}},
  \code{\link{image}}. 
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
