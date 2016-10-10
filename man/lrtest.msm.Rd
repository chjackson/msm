\name{lrtest.msm}
\alias{lrtest.msm}
\title{Likelihood ratio test}
\description{Likelihood ratio test between two or more fitted multi-state models}
\usage{
lrtest.msm(...)
}
\arguments{

  \item{...}{Two or more fitted multi-state models, as returned by
    \code{\link{msm}}, ordered by increasing numbers of parameters.}

}
\value{
  A matrix with three columns, giving the likelihood ratio statistic,
  difference in degrees of freedom and the chi-squared p-value for a
  comparison of the first model supplied with each subsequent model.
}
\section{Warning}{
  The comparison between models will only be valid if
  they are fitted to the same dataset. This may be a problem if
  there are missing values and R's default of 'na.action = na.omit'
  is used.

  The likelihood ratio statistic only has the indicated chi-squared
  distribution if the models are nested. An alternative for comparing
  non-nested models is Akaike's information criterion.  This can be 
  computed for one or more fitted \code{msm} models \code{x,y,...} using
  \code{\link{AIC}(x,y,...)}.
}
\seealso{\code{\link{logLik.msm}},\code{\link{msm}}}
\keyword{models}
