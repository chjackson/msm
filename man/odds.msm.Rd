\name{odds.msm}
\title{Calculate tables of odds ratios for covariates on
  misclassification probabilities}
\alias{odds.msm}
\description{
  Odds ratios are computed by exponentiating the
  estimated covariate effects on the logit-misclassification
  probabilities.
}
\usage{
odds.msm(x, odds.scale = 1, cl = 0.95)
}
\arguments{
  \item{x}{Output from \code{\link{msm}} representing a fitted
    multi-state model.}

  \item{odds.scale}{Vector with same elements as number of covariates
    on misclassification probabilities. Corresponds to the increase in
    each covariate used to calculate its odds ratio. Defaults to all 1.}

  \item{cl}{Width of the symmetric confidence interval to present. 
    Defaults to 0.95.}

}
\value{

  A list of tables containing odds ratio estimates, one table for each covariate.
  Each table has three columns, containing the odds ratio, and an
  approximate upper 95\% and lower 95\% confidence limit respectively
  (assuming normality on the log scale), for each misclassification
  probability.

}
\seealso{
  \code{\link{msm}}, \code{\link{hazard.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
