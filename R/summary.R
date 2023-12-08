#' Summarise a fitted multi-state model
#' 
#' Summary method for fitted \code{\link{msm}} models. This is simply a wrapper
#' around \code{\link{prevalence.msm}} which produces a table of observed and
#' expected state prevalences for each time, and for models with covariates,
#' \code{\link{hazard.msm}} to print hazard ratios with 95\% confidence
#' intervals for covariate effects.
#' 
#' @name summary.msm
#' @aliases summary.msm print.summary.msm
#' @param object A fitted multi-state model object, as returned by
#' \code{\link{msm}}.
#' @param hazard.scale Vector with same elements as number of covariates on
#' transition rates. Corresponds to the increase in each covariate used to
#' calculate its hazard ratio. Defaults to all 1.
#' @param ... Further arguments passed to \code{\link{prevalence.msm}}.
#' @return A list of class \code{summary.msm}, with components:
#' 
#' \item{prevalences}{Output from \code{\link{prevalence.msm}}.}
#' 
#' \item{hazard}{Output from \code{\link{hazard.msm}}.}
#' 
#' \item{hazard.scale}{Value of the \code{hazard.scale} argument.}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}},\code{\link{prevalence.msm}},
#' \code{\link{hazard.msm}}
#' @keywords models
#' @export
summary.msm <- function(object, # fitted model
                        hazard.scale = 1,
                        ...
                        )
{
    if (!inherits(object, "msm")) stop("expected object to be a msm model")
    prevalences <- prevalence.msm(object, ...)
    if (object$qcmodel$ncovs > 0) {
        if (missing (hazard.scale))
            hazard.scale <- rep(1, object$qcmodel$ncovs)
        hazard <- hazard.msm(object)
    }
    else {hazard <- hazard.scale <- NULL}
    ret <- list(prevalences=prevalences,
                hazard=hazard,
                hazard.scale=hazard.scale)
    class(ret) <- "summary.msm"
    ret
}

print.summary.msm <- function(x,...)
{
    if (!is.null(x$prevalences)) {
        cat("\nObserved numbers of individuals occupying states at each time\n\n")
        print(x$prevalences$Observed)
        cat("\nExpected numbers of individuals occupying states at each time\n\n")
        print(x$prevalences$Expected)
        cat("\nObserved prevalences of states (percentages of population at risk)\n\n")
        print(x$prevalences$"Observed percentages")
        cat("\nExpected prevalences of states (percentages of population at risk)\n\n")
        print(x$prevalences$"Expected percentages")
    }
    i <- 1
    for (cov in names(x$hazard)) {
        cat ("\nTransition hazard ratios corresponding to covariate effects\n\n" )
        cat (cov, " ( unit of",x$hazard.scale[i],")\n")
        print(round(x$hazard[[cov]], 2))
        i <- i+1
    }
    invisible()
}
