#' Generic function to find quantiles of a distribution
#' 
#' Generic function to find the quantiles of a distribution, given the
#' equivalent probability distribution function.
#' 
#' This function is intended to enable users to define \code{"q"} functions for
#' new distributions, in cases where the distribution function \code{pdist} is
#' available analytically, but the quantile function is not.
#' 
#' It works by finding the root of the equation \eqn{h(q) = pdist(q) - p = 0}.
#' Starting from the interval \eqn{(-1, 1)}, the interval width is expanded by
#' 50\% until \eqn{h()} is of opposite sign at either end.  The root is then
#' found using \code{\link{uniroot}}.
#' 
#' This assumes a suitably smooth, continuous distribution.
#' 
#' An identical function is provided in the \pkg{flexsurv} package.
#' 
#' @param pdist Probability distribution function, for example,
#' \code{\link{pnorm}} for the normal distribution, which must be defined in
#' the current workspace.  This should accept and return vectorised parameters
#' and values.  It should also return the correct values for the entire real
#' line, for example a positive distribution should have \code{pdist(x)==0} for
#' \eqn{x<0}.
#' @param p Vector of probabilities to find the quantiles for.
#' @param special Vector of character strings naming arguments of the
#' distribution function that should not be vectorised over. Used, for example,
#' for the \code{rate} and \code{t} arguments in \code{\link{qpexp}}.
#' @param ...  The remaining arguments define parameters of the distribution
#' \code{pdist}.  These MUST be named explicitly.
#' 
#' This may also contain the standard arguments \code{log.p} (logical; default
#' \code{FALSE}, if \code{TRUE}, probabilities p are given as log(p)), and
#' \code{lower.tail} (logical; if \code{TRUE} (default), probabilities are P[X
#' <= x] otherwise, P[X > x].).
#' 
#' If the distribution is bounded above or below, then this should contain
#' arguments \code{lbound} and \code{ubound} respectively, and these will be
#' returned if \code{p} is 0 or 1 respectively.  Defaults to \code{-Inf} and
#' \code{Inf} respectively.
#' @return Vector of quantiles of the distribution at \code{p}.
#' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
#' @keywords distribution
#' @examples
#' 
#' qnorm(c(0.025, 0.975), 0, 1)
#' qgeneric(pnorm, c(0.025, 0.975), mean=0, sd=1) # must name the arguments
#' 
#' @export qgeneric
qgeneric <- function(pdist, p, special=NULL, ...)
{
    args <- list(...)
    args.special <- if (!is.null(special)) args[special] else NULL
    args <- args[!names(args) %in% special]
    if (is.null(args$log.p)) args$log.p <- FALSE
    if (is.null(args$lower.tail)) args$lower.tail <- TRUE
    if (is.null(args$lbound)) args$lbound <- -Inf
    if (is.null(args$ubound)) args$ubound <- Inf
    if (args$log.p) p <- exp(p)
    if (!args$lower.tail) p <- 1 - p
    ret <- numeric(length(p))
    ret[p == 0] <- args$lbound
    ret[p == 1] <- args$ubound
    args[c("lower.tail","log.p","lbound","ubound")] <- NULL
    ## Other args assumed to contain params of the distribution.
    ## Replicate all to their maximum length, along with p 
    maxlen <- max(sapply(c(args, p=list(p)), length))
    for (i in seq_along(args))
        args[[i]] <- rep(args[[i]], length.out=maxlen)
    p <- rep(p, length.out=maxlen)

    ret[p < 0 | p > 1] <- NaN
    ind <- (p > 0 & p < 1)
    if (any(ind)) {
        hind <- seq_along(p)[ind]
        h <- function(y) {
            args <- lapply(args, function(x)x[hind[i]])
            p <- p[hind[i]]
            args$q <- y
            args <- c(args, args.special)
            (do.call(pdist, args) - p)
        }
        ptmp <- numeric(length(p[ind]))
        for (i in 1:length(p[ind])) {
            interval <- c(-1, 1)
            while (h(interval[1])*h(interval[2]) >= 0) {
              interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
            }
            ptmp[i] <- uniroot(h, interval, tol=.Machine$double.eps)$root
        }
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) warning("NaNs produced")
    ret
}
