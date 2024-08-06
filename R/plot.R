### Estimated survival probability from each state



#' Plots of multi-state models
#' 
#' This produces a plot of the expected probability of survival against time,
#' from each transient state. Survival is defined as not entering an absorbing
#' state.
#' 
#' Note that while this function is only relevant to models with absorbing
#' states, models in \pkg{msm} can have any transition structure and do not
#' necessarily have to have an absorbing state.
#' 
#' 
#' @param x Output from \code{\link{msm}}, representing a fitted multi-state
#' model object.
#' @param from States from which to consider survival. Defaults to the complete
#' set of transient states.
#' @param to Absorbing state to consider. Defaults to the highest-labelled
#' absorbing state.
#' @param range Vector of two elements, giving the range of times to plot for.
#' @param covariates Covariate values for which to evaluate the expected
#' probabilities.  This can either be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' @param legend.pos Vector of the \eqn{x} and \eqn{y} position, respectively,
#' of the legend.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param lwd Line width. See \code{\link{par}}.
#' @param ... Other arguments to be passed to the generic \code{\link{plot}}
#' and \code{\link{lines}} functions.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}
#' @keywords models
#' @export
plot.msm <- function(x, from=NULL, to=NULL, range=NULL, covariates="mean", legend.pos=NULL, xlab="Time", ylab="Fitted survival probability", lwd=1,...)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (is.null(from))
        from <- transient.msm(x)
    else {
        if (!is.numeric(from)) stop("from must be numeric")
        if (any (! (from %in% 1:x$qmodel$nstates ) ) )
            stop("from must be a vector of states in 1, ..., ", x$qmodel$nstates)
    }
    if (is.null(to)){
        if (length(absorbing.msm(x))==0)
            stop("\"to\" not specified, and no absorbing state. See help(plot.msm)")
        to <- max(absorbing.msm(x))
    }
    else {
        if (!is.numeric(to)) stop("to must be numeric")
        if (! (to %in% absorbing.msm(x) ) ) stop("to must be an absorbing state")
    }
    if (is.null(range))
        rg <- range(model.extract(x$data$mf, "time"))
    else {
        if (!is.numeric(range) || length(range)!= 2) stop("range must be a numeric vector of two elements")
        rg <- range
    }
    timediff <- (rg[2] - rg[1]) / 50
    times <- seq(rg[1], rg[2], timediff)
    pr <- numeric()
    cols <- rainbow(length(from))
    for (t in times)
        pr <- c(pr, pmatrix.msm(x, t, times[1], covariates)[from[1], to])
    plot(times, 1 - pr, type="l", xlab=xlab, ylab=ylab, lwd=lwd,
         ylim=c(0,1), lty = 1, col=cols[1],...)
    lt <- 2
    for (st in from[-1]){
        pr <- numeric()
        for (t in times)
            pr <- c(pr, pmatrix.msm(x, t, times[1], covariates)[st, to])
        lines(times, 1 - pr, type="l", lty = lt, lwd=lwd, col=cols[lt],...)
        lt <- lt+1
    }
    if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(max(times) - 15*timediff, 1)
    legend(legend.pos[1], legend.pos[2], legend=paste("From state",from), lty = seq(lt-1), col=cols, lwd=lwd)
    invisible()
}

### Plot KM estimate of time to first occurrence of each state



#' Kaplan Meier estimates of incidence
#' 
#' Compute and plot Kaplan-Meier estimates of the probability that each
#' successive state has not occurred yet.
#' 
#' If the data represent observations of the process at arbitrary times, then
#' the first occurrence of the state in the data will usually be greater than
#' the actual first transition time to that state.  Therefore the probabilities
#' plotted by \code{\link{plotprog.msm}} will be overestimates.
#' 
#' @param formula A formula giving the vectors containing the observed states
#' and the corresponding observation times. For example,
#' 
#' \code{state ~ time}
#' 
#' Observed states should be in the set \code{1, \dots{}, n}, where \code{n} is
#' the number of states.
#' @param subject Vector of subject identification numbers for the data
#' specified by \code{formula}. If missing, then all observations are assumed
#' to be on the same subject. These must be sorted so that all observations on
#' the same subject are adjacent.
#' @param data An optional data frame in which the variables represented by
#' \code{state}, \code{time} and \code{subject} can be found.
#' @param legend.pos Vector of the \eqn{x} and \eqn{y} position, respectively,
#' of the legend.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param lwd Line width. See \code{\link{par}}.
#' @param xlim x axis limits, e.g. c(0,10) for an axis ranging from 0 to 10.
#' Default is the range of observation times.
#' @param mark.time Mark the empirical survival curve at each censoring point,
#' see \code{\link[survival]{lines.survfit}}.
#' @param ... Other arguments to be passed to the \code{\link{plot}} and
#' \code{\link[survival]{lines.survfit}} functions.
#' @seealso \code{\link[survival]{survfit}},
#' \code{\link[survival]{plot.survfit}}
#' @keywords models
#' @export
plotprog.msm <- function(formula, subject, data, legend.pos=NULL, xlab="Time", ylab="1 - incidence probability", lwd=1, xlim=NULL,
                         mark.time=TRUE, ...) {
    data <- na.omit(data)
    mf <- model.frame(formula, data=data)
    state <- mf[,1]
    time <- mf[,2]
    if (!is.null(data))
        subject <- eval(substitute(subject), as.list(data), parent.frame())
    subject <- match(subject, unique(subject))
    rg <- range(time)
    if (is.null(xlim)) xlim=rg
    plot(0, xlim=xlim, ylim=c(0,1), type="n", xlab=xlab, ylab=ylab, ...)
    states <- sort(unique(state))[-1]
    cols <- rainbow(length(states))
    for (i in states) {
        dat <- cbind(subject, time, state)
        st <- as.data.frame(
                            do.call("rbind", by(dat, subject, function(x)
                                            {
                                                c(anystate = if(any(x[,"state"]>=i)) 1 else 0,
                                                  mintime = if(any(x[,"state"]>=i)) min(x[x[,"state"] >= i, "time"]) else max(x[,"time"]))
                                            }
                                                ))) # slow
        lines(survfit(Surv(st$mintime,st$anystate) ~ 1),
              col=cols[i-1], lty=i-1, lwd=lwd, mark.time=mark.time, ...)
    }
    timediff <- (rg[2] - rg[1]) / 50
    if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(max(time) - 25*timediff, 1)
    legend(legend.pos[1], legend.pos[2], lty=states-1, lwd=lwd, col=cols,
           legend=paste("To state", states, c(rep("or greater", length(states)-1), "")))
    invisible()
}


### Likelihood surface plots



#' Explore the likelihood surface
#' 
#' Plot the log-likelihood surface with respect to two parameters.
#' 
#' Draws a contour or perspective plot.  Useful for diagnosing irregularities
#' in the likelihood surface.  If you want to use these plots before running
#' the maximum likelihood estimation, then just run \code{msm} with all
#' estimates fixed at their initial values.
#' 
#' \code{contour.msm} just calls surface.msm with \code{type = "contour"}.
#' 
#' \code{persp.msm} just calls surface.msm with \code{type = "persp"}.
#' 
#' \code{image.msm} just calls surface.msm with \code{type = "image"}.
#' 
#' As these three functions are methods of the generic functions
#' \code{contour}, \code{persp} and \code{image}, they can be invoked as
#' \code{contour(x)}, \code{persp(x)} or \code{image(x)}, where \code{x} is a
#' fitted \code{msm} object.
#' 
#' @aliases surface.msm persp.msm contour.msm image.msm
#' @param x Output from \code{\link{msm}}, representing a fitted msm model.
#' @param params Integer vector with two elements, giving the indices of the
#' parameters to vary. All other parameters will be fixed. Defaults to
#' \code{c(1,2)}, representing the first two log transition intensities. See
#' the \code{fixedpars} argument to \code{msm} for a definition of these
#' indices.
#' @param np Number of grid points to use in each direction, by default 10.  An
#' \code{np x np} grid will be used to evaluate the likelihood surface. If 100
#' likelihood function evaluations is slow, then reduce this.
#' @param type Character string specifying the type of plot to produce.
#' \tabular{ll}{ \code{"contour"} \tab Contour plot, using the R function
#' \code{\link{contour}}. \cr \code{"filled.contour"} \tab Solid-color contour
#' plot, using the R function \code{\link{filled.contour}}. \cr \code{"persp"}
#' \tab Perspective plot, using the R function \code{\link{persp}}. \cr
#' \code{"image"} \tab Grid color plot, using the R function
#' \code{\link{image}}. \cr }
#' @param point Vector of length \code{n}, where \code{n} is the number of
#' parameters in the model, including the parameters that will be varied here.
#' This specifies the point at which to fix the likelihood.  By default, this
#' is the maximum likelihood estimates stored in the fitted model \code{x},
#' \code{x$estimates}.
#' @param xrange Range to plot for the first varied parameter.  Defaults to
#' plus and minus two standard errors, obtained from the Hessian at the maximum
#' likelihood estimate.
#' @param yrange Range to plot for the second varied parameter.  Defaults to
#' plus and minus two standard errors, obtained from the Hessian at the maximum
#' likelihood estimate.
#' @param ... Further arguments to be passed to the plotting function.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{msm}}, \code{\link{contour}},
#' \code{\link{filled.contour}}, \code{\link{persp}}, \code{\link{image}}.
#' @keywords models
#' @export
surface.msm <- function(x, params=c(1,2), np=10, type=c("contour","filled.contour","persp","image"),
                        point=NULL, xrange=NULL, yrange=NULL,...)
{
    type <- match.arg(type)
    if (is.null(point))
        point <- x$paramdata$opt$par
    se <- sqrt(diag(x$covmat[x$paramdata$optpars,x$paramdata$optpars]))
    i1 <- params[1]; i2 <- params[2]
    if (is.null(xrange)) {
        pmin <- point[i1] - 2*se[i1]
        pmax <- point[i1] + 2*se[i1]
        p1 <- seq(pmin, pmax, length.out=np)
    }
    else p1 <- seq(xrange[1], xrange[2], length.out=np)
    if (is.null(yrange)){
        pmin <- point[i2] - 2*se[i2]
        pmax <- point[i2] + 2*se[i2]
        p2 <- seq(pmin, pmax, length.out=np)
    }
    else p2 <- seq(yrange[1], yrange[2], length.out=np)

    z <- matrix(nrow=np, ncol=np)
    for (i in 1:np) {
        for (j in 1:np) {
            point[i1] <- p1[i]; point[i2] <- p2[j]
            z[i,j] <- -0.5*Ccall.msm(point, "lik", expand.data(x), x$qmodel, x$qcmodel, x$cmodel, x$hmodel, x$paramdata)
        }
    }

    switch(type,
           contour = contour(p1, p2, z, ...),
           filled.contour = filled.contour(p1, p2, z, ...),
           image = image(p1, p2, z, ...),
           persp = persp(p1, p2, z, zlab="Log-likelihood",...)
           )
    invisible()
}

#' @rdname surface.msm
#' @export 
contour.msm <- function(x, ...)
{
    surface.msm(x, type="contour",...)
}

#' @rdname surface.msm
#' @export 
persp.msm <- function(x, ...)
{
    surface.msm(x, type="persp",...)
}

#' @rdname surface.msm
#' @export 
image.msm <- function(x, ...)
{
    surface.msm(x, type="image",...)
}
