#' Convert a hmodel object to HMM constructor function calls
#' 
#' @param hmodel A list of class \code{hmodel}, as returned in the \code{hmodel} component of the
#' fitted model object from \code{\link{msm}}.
#'
#' @param hmmdist \code{TRUE} or \code{FALSE} (see "Value" section).
#'
#' @returns 
#' 
#' If \code{hmmdist=TRUE}, returns a list of objects of class \code{hmmdist}. 
#'  These are the kind of objects
#' returned by HMM constructor functions such
#' as \code{\link{hmmNorm}}, \code{\link{hmmPois}} etc.  Therefore the list can be
#' passed as the \code{hmodel} argument to \code{\link{msm}}.
#' 
#' If \code{hmmdist=FALSE}, returns a list comprised of the corresponding input
#' arguments for the constructor functions, i.e. parameter values of HMM emission
#' distributions.   The list has one element per state.  Each of these elements 
#' has one element per parameter (for univariate HMMs), or one element per outcome
#' distribution, which in turn has one element per parameter (for multivariate HMMs).
#'
#' @author Will Hulme \code{https://github.com/wjchulme} and Chris Jackson. 
#' 
#' @export
hmodel2list <- function(hmodel, hmmdist = TRUE){

  if(!hmodel$hidden) stop("`hmodel` is not an msm hidden Markov model object")

  ## makes a state-specific vector of parameters extracted from hmodel
  ## into a list of parameters (treating hmmCat as a special case)

  makeargslist <- function(params, label){
    ## params = named vector of parameters for the distribution function
    ## label = label (character) of the distribution function
    if(!(label %in% .msm.HMODELFNS$label))
      stop("Distribution ", label, " not currently supported for hmodel2list")
    if(label=="categorical")
      list(prob = params[names(params) %in% c("p", "p0", "pbase")],
           basecat = params[names(params)=="basecat"])
    else if(label=="identity")
      list(x = params[names(params) == "which"])
    else
      as.list(params)
  }

  labellist <- as.list(na.omit(as.vector(hmodel$labels)))
  
  paramlist <- split(hmodel$pars, list(hmodel$parout, hmodel$parstate))
  paramlist <- paramlist[sapply(paramlist, length)>0] # for e.g. hmmMV mixed with hmmIdent
  
  paramnestedlist <- mapply(makeargslist, paramlist, labellist,
                            SIMPLIFY=FALSE, USE.NAMES=FALSE)
  distlist <- lapply(labellist, function(label){
    match.fun(.msm.HMODELFNS$hmmname[.msm.HMODELFNS$label==label])
  })

  if(hmodel$mv){

    hmmdistlist <- invoke_map_base(distlist, paramnestedlist)
    
    hmmdistnestedlist <- split(hmmdistlist, 
                               rep(seq_len(hmodel$nstates), times=hmodel$nout))

    msmlist <- lapply(hmmdistnestedlist, 
                      function(hmmdist){lift_dl_base(msm::hmmMV)(hmmdist)})

    if(hmmdist)
      msmlist
    else
      split(paramnestedlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
  } else {
    if(hmmdist)
      invoke_map_base(distlist, paramnestedlist)
    else
      paramnestedlist
  }

}

## base R rewrites of functions originally in purrr

invoke_map_base <- function(fnlist, paramlist){
  reslist <- vector(length(fnlist), mode="list")
  for (i in seq_along(reslist))
    reslist[[i]] <- do.call(fnlist[[i]], paramlist[[i]])
  reslist
}

lift_dl_base <- function (..f, ..., .unnamed = FALSE) {
    force(..f)
    defaults <- list(...)
    function(.x = list(), ...) {
      if (.unnamed) {
        .x <- unname(.x)
      }
      do.call("..f", c(.x, defaults, list(...)))
    }
}
