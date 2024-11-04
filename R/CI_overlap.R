#' Confidence Measure Overlap
#'
#' @param lw1 is the lower bound of the first confidence interval (vector or single value)
#' @param lw2 is the lower bound of the comparison confidence interval (vector or single value)
#' @param up1 is the upper bound of the first confidence interval (vector or single value).
#'    \code{up1} must have the same length as \code{lw1}.
#' @param up2 is the upper bound of the comparison confidence interval (vector or single value).
#'    \code{up2} must have the same length as \code{lw2}.
#' @returns the measure of confidence interval overlap.
#' @details
#' The length of \code{lw1} or \code{lw2} must be a multiple of the other.
#'
#'
CI_overlap=function(lw1,up1,lw2,up2)
{
  stopifnot((length(lw1)==length(up1))&(length(lw2)==length(up2)))
  lw.mx=base::pmax(lw1,lw2)
  up.mn=base::pmin(up1,up2)
  dif.upmn.lwmx=up.mn-lw.mx
  ci.len1=up1-lw1
  ci.len2=up2-lw2
  ci.overlap=dif.upmn.lwmx*0.5*((1/ci.len1)+(1/ci.len2))
  if((length(ci.overlap)==1)&&(dif.upmn.lwmx<0)){
    ci.overlap=0
  }else{
    ci.overlap[dif.upmn.lwmx<0]=0
  }
    return(ci.overlap)
}
