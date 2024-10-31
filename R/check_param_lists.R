
#' Internal Function: Check lengths of Lists for Necessary Parameters
#' @param eps.list list of epsilons
#' @param del.list list of deltas
#' @param alph.list list of alphas
#' @param bd.sd.list list of bounds.sd tuples
#' @param bd.mean.list list of bound.mean values
#' @param ppart.list list of pvalue.partitions
#' @param num.coefs number of coefficients to be made privacy-preserving
#' @returns list of the lists, ensuring that each list has
#'
#' @keywords internal
#' @noRd
check_lists_coeftable=function(eps.list,
                               del.list,
                               alph.list,
                               bd.sd.list,
                               bd.mean.list,
                               ppart.list,
                               num.coefs){

  list.names=c("epsilon",
               "delta",
               "alpha",
               "bounds.sd",
               "bound.mean",
               "pvalue.partition")
  all.lists=list(eps.list,
                 del.list,
                 alph.list,
                 bd.sd.list,
                 bd.mean.list,
                 ppart.list)
  out=lapply(seq(1,6),function(idx)check_length_list(all.lists[[idx]],list.names[[idx]],num.coefs))

}

#' Internal Check List of Parameters
#' @param ls is list
#' @param ls.name is name of list
#' @param m is expected length of list
#' @returns list of correct length
#' @keywords internal
#' @noRd
check_length_list=function(ls,ls.name,m){
  warn=FALSE
  if(is.list(ls)==FALSE){
    ls=list(ls)
  }
  if(length(ls)<m){
    if(length(ls)>1){
      warn=TRUE
    }
    if((ls.name=="epsilon"|ls.name=="delta")&length(ls)==1){
      ls=rep(list(unlist(ls)[1]/m),m)
    }else{
      ls=rep(ls,m)
    }
  }
  if(warn==TRUE){
    message(paste0("List of ",ls.name,
                   " values is less than the number of model coefficients to be addressed. ",
                   "The following list is used instead: ",
                   paste0(ls[1:m],collapse=", ")))
  }
  return(ls[1:m])
}
