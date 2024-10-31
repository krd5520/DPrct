#' Internal-use function to parse the returntype variable for the three algorithm code functions.
#' @inheritParams hybrid_synth
#' @param synth.priv.cost is a named vector with the epsilon and delta totals
#'    for the generation of the synthetic data
#' @param treat.stat.epsilon a vector or single value for the epsilon used for the
#'    mechanisms to return the privacy preserving statistics about the treatment effect.
#' @param treat.stat.delta a vector or single value for the delta used for the
#'    mechanisms to return the privacy preserving statistics about the treatment effect.
#' @returns list of objects to be returned from \code{hybrid_dp}, \code{crude_fully_dp}, or \code{advanced_fully_dp}
#' @keywords internal
#' @export
output_returntypes<-function(returntypes,
                             synth.data=NULL,
                             san.treatment.effect=NA,
                             san.model=NULL,
                             treatment.var=NULL,
                             confidential.model=NULL,
                             start.time=NULL,
                             time.values=NULL,
                             synth.priv.cost=NULL,
                             ...){
  dup.returntypes=duplicated(returntypes)
  if(sum(dup.returntypes)>0){
    returntypes=returntypes[!dup.returntypes]
    warning(paste0("Duplicate return types listed. Returning: ",
                   paste(returntypes,collapse=", ")))
  }
  output<-list()
  eps.account=0
  del.account=0
  for(rtype in returntypes){
    if(rtype=="synth.data"){ #if returning synthetic data, check it is supplied
      output=append(output,list("synth.data"=synth.data))
      if(is.null(attr(synth.data,"priv.cost"))==FALSE){
        eps.account=eps.account+attr(synth.data,"priv.cost")[1]
        del.account=eps.account+attr(synth.data,"priv.cost")[2]
      }else{
        names(synth.priv.cost)=NULL
        eps.account=synth.priv.cost[1]+eps.account
        del.account=del.account+synth.priv.cost[2]
      }
    }


    if(rtype=="san.model"){
      output=append(output,list("san.model"=san.model))
    }
    if(rtype=="confidential.model"){
      output=append(output,list("confidential.model"=confidential.model))
    }

  }
  if(length(returntypes)==1){ #if one element or less, force to not be a list
    output=unlist(output)
  }
  if("comp.time"%in%returntypes){
    stop.time=proc.time()
    output=append(output,list("comp.time"=c("total_time"=(stop.time-start.time)[[3]],time.values)))

  }
  attr(output,"priv.cost")=c("epsilon"=eps.account,"delta"=del.account)
  return(output)
}
