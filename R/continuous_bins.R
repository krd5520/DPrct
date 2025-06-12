#' Internal-use function to discretize continuous columns in a certain number of bins
#'
#' @param cont.col a vector of continuous values to be split into bins
#' @param num.bins a integer value for the number of bins to split cont.col values into. If supplied.
#' @param bin.breaks a vector of numeric values to define the boundaries of the bins
#' @param bin.lab a vector of names for the bins. If \code{bin.lab=NULL}, then midpoints are used.
#' @param cont.limits a vector of the minimum theoretical value and max theoretical value (or truncation values)
#' @return a vector of categorical variables taking the value of the midpoint of the interval
#'
#' @examples
#' continuous_bins(sample(1:100,300,replacement=T),num.bins=3)
#'
#' @noRd
#' @keywords internal

continuous_bins<-function(cont.col,num.bins=NA,bin.breaks=NULL,bin.lab=NULL,cont.limits=NULL){
  #check that number of bins or vector of breaks is provided
  if(base::is.na(num.bins)&base::is.null(bin.breaks)){
    stop("num.bins or bin.breaks must be inputted.")
  }
  cont.col=unlist(cont.col)

  na.idx=is.na(cont.col)
  count.col=as.numeric(cont.col)
  if(base::is.null(bin.breaks)==FALSE){ #if num.bins not provided, use bin.breaks
    if(base::is.null(bin.lab)==TRUE){ #but no bin.labels supplied. Then midpoints are labels
      bin.lab=base::as.character(0.5*(bin.breaks[base::seq(2,base::length(bin.breaks))]-
                                        bin.breaks[base::seq(1,base::length(bin.breaks)-1)]))
    }

    #discretize the continuous data
    out.col<-cut(cont.col[!na.idx],breaks=bin.breaks,label=base::as.character(bin.lab),include.lowest=T)
  }else{ #if num.bins is given
    if(is.list(cont.limits)==TRUE){
      cont.limits=unlist(cont.limits)
    }
    #stopifnot(length(cont.limits)==2)
    out.limits=sum(cont.col[!na.idx]<cont.limits[1])+sum(cont.col[!na.idx]>cont.limits[2])
    if(out.limits>0){
      warning(paste("Limits on continuous variable (",paste(cont.limits,collapse = ", "),") do not contain all observed values.",out.limits,"observations have been clipped to fit limits."))
      cont.col[((!na.idx)&(cont.col[!na.idx]<cont.limits[1]))]=cont.limits[1]
      cont.col[((!na.idx)&(cont.col[!na.idx]>cont.limits[2]))]=cont.limits[2]
    }
    bin.breaks=seq(cont.limits[1],cont.limits[2],length.out=num.bins+1)
    if(base::is.null(bin.lab)==T){
      bin.lab=as.character(bin.breaks[-1]-((cont.limits[2]-cont.limits[1])/(2*num.bins)))
    }
    cat.col=base::cut(cont.col[!na.idx],breaks=bin.breaks,label=bin.lab,include.lowest=T) #cut into intervals
    #if(base::is.null(bin.lab)==TRUE){ #if no label provided.
    #  intervals=base::levels(cat.col) #get levels for intervals
    #  #get the boundary points for the intervals
    #  boundaries= base::as.numeric(
    #    unique(base::gsub("\\(|\\]|\\[|\\)","",
    #               base::unlist(base::strsplit(intervals,",")))))
    #  boundaries[1]=base::min(cont.col) #change first lower bound to be minimum cont.col value
    #  bin.labs<-0.5*(boundaries[base::seq(2,num.bins+1)]+boundaries[base::seq(1,num.bins)])
    #}
    out.col=rep(NA,length(cont.col))
    out.col[!na.idx]=as.character(cat.col)
    out.col=factor(out.col,levels=as.character(bin.lab))
  }

  return(out.col)
}
