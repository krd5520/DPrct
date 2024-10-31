
#' Approximate DP Estimate, Std.Error, z value, and Pr(>|z|) for a Model Coefficient
#' @param effect.vec vector of coefficients from the iterated fitted model
#' @inheritParams dp_pvalue
#' @inheritParams dp_confidence_interval
#' @param include.san.pvalue logical to indicate if z value and P(>|z|) are returned
#' @param include.ci logical to indicate if sanitized confidence interval is returned
#' @param san.sd (optional) a sanitized standard error for use in getting the
#'      sanitized range and then sanitized mean. If NULL, then \link{dp_estimate_sd}
#'      is used instead.
#'@returns if  \code{bound.mean!=NA} the vector of sanitized statistics, which include:
#'    the point estimate; estimated standard error; lower and upper confidence
#'    interval bounds and the confidence level (if \code{include.ci==TRUE});
#'    the z-value and p-value (if \code{include.san.pvalue==TRUE}); and the epsilon
#'    and delta privacy parameters. Otherwise, a vector of NA is returned.
#'
#' @seealso [dp_confidence_interval(), dp_estimate_sd(), dp_pvalue()]
#' @export
dp_coef_stats=function(effect.vec,
                       san.sd=NA,
                       bounds.sd=c(2^(-15),2^15),
                       bound.mean,
                       alphas=0.05,
                       epsilon,
                       delta=0,
                       pvalue.partitions=NA,
                       include.san.pvalue=TRUE,
                       include.ci=TRUE){
  if(is.na(bound.mean)==TRUE){ #no bound.mean then no privacy-preserving version
    if(include.san.pvalue==TRUE){
      cname=c("san.coef","san.sd","san.lower","san.upper","conf.level",
              "san.z.value","san.pvalue","epsilon","delta")
    }else if (include.ci==TRUE){
      cname=c("san.coef","san.sd","san.lower","san.upper","conf.level",
              "epsilon","delta")
    }else{
      cname=c("san.coef","san.sd","epsilon","delta")
    }
    san.stat.output=rep(NA,length(cname))
    names(san.stat.output)=cname
  }else{ # bound.mean value is supplied


    #2 privacy steps: san.range, san.point (san.point has no delta value)
    #if including pvalue, then there are +1 privacy steps
    #if no sanitized standard deviation included, then +1 privacy steps
    numprivsteps=2+as.numeric(include.san.pvalue==TRUE)+as.numeric(is.na(san.sd)==TRUE)

    #check that the epsilon and deltas are positive and non-negative
    stopifnot(sum(epsilon>0)==length(epsilon))
    stopifnot(sum(delta>=0)==length(delta))

    #if only 1 value given for epsilon (or delta), split it evenly across steps
    if(length(epsilon)==1){
      epsilon=rep(epsilon/numprivsteps,numprivsteps)
    }
    if(length(delta)==1){
      delta=rep(delta/(numprivsteps-1),numprivsteps-1)
    }

    #if only 1 alpha given, split it evenly across the confidence interval steps
    #if including p-value then let it have the same significance level as the CI
    if(include.san.pvalue==TRUE){
    if(length(alphas)==1){
        alphas=c(rep(alphas/3,3),alphas)
    }else if(length(alphas)==2){
      alphas=c(rep(alphas[1]/3,3),alphas[2])
    }else if(length(alphas)==3){
      alphas=c(alphas,sum(alphas))
    }
    }else{
      if(length(alphas)<3){
        alphas=rep(alphas[1]/3,3)
      }
    }


    #check that the CI alphas are between 0 and 1
    stopifnot(sum(alphas[1:3])>0&sum(alphas[1:3])<1)

    if(is.na(san.sd)==TRUE){ #getting sanitized standard deviation
      san.stats=dp_confidence_interval(effect.vec,
                                       epsilon[1:3],
                                       delta[1:2],
                                       alphas=alphas[1:3],
                                       bounds.sd=bounds.sd,
                                       bound.mean=bound.mean,
                                       return.point.sd=TRUE)
    }else{ #using supplied sanitized standard deviation
      san.stats=dp_confidence_interval(effect.vec,
                                       epsilon.vec=epsilon[1:2],
                                       delta.vec=delta[1],
                                       alphas=alphas[1:3],
                                       x.sd=san.sd,
                                       bound.mean=bound.mean,
                                       return.point.sd=TRUE)

    }

    san.stat.output=unlist(san.stats)
    if(include.ci==TRUE){
      #reorder so it is c(point, sd, lower, upper)
      san.stat.output=c(san.stat.output[3:4],san.stat.output[1:2],
                        "conf.level"=1-sum(alphas[1:3]))
    }else{
      san.stat.output=c(san.stat.output[3:4])
    }




    if(include.san.pvalue==TRUE){
      pvalue.test=function(x,est.sd){ #one-sided test
        #H0: |coef|=0 vs H1: |coef|>0
        stats::pnorm(base::abs(mean(x)),
                     sd=est.sd/base::sqrt(length(x)),
                     lower.tail=F)
      }

      san.pvalue=dp_pvalue(effect.vec,
                           epsilon=epsilon[numprivsteps],
                           delta=delta[numprivsteps-1],
                           test=pvalue.test,
                           pvalue.partitions=pvalue.partitions,
                           alpha=alphas[4]/2) #alpha/2 for two-sided test
      san.pvalue=unlist(san.pvalue)
      #reorder san.pvalue output to be test statistic then pvalue
      san.stat.output=c(san.stat.output,san.pvalue[2],san.pvalue[1])
    }
    san.stat.output=c(san.stat.output,
                      c("epsilon"=sum(epsilon),
                        "delta"=sum(delta)))
  } #end else bound.mean is supplied

  return(san.stat.output)
}
