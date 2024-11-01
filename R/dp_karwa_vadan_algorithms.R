#' Sanitized max and minimum satisfying Approximate-DP with known variance.
#' @description Algorithm 1 by Karwa and Vadhan(2017)
#'
#' @param x vector of continuous numeric confidential values.
#' @param sd a public standard deviation of the population \code{x} is sampled from.
#'      A sanitized estimate of it may be used (see details.)
#' @param epsilon a positive privacy parameter
#' @param delta a privacy parameter. Default is 0
#' @param bound.mean finite values such that the mean of the population \code{x}
#'     is sampled from is bounded by (-\code{bound.mean}, \code{bound.mean})
#' @param range.prob numeric value so the probability that
#'     \eqn{P(\bigcup_{i=1}^n X_{min}\leq X_i\leq X_{max})\geq 1-}$\code{range.prob}
#'     where $X_i$ are i.i.d. normal random variable with mean and variance equal
#'     to that of the population \code{x} is sampled from and $n$ is the length of \code{x}.
#' @returns a vector with sanitized minimum and maximum of \code{x}.
#'
#' @details
#' If \code{sd} is unknown or not public, the \code{dp_estimate_sd()} function
#' can be used to get a sanitized estimate of the standard deviation to use in this function.
#'
#' @examples
#' set.seed(1)
#' x=rnorm(500,70,5)
#' print(c(min(x),max(x)))
#' dp_range(x=x,sd=5,epsilon=1,delta=0.001,bound.mean=1e6,range.prob=0.025)

#'
#'
#' @references{
#' Karwa, V. and Vadhan, S. (2017). Finite sample differentially private confidence intervals.
#'
#' }
#'
#'
#' @family ApproximateDPmeasures
#' @importFrom VGAM rlaplace
#'
#'@export

dp_range<-function(x,sd,epsilon,delta=0,bound.mean,range.prob){
  n=length(x)
  stopifnot(is.numeric(range.prob)&&range.prob>0&range.prob<=1)

  bound=ceiling(abs(bound.mean)/sd) #r in paper
  iter.bound.breaks=seq(-bound,bound) #j=-r,-r+1,...,r
  #breaks for 2r+1 bins in [-R-sd/2,R+sd/2] where R is the bound.mean
  breaks.vec=c(-Inf,sd*(iter.bound.breaks-0.5),Inf)
  bin.labs=as.character(c(iter.bound.breaks,
                          iter.bound.breaks[length(iter.bound.breaks)]+1))

  #break the data into bins and get frequencies
  hist.df=data.frame(table(
    continuous_bins(x,bin.breaks=breaks.vec,bin.lab=bin.labs)))
  hist.df=dp_perturbed_hist(hist.df=hist.df,epsilon=epsilon,delta=delta)
  #l-hat in paper. This is the bin that has the highest sanitized proportion
  biggest.san.bin=as.numeric(as.character(hist.df[which.max(hist.df$san.prop),1]))

  #get sanitized min and max
  half.range=4*sd*base::sqrt(base::log(n/range.prob))
  san.range.center=sd*biggest.san.bin
  san.range=c("san.min"=san.range.center-half.range,
                "san.max"=san.range.center+half.range)
  attr(san.range,"priv.cost")=c("epsilon"=epsilon,"delta"=delta)
  return(san.range)
}

#' Sanitized estimate of standard deviation satisfying Approximate-DP
#' @description Algorithm 2 by Karwa and Vadhan (2017)
#'
#' @param x vector of continuous numeric confidential values.
#' @param bounds.sd a vector of two elements: the minimum and the maximum the standard deviation of the population \code{x} is sampled from
#' @param epsilon a positive privacy parameter
#' @param delta a nonnegative privacy parameter. Defaul is \code{delta=0}
#' @returns a sanitized estimated variance of \code{x}
#'
#' @examples
#' set.seed(1)
#' x=rnorm(500,70,5)
#' dp_estimate_sd(x=x,epsilon=1,delta=0.001,bounds.sd=c(2^(-15),2^15))
#'
#'
#' @family ApproximateDPmeasures
#' @references{
#' Karwa, V. and Vadhan, S. (2017). Finite sample differentially private confidence intervals.
#'
#' }
#'
#'
#' @importFrom VGAM rlaplace
#'
#' @export
#'

dp_estimate_sd<-function(x,epsilon,delta=0,bounds.sd){
  n=length(x)

  #j=jmin,....jmax
  iter.bound.breaks=seq(floor(base::log2(bounds.sd[1]))-2,
                        ceiling(base::log2(bounds.sd[2])+1))
  #breaks for 2r+1 bins in [-R-sd/2,R+sd/2] where R is the bound.mean
  breaks.vec=c(-Inf,2^iter.bound.breaks,Inf)
  length(breaks.vec)

  half.n=n%/%2 #half rounded down
  x=x[sample(n,n)] #randomize the order of the x
  x.diffs=x[seq(1,half.n)]-x[seq(half.n+1,2*half.n)]

  bin.lab=c(iter.bound.breaks,
            iter.bound.breaks[length(iter.bound.breaks)]+1)

  #break the data into bins and get frequencies
  hist.df=data.frame(table(
    continuous_bins(abs(x.diffs),bin.breaks=breaks.vec,bin.lab=bin.lab)))
  hist.df=dp_perturbed_hist(hist.df=hist.df,epsilon=epsilon,delta=delta)

  #l-hat in paper. This is the bin that has the highest sanitized probability
  biggest.san.bin=as.numeric(as.character(hist.df[which.max(hist.df$san.prop),1]))
  san.sd=2^(biggest.san.bin+2)
  attr(san.sd,"priv.cost")=c("epsilon"=epsilon,"delta"=delta)
  return(san.sd)

}

#' Confidence Intervals satisfying Approximate-DP
#' @description Algorithm 4 by Karwa and Vadhan, (2017)
#' @inheritParams dp_range
#' @inheritParams dp_estimate_sd
#' @param epsilon.vec a positive value or a vector of them for the privacy parameters. (See details).
#' @param delta.vec a nonnegative value or a vector of them for the privacy parameters. Default is 0. (See details).
#' @param alphas a single value between (0,1) to be divided evenly by into 3 or
#'      three values whose sum is between (0,1) such that the confidence interval is
#'     \code{sum(alphas)}-level confidence interval. The first value is used in
#'     the \eqn{1-\alpha/2} quantile of a standard normal. The second is used in
#'     \eqn{b_1 log(1/\alpha_1)} where \eqn{b_1} is the scale of the Laplace noise
#'     for the point estimate. The third is the \code{range.prob} for the
#'     \code{dp_range()} function.
#' @param return.point.sd logical indicator for whether the privacy-preserving point
#'      and standard deviation estimate of treatment effect should be returned.
#' @param x.sd (optional) the publicly known (or privacy-preserving estimate of)
#'      standard deviation of the population \code{x} is sampled from. If it is NA,
#'      then the \code{dp_estimate_sd} function will be used to get a DP estimate of it.
#'      In this case, \code{bounds.sd} must be a tuple of the bounds on the
#'      standard deviation. (See details for specifications of \code{epsilon.vec}
#'      and \code{delta.vec} when \code{x.sd=NA})
#' @param san.point (optional) the sanitized point estimate of the coefficient.
#'      If \code{san.point=NA} (default), then \code{\link{dp_range}} will be used
#'      with the \code{\link{laplace_mechanism}} to get a sanitized point estimate.
#' @returns a vector of the lower and upper bounds of the
#'      \code{sum(alphas)}-level confidence interval that satisfies
#'      (\code{sum(epsilon.vec)}, \code{sum(delta.vec)})-DP, or a list of the
#'      privacy-preserving confidence interval with the privacy-preserving point
#'      and standard deviation estimates for the treatment.
#'
#' @details
#' For the \code{epsilon.vec} and \code{delta.vec}, If the \code{x.sd=NA} is
#'     provided, \code{epsilon.vec} should be 3 values (or 1 value to be split evenly
#'     across 3 privacy steps) and the \code{delta.vec} should be 2 values (or 1
#'     value to be split evenly). Otherwise, \code{epsilon.vec} should be 2 values
#'     (or 1 to be split evenly) and \code{delta.vec} should be 1.
#'
#' @family ApproximateDPmeasures
#' @importFrom VGAM rlaplace
#' @importFrom stats qnorm
#'
#' @examples
#'set.seed(1)
# x=rnorm(500,70,5)
# dp_confidence_interval(x=x,epsilon.vec=1,delta.vec=0.001,bound.mean=1e5)
#'
#'
#'
#' @references
#' Karwa, V. and Vadhan, S. (2017). Finite sample differentially private confidence intervals.
#'
#'
dp_confidence_interval=function(x,epsilon.vec,delta.vec=0,alphas=0.05,san.point=NA,
                                bounds.sd=c(2^(-15),2^15),x.sd=NA,bound.mean=NA,
                                return.point.sd=FALSE){

  num.priv.steps=ifelse(is.na(x.sd)==TRUE,1,0)
  num.priv.steps=num.priv.steps+
    ifelse(is.na(san.point)==TRUE,2,0)
  num.priv.steps=max(num.priv.steps,1)
  # if epsilon, delta, or alpha are one value. Split the value evenly across the steps.
  if(length(epsilon.vec)==1){
    epsilon.vec=rep(epsilon.vec/num.priv.steps,num.priv.steps)
  }
  if(length(delta.vec)==1){
    delta.vec=rep(delta.vec/(num.priv.steps-1),max(num.priv.steps-1,1))
  }
  if(length(alphas)==1){
    alphas=rep(alphas/3,3)
  }
  stopifnot(sum(alphas)>0&sum(alphas)<1) #significance level must be in (0,1)

  n=length(x)
  if(is.na(x.sd)==TRUE){ #if no standard deviation provided
    san.sd=dp_estimate_sd(x=x,epsilon=epsilon.vec[1],delta=delta.vec[1],bounds.sd=bounds.sd)
  }else{
    san.sd=x.sd
    if(is.null(attr(san.sd,"priv.cost"))==TRUE){
      epsilon.vec=c(0,epsilon.vec)
      delta.vec=c(0,delta.vec)
    }else{
      epsilon.vec=c(attr(san.sd,"priv.cost")[1],epsilon.vec)
      delta.vec=c(attr(san.sd,"priv.cost")[2],delta.vec)
    }
  }


  if(is.na(san.point)==TRUE){
    san.range=dp_range(x=x,sd=san.sd,epsilon=epsilon.vec[2],delta=delta.vec[2],
                       bound.mean=bound.mean,range.prob=alphas[3])

    #trim values to be within the dp range
    x[x<san.range[1]]=san.range[1]
    x[x>san.range[2]]=san.range[2]

    scale.param=(san.range[2]-san.range[1])/(epsilon.vec[3]*n) #scale param for laplace noise
    san.point=mean(x)+VGAM::rlaplace(1,0,scale.param) #sanitized point estimate
  }else{
    if(is.null(attr(san.point,"priv.cost"))==TRUE){
      epsilon.vec=c(0,epsilon.vec)
      delta.vec=c(0,delta.vec)
    }else{
      epsilon.vec=c(attr(san.point,"priv.cost")[1],epsilon.vec)
      delta.vec=c(attr(san.point,"priv.cost")[2],delta.vec)
    }
  }

  #half the width of the confidence interval
  half.width=((san.sd/sqrt(n))*stats::qnorm(1-(alphas[1]/2)))+
    (scale.param*base::log(1/alphas[2]))
  #upper and lower bounds of CI
  san.lower=san.point-half.width
  san.upper=san.point+half.width
  output=c("san.lower"=san.lower,"san.upper"=san.upper)
  attr(output,"significance.level")=sum(alphas)
  if(return.point.sd==TRUE){
    output=list("san.CI"=output,"san.coef"=san.point,"san.sd"=san.sd)
  }
  attr(output,"priv.cost")=c("epsilon"=sum(epsilon.vec),"delta"=sum(delta.vec))
  return(output)
}

