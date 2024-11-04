#' Randomized Response changes a binary (0,1) variable with some probability
#' @param binary.vec vector of 0,1 values (or TRUE FALSE values)
#' @param epsilon a positive numeric privacy parameter. If supplied,
#'    the randomized output satifies \code{epsilon}-differential privacy
#' @param prob.flip (optional) a numeric value between 0 and 1. If \code{epsilon=NA},
#'    then \code{prob.flip} cannot be \code{NA} and it determines the probability
#'    any binary entry is flipped.
#' @param rseed a numeric value if supplied sets the random seed
#' @returns randomized binary vector using randomized response.
#'
#' @examples
#' set.seed(1)
#' x=c(rep(0,20),rep(1,20))
#' randomized_response(x,epsilon=1)
#'
#'
#' @family DPmechanisms
#'
#' @export
randomized_response<-function(binary.vec,epsilon=NA,prob.flip=NA,rseed=NA){
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  if(is.logical(binary.vec)){
    binary.vec=as.numeric(binary.vec)
  }
  n=length(binary.vec)
  stopifnot(is.numeric(binary.vec)&&sum(binary.vec==1)+sum(binary.vec==0)==n)
  if(is.na(epsilon)&is.na(prob.flip)){stop("Need epsilon and prob.flip.")}
  if(is.na(epsilon)==FALSE){
    prob.flip=(1-((base::exp(epsilon)-1)/(base::exp(epsilon)+1)))/2
  }
  changes=sample(c(0,1),n,replace = T,
                       prob = c(1-prob.flip,prob.flip))
  return((binary.vec+changes)%%2)
}


#' Gaussian Mechanism
#' @description Gaussian Mechanism infuses noise to satisfy approximate differential privacy.
#' @param x vector of numeric confidential values
#' @param epsilon a positive value for the privacy cost.
#' @param delta a numeric positive value for the privacy cost
#' @param f a function applied to \code{x} which will have its result protected
#'    with Gaussian noise. If \code{f=NULL}, then it is the identity function (i.e. f(x)=x).
#' @param sensitivity the L2 sensitivity of \code{f}
#' @param rseed a numeric value if supplied sets the random seed
#' @returns randomized output that is generated with Gaussian noise to satisfy
#'    approximate differential privacy.
#'
#' @examples
#' set.seed(1)
#' x=rbinom(50,20,0.5)/20
#' gaussian_mechanism(x,epsilon=1,delta=0.001,sensitivity=1)
#'
#'
#' @family DPmechanisms
#'
#' @importFrom stats rnorm
#' @export
gaussian_mechanism<-function(x,epsilon,delta,f=NULL,sensitivity=NULL,rseed=NA){
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  stopifnot(is.numeric(sensitivity))
  if(is.null(f)==FALSE){
    x=f(x)
  }
  noise.sd=sqrt(2*base::log(1.25/delta)*sensitivity/epsilon)
  return(x+stats::rnorm(length(x),mean=0,sd=noise.sd))
}

#' Laplace Mechanism
#' @description Laplace Mechanism infuses noise to satisfy pure differential privacy.
#' @param x vector of numeric confidential values
#' @param epsilon a positive value for the privacy cost.
#' @param f a function applied to \code{x} which will have its result protected
#'    with Gaussian noise. If \code{f=NULL}, then it is the identity function (i.e. f(x)=x).
#' @param sensitivity the L2 sensitivity of \code{f}
#' @param rseed a numeric value if supplied sets the random seed
#' @returns randomized output that is generated with Laplace noise to satisfy
#'    pure differential privacy.
#'
#' @examples
#' set.seed(1)
#' x=rbinom(50,20,0.5)/20
#' laplace_mechanism(mtcars$mpg,epsilon=1,sensitivity=1)
#'
#'
#' @family DPmechanisms
#'
#' @importFrom VGAM rlaplace
#' @export
laplace_mechanism<-function(x,epsilon,f=NULL,sensitivity=NULL,rseed=NA){
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  stopifnot(is.numeric(sensitivity))
  if(is.null(f)==FALSE){
    x=f(x)
  }
  scale.param=sensitivity/epsilon
  if(scale.param<=0){
    stop(paste("The scale param is not positive,",scale.param," Using the sensivity:",sensitivity, " and epsilon",epsilon))
  }
  return(x+VGAM::rlaplace(length(x),0,scale.param))
}
