#' P-Values from Approximate DP
#' @description Algorithm 1 by Kazan et al.(2023) and uses functions from binomialDP package (Tran,2020).
#' @param x is confidential data to be tested
#' @param test is a function that tests data \code{x} and returns a p-value
#' @param alpha is the significance level of the tests. Default is 0.05
#' @param epsilon is the positive privacy parameter
#' @param delta is a non-negative privacy parameter. Default is \code{delta=0}.
#' @param pvalue.partitions is the number of subsets of \code{x} to use to produce p-values.
#'    If \code{pvalue.partitions} is NA, then \code{sub.sizes} cannot be NULL and
#'    \code{pvalue.partitions} will be set to \code{length(sub.sizes)}
#' @param sub.sizes a vector of the size for each subset of the partition of \code{x}.
#'    If \code{pvalue.partitions} is not NA and \code{sub.sizes} is not NULL,
#'    \code{pvalue.partitions} will be used to override \code{sub.sizes}.
#' @param ... any additional parameters are inputted into the \code{test} function.
#' @returns a list of the sanitized p-value and the test statistic \code{z}.
#'
#' @examples
#' set.seed(1)
#' test.func=function(x){pnorm(mean(x),lower.tail=FALSE)}
#' dp_pvalue(x=rnorm(500,1.2,1),epsilon=1,test=test.func,alpha=0.05,pvalue.partitions=12)
#'
#' @family ApproximateDPmeasures
#'
#' @references \code{dp_pvalue} is Algorithm 1 in
#'     Kazan, Z., Shi, K., Groce, A.,and Bray, A. P. (2023). The test of tests:
#'     A framework for differentially private hypothesis testing. In Krause, A.,
#'     Brunskill, E., Cho, K., Engelhardt, B., Sabato, S.,and Scarlett, J., editors,
#'     Proceedings of the 40th International Conference on Machine Learning,
#'     volume 202 of Proceedings of Machine Learning Research, pages 16131–16151. PMLR.
#'
#'     \code{dp_pvalue} used binomialDP package from GitHub.
#'     Tran, T. (2020). Differentiallly private inference for binomial data.
#'     \url{https://github.com/tranntran/binomialDP.}
#'
#'     binomialDP uses results from:
#'     Awan, J. and Slavkovi´c, A. (2018). Differentially private uniformly most
#'        powerful tests for binomial data. Advances in Neural Information
#'        Processing Systems, 31:4208–4218.
#'
#'
#' @export

dp_pvalue=function(x,test,epsilon,delta=0,alpha=0.05,
                   pvalue.partitions=NA,sub.sizes=NULL,...){

  ## The necessary functions from binomialDP are copied below
  ## code from https://github.com/tranntran/binomialDP/tree/master/R
  ## SHA: bd338e518b405e12312586c8b96259f7765e639d
    pvalRight <- function(Z, size, theta, b, q){
      reps = length(Z)
      pval = rep(0,reps)
      values = seq(0,size)

      B = stats::dbinom(values, size = size, prob = theta)

      for(r in 1:reps){
        F = ptulap(t = values-Z[r], m = 0, b, q)
        pval[r]=t(F)%*%B
      }
      return(pval)
    }
    rtulap <- function (n, m = 0, b = 0, q = 0) {
      if(q>=0){
        alpha = .95
        lcut = q/2
        rcut = q/2

        # Approximates M such that given n Bernoulli trials with success rate prob,
        # is such that alpha of the times, there are at least n successes among M. n
        # - the number of trials that we want more than success of prob - Bernoulli
        # probability success of each iid trial alpha - probability that among the M
        # trials, at least n successes.
        approx.trials <- function (n, prob = 1, alpha = 0) {
          # Solve a quadratic form for this:
          a = prob^2
          b = -((2 * n * prob) + ((stats::qnorm(alpha)^2) * prob * (1 - prob)))
          c = n^2
          return ((-b + base::sqrt(b^2 - (4 * a * c))) / (2 * a))
        }

        # Calculate actual amount needed
        q = lcut + rcut
        n2 = approx.trials(n, prob=(1 - q), alpha=alpha)

        # Sample from the original Tulambda distribution
        geos1 = stats::rgeom(n2, prob=(1 - b))
        geos2 = stats::rgeom(n2, prob=(1 - b))
        unifs = stats::runif(n2, min=(-1/2), max=(1/2))
        samples = m + geos1 - geos2 + unifs

        # Cut the tails based on the untampered CDF (ie no cuts)
        probs = ptulap(samples, m = m, b = b)
        is.mid = (lcut <= probs) & (probs <= (1 - rcut))

        # Abuse the NA property of R wrt arithmetics
        mids = samples[is.mid]
        while ({len = length(mids); len} < n) {
          diff = n - len
          mids = c(mids, rtulap(diff, m=m, b=b,q=q))
        }
        return (mids[1:n])
      }
      geos1 = stats::rgeom(n2, prob=(1 - b))
      geos2 = stats::rgeom(n2, prob=(1 - b))
      unifs = stats::runif(n2, min=(-1/2), max=(1/2))
      samples = m + geos1 - geos2 + unifs
      return(samples)
    }

    ptulap <- function (t, m = 0, b = 0, q = 0) {
      lcut = q/2
      rcut = q/2
      # Normalize
      t = t - m

      # Split the positive and negsative t calculations, and factor out stuff
      r = round(t)
      g = -base::log(b)
      l = base::log(1 + b)
      k = 1 - b
      negs = base::exp((r * g) - l + base::log(b + ((t - r + (1/2)) * k)))
      poss = 1 - base::exp((r * (-g)) - l + base::log(b + ((r - t + (1/2)) * k)))

      # Check for infinities
      negs[is.infinite(negs)] <- 0
      poss[is.infinite(poss)] <- 0

      # Truncate wrt the indicator on t's positivity
      is.leq0 = t <= 0
      trunc = (is.leq0 * negs) + ((1 - is.leq0) * poss)

      # Handle the cut adjustment and scaling
      q = lcut + rcut
      is.mid = (lcut <= trunc) & (trunc <= (1 - rcut))
      is.rhs = (1 - rcut) < trunc
      return (((trunc - lcut) / (1 - q)) * is.mid + is.rhs)
    }


  n=length(x)

  # if m is given, it determines subset sizes
  # otherwise m is the length of sub.sizes
  if(base::is.na(pvalue.partitions)==FALSE){
    stopifnot(is.numeric(pvalue.partitions)&&pvalue.partitions>0)
    remainder=n%%pvalue.partitions
    sub.sizes=rep(n%/%pvalue.partitions,pvalue.partitions)
    sub.sizes[1:remainder]=sub.sizes[1:remainder]+1
  }else{ #check subsample sizes sum to the number of entries.
    stopifnot(sum(sub.sizes)==n)
    pvalue.partitions=length(sub.sizes)
  }

  #each row of x gets a value to indicate which subsample it is in.
  partition.idx=sample(rep(1:pvalue.partitions,sub.sizes))

  #get the p-value for each subsample xj using test function. If test does not work,
  # get a randomly generated uniform random variable instead.
  get_pj=function(xj,tau=test,...){
    pj=try(tau(xj,...),silent=T)
    if(is.numeric(pj)==F){
      pj=stats::runif(1,0,1)
    }
    return(pj)
  }
  ### End get_pj function

  #get confidential p-values for each subsample
  sub.conf.pvalue=sapply(1:pvalue.partitions,
                         function(idx)get_pj(x[partition.idx==idx],tau=test))
  #count the number of pvalues below significance level
  num.sig.pvalues=sum(sub.conf.pvalue<alpha)
  #tulap parameters to meet (epsilon,delta)-DP
  b = base::exp(-epsilon)
  q = (2*delta*b)/(1-b+(2*delta*b))
  z=rtulap(1,m=sum(sub.conf.pvalue<=alpha),b=b,q=q)
  san.pvalue=pvalRight(z,size=pvalue.partitions,theta=alpha,b=b,q=q)
  output=list("san.pvalue"=san.pvalue,"z"=z)
  attr(output,"priv.cost")=c("epsilon"=epsilon,"delta"=delta)
  return(list("san.pvalue"=san.pvalue,"san.z.value"=z))
}
