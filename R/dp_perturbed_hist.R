#' Generate sanitized proportions that satisfy (epsilon,delta)-DP
#'
#' @param hist.df a data.frame of the counts for the number of occurrences of a category (or combination of categories) which is column \code{hist.df$Freq}vector of continuous values to be split into bins
#' @param epsilon a positive privacy parameter value
#' @param delta a non-negative value for a privacy parameter. If \code{nrow(hist.df)<=}
#' @return a data.frame which is the hist.df with column \code{hist.df$san.prop} which is a sanitized proportion of the data that is in a bin which is generated by a Laplace mechanism satisfying (\code{epsilon},\code{delta})-differential privacy.
#'
#' @examples
#' freq.df=multivariate_histogram(mtcars[,2:4],continuous.vars=c("disp","hp"),num.bin=3)
#' dp_perturbed_hist(freq.df,0.9,delta=0.01)
#'
#' @references
#' When \code{delta==0} or \code{nrow(hist.df)<=2/delta}, perturbed histogram algorithm is from:
#' Karwa, V. and Vadhan, S. (2017). Finite sample differentially private confidence intervals.
#'
#' When \code{delta>0} and \code{nrow(hist.df)>2/delta}, perturbed histogram algorithm is from:
#' Bun, M., Nissim, K., and Stemmer, U. (2016). Simultaneous private learning of multiple
#' concepts. In Proceedings of the 2016 ACM Conference on Innovations in Theoretical Computer Science,
#' ITCS ’16, page 369–380, New York, NY, USA. Association for Computing Machinery.
#'
#'
#' @family syntheticData
#' @importFrom VGAM rlaplace
#' @noRd


dp_perturbed_hist<-function(hist.df,epsilon,delta=0,possible.combos=NULL){

 # warn.message=NULL
  hist.df=hist.df[hist.df$Freq>0,,drop=F] #remove bins with 0 counts
  nobs=base::sum(hist.df$Freq,na.rm=T) #number of observations
  if(nobs<10){
    message(paste("The number of histogram observations is less then 10:",nobs,"\n"))
  }
  #if using delta>0, and number of bins is high enough
  gen_unrealized=function(missing.combos,threshold,sc.param,rseed=NA){
    if(is.na(rseed)==FALSE){
      set.seed(rseed)
    }
    # count.threshold=nobs*threshold
    # if((1/epsilon)<1){
    #   message("epsilon is such that the discrete laplace mechanism can be used")
    #   p.above=extraDistr::::pdlaplace(count.threshold,mu=0, scale=(1/epsilon),lower.tail=F)
    # }else{
    #   message("epsilon is such that the discrete laplace mechanism can not be used")
    #   p.above=extraDistr::::pdlaplace(count.threshold,mu=0, scale=(1/epsilon),lower.tail=F)
    # }
    p.above=VGAM::plaplace(threshold,scale=sc.param)
    rcount.above.threshold=rbinom(1,missing.combos,p.above)
    #print(rcount.above.threshold)
    unif.prob1=stats::runif(rcount.above.threshold/2,1-p.above,1)
    unif.prob2=stats::runif(rcount.above.threshold/2,1-p.above,1)
    san.prop1=VGAM::qlaplace(unif.prob1,0,sc.param)
    san.prop2=VGAM::qlaplace(unif.prob2,0,sc.param)
    return(c(san.prop1,san.prop2))
  }


  if(is.null(possible.combos)==TRUE){
    #message("No possible.combos provided. It is assumed all potential combinations of variables are represented in hist.df.")
    #if(possible.combos==base::nrow(hist.df)){
    if((delta>0)&(base::nrow(hist.df)>(2/delta))){ #use Bun et al. 2016
      threshold=((2*base::log(2/delta))/(epsilon/nobs))+(1/nobs)
    }else{ #pure-DP and low number of bins don't use a threshold
      threshold=0
    }
    sanprop.zero.to.add=0
    used.sanprop.add=F
  }else{
    used.sanprop.add=T
    missing.combos=possible.combos-nrow(hist.df)
    if((delta>0)&(possible.combos>(2/delta))){ #use Bun et al. 2016
      threshold=((2*base::log(2/delta))/(epsilon/nobs))+(1/nobs)
    }else{ #pure-DP and low number of bins don't use a threshold
      threshold=0
    }
    sc.param=2/(nobs*epsilon)
    sanprop.zero.to.add=gen_unrealized(missing.combos = missing.combos,
                                       threshold=threshold,sc.param=sc.param)
    # if(missing.combos<=(2^14)){
    #   print("missing combos <=2^14")
    #   bins.zero.san.prop= VGAM::rlaplace(missing.combos,0,sc.param)
    # }else{
    #   reps.samp=missing.combos%/%(2^14)
    #   bins.zero.san.prop=VGAM::rlaplace(missing.combos%%(2^14),0,sc.param)
    #   bins.zero.san.prop[bins.zero.san.prop>threshold]
    #   for(i in seq(1,reps.samp)){
    #     temp.san.prop=VGAM::rlaplace(2^14,0,sc.param)
    #     bins.zero.san.prop=c(bins.zero.san.prop,temp.san.prop[temp.san.prop>threshold])
    #   }
    #}
    #print(summary(sanprop.zero.to.add))
    if(sum(sanprop.zero.to.add<threshold)>0){
      stop("Sanitized proportion generation has made an unknown error. There are produced sanitized proportions that are below the threshold.")
    }
   # sanprop.zero.to.add=bins.zero.san.prop[bins.zero.san.prop>threshold]
    #message(paste("There are",length(sanprop.zero.to.add)," zero bins to add."))
  }
  num.bins=base::nrow(hist.df) #number of bins
  hist.df$san.prop=hist.df$Freq/nobs

  bins.pos=hist.df$san.prop>0


  #add laplace noise
  hist.df$san.prop[bins.pos]=(hist.df$san.prop[bins.pos]+
                                VGAM::rlaplace(sum(bins.pos),0,2/(nobs*epsilon)))
  ### if use discrete laplace  extraDistr::rdlaplace(sum(bins.pos),0,2/(nobs*epsilon)))

  if(sum(hist.df$san.prop>threshold)+sum(sanprop.zero.to.add)==0){
    stop("After noise added, there are no sanitized proportions above the threshold. This is not good. Try again?")
    #min.san.prop=min(c(hist.df$san.prop,sanprop.zero.to.add))
    #hist.df$san.prop=hist.df$san.prop+abs(min.san.prop)
  }else{
    hist.df$san.prop=base::pmax(hist.df$san.prop,threshold)
  }
  #if(is.null(warn.message)==FALSE){
  #  warning(warn.message)
  #}
  total=(sum(hist.df$san.prop)+sum(sanprop.zero.to.add))
  hist.df$san.prop=hist.df$san.prop/total #normalize
  if(used.sanprop.add==TRUE){
    total=(sum(hist.df$san.prop)+sum(sanprop.zero.to.add))
    hist.df$san.prop=hist.df$san.prop/total #normalize
    sanprop.zero.to.add=sanprop.zero.to.add/total
    return(list("hist.df"=hist.df,"san.prop"=sanprop.zero.to.add))
  }else{
    hist.df$san.prop=hist.df$san.prop/sum(hist.df$san.prop)
    return(hist.df)
  }
}

