#' Simulate Responses Values Given a Model and New Data
#'
#' @param mod glm model object
#' @param newdata a data.frame of the new data of predictors to generate responses from
#' @param predictor.formula a formula object that is the formula used in the \code{mod} object.
#'    If \code{NULL}, the formula is extracted from \code{mod}.
#' @param nsim a integer for the number of responses to simulate from the \code{newdata}
#' @param cov.mat the covariance matrix for the estimated model coefficients. If NULL
#'    (default), then \code{stats::vcov(mod)} is used.
#' @returns a vector (or matrix) with \code{nsim} columns and \code{nrows(newdata)} rows
#'    of response randomly generated from the \code{mod} using \code{newdata} as the
#'    predictors. Each column is a new simulation.
#' @importFrom stats as.formula coefficients vcov model.matrix
#' @importFrom mvtnorm rmvnorm
#'
#' @export
#'
simulate_response_glm=function(mod,newdata,predictor.formula=NULL,
                               nsim=1,cov.mat=NULL){

  #get predictor formula from mod if no predictor.formula given
  if(is.null(predictor.formula)==TRUE){
    predictor.formula=stats::as.formula(paste0("~",as.character(stats::as.formula(mod$formula))[3]))
  }else if(is.character(predictor.formula)==TRUE){
    predictor.formula=stats::as.formula(predictor.formula)
  }

  #get predictor variables and data
  pred.vars=attr(mod$terms,"term.labels")

  pred.newdata=data.frame(newdata[,pred.vars])
  colnames(pred.newdata)=pred.vars

  #get estimated model coefficients
  mod.coefs=stats::coefficients(mod)
  has.coef.name=names(mod.coefs)[!is.na(mod.coefs)] #only columns where coefficient!=NA


  #model predictor matrix
  modMat=stats::model.matrix(predictor.formula, pred.newdata)
  modMat=modMat[,colnames(modMat)%in%has.coef.name] #remove columns with no coefficient
  modMat.cnames=colnames(modMat)
  not.in.modMat=has.coef.name[!(has.coef.name %in% modMat.cnames)]
  if(length(not.in.modMat)>0){ #if some variables don't appear
    zeros.mat=matrix(0,ncol=length(not.in.modMat),nrow=nrow(modMat))
    modMat=cbind(modMat,zeros.mat)
    colnames(modMat)=c(modMat.cnames,not.in.modMat)
    modMat=modMat[,has.coef.name]
  }
  #if no cov.mat provided use the one from mod
  if(is.null(cov.mat)==TRUE){
    cov.mat=stats::vcov(mod)
    cov.mat=cov.mat[rownames(cov.mat)%in%has.coef.name,colnames(cov.mat)%in%has.coef.name] #remove columns and rows for NA coefs
    #NOTE: if models are generated from san.X
    # cov.mat=(mean(resid^2)/(n-num.predictors-1))solve(t(san.X)%*%san.X/n)
  }

  #generate random coefficients
  coef.rv=mvtnorm::rmvnorm(nsim, mod.coefs[names(mod.coefs)%in% has.coef.name], cov.mat)



  #get model value before inverse link function (still need to get response)
  sim.model=tcrossprod(modMat,coef.rv)
  mod.family=mod$family
  sim.response=mod.family$linkinv(sim.model) #get response

  #if there is more than 1 simulation label the simulations
  if(nsim>1){colnames(sim.response)=paste("sim",seq(1,nsim),sep="_")
  }else{
    as.numeric(sim.response)
  }


  return(sim.response)
}
