#' Internal-use simulate B estimated model coefficients responses from confidential model
#'
#'@param conf.model a glm or lm object of the model fitted using the confidential data
#'@param confidential.data a data.frame of the confidential model variables
#'@param predictor.vars a vector of strings, indices, or logicals to indicate the
#'      covariate and treatment variable columns of the \code{confidential.data}.
#'@param response.var a string, index, or vector of logicals to indicate the response
#'    variable of the model
#'@param num.iters the number of iterations for generating proxy responses and
#'    fitting the model.
#'@returns a matrix with columns for the intercept, treatment effect, and covariate
#'    coefficients. Rows are each iteration.
#'
#'@importFrom stats simulate
#'
#'@keywords internal
#'@noRd
iter_hybrid=function(conf.model,
                     confidential.data,
                     num.iters,
                     predictor.vars,
                     response.var){
  yb=stats::simulate(conf.model,nsim=num.iters,type="response") #simulate proxy responses

  #########
  #FUNCTION:attach proxy response to predictor data, fit model, return coefficients
  get_betas=function(coly,
                     conf.predictors=confidential.data[,predictor.vars],
                     conf.mod=conf.model,
                     y.var=response.var
  ){
    iter.data=cbind(conf.predictors,coly)
    colnames(iter.data)=c(colnames(conf.predictors),y.var)
    modb=stats::glm(formula=conf.mod$formula,family = conf.mod$family,data=iter.data)
    return(modb$coefficients)
  }
  #######

  betas=t(sapply(1:num.iters,function(x)get_betas(yb[,x])))
  #use get_betas for each iteration.
  return(betas[,!is.na(conf.model$coefficients)])
}


#' Internal-use to simulate B proxy responses and sanitize the coefficient covariance matrix from confidential model
#'
#'@param conf.model a glm or lm object of the model fitted using the confidential data
#'@param confidential.data a data.frame of the confidential model variables
#'@param model.vars a vector of strings for model variables with the first being the response.variable
#'@param response.var a string, index, or vector of logicals to indicate the response
#'    variable of the model
#'@param num.iters the number of iterations for generating proxy responses and
#'    fitting the model.
#'@returns a list with "iter.beta" and "cov.mat". "iter.beta" is a matrix with columns
#'    for the estimated model coefficients and a row for each iteration. "cov.mat" is
#'    the sanitized covariance matrix for the estimated model coefficients that is made
#'    from a sanitized MSE and the synthetic covariate and treatment data.
#'
#'@importFrom stats simulate
#'
#'@keywords internal
#'@noRd
dp_iter_hybrid=function(conf.model,
                        confidential.data,
                        synth.data,
                        num.iters,
                        model.vars,
                        response.var,
                        mse.epsilon,
                        mse.delta,
                        mse.bd.sd,
                        mse.bd.mean=100){
  yb=simulate_response_glm(conf.model,newdata=synth.data[,model.vars[-1]],nsim=num.iters)

  predictor.formula=stats::as.formula(
    paste0("~",as.character(
      as.formula(conf.model$formula))[3]))
  pred.newdata=data.frame(synth.data[,model.vars[-1]])
  colnames(pred.newdata)=model.vars[-1]
  modMat=stats::model.matrix(predictor.formula, pred.newdata)

  #get estimated model coefficients
  mod.coefs=stats::coefficients(conf.model)
  has.coef.name=names(mod.coefs)[!is.na(mod.coefs)] #only columns where coefficient!=NA
  num.coefs=length(has.coef.name)

  modMat=modMat[,colnames(modMat)%in%has.coef.name]
  #########
  #FUNCTION:attach proxy response to predictor data, fit model, return coefficients
  get_betas_and_residuals=function(coly,
                                   synth.predictors=synth.data[,model.vars[-1]],
                                   conf.mod=conf.model,
                                   y.var=response.var,
                                   return.invxtx=FALSE,
                                   treat.se.normal=FALSE){

    #bind response to predictors
    iter.data=cbind(synth.predictors,coly)
    colnames(iter.data)=c(colnames(synth.predictors),y.var)
    #fit model
    modb=stats::glm(formula=conf.mod$formula,family = conf.mod$family,data=iter.data)
    if(as.numeric(return.invxtx)==1){
      inv.xtx=stats::vcov(modb)*summary(modb)$df.residuals/sum(modb$residuals^2)
    }else{
      inv.xtx=NULL
    }
    #return coefficient estimates and residuals
    return(list("betas"=modb$coefficients,"residuals"=modb$residuals,"inv.xtx"=inv.xtx))
  }
  #######

  #parallelize apply to each iteration
  iter.out=parallel::mclapply(1:num.iters,
                              function(x)get_betas_and_residuals(yb[,x],return.invxtx = x))
  #combine coefficients into matrix
  #column for each coefficient (including intercept), row for each iteration
  betas=t(sapply(1:num.iters,function(idx)iter.out[[idx]]$betas))

  if(treat.se.normal==FALSE){
  #stack residuals into a vector
  residuals=unlist(lapply(1:num.iters,function(idx)iter.out[[idx]]$residuals))

  #get sanitized estimate of residual standard deviation
  san.mse=(dp_estimate_sd(residuals,
                          epsilon=mse.epsilon,delta=mse.delta,
                          bounds.sd = mse.bd.sd))^2
  }else{
    se=unlist(lapply(1:num.iters,function(idx)base::sqrt(stats:var(iter.out[[idx]]$residuals))/length(iter.out[[idx]]$residuals)))
    ci.out=dp_confidence_interval(se,epsilon.vec=mse.epsilon,delta.vec=mse.delta,alphas=0.05,san.point=NA,
                                    bounds.sd=mse.bd.sd,x.sd=NA,bound.mean=mse.bd.mean,
                                    return.point.sd=TRUE)
    san.mse=ci.out[[2]]
  }
  nr=nrow(confidential.data) #number of observations

  if(sum(colSums(modMat==0)==nr)>0){
    warning("some columns in modMat have only 0 values")
  }
  #warning(paste("colnames modMat",paste0(colnames(modMat),collapse=", "),"removed colnames are",paste0(names(mod.coefs)[is.na(mod.coefs)],collapse=", ")))
  #covariance matrix of coefficients is sigma^2
  inv.xtx=iter.out[[1]]$inv.xtx

  cov.mat=(san.mse/(nr-num.coefs))*nr*(inv.xtx)#t(modMat)%*%modMat/nr)
  #colnames(cov.mat)=colnames(modMat)
  #rownames(cov.mat)=colnames(modMat)


  return(list("iter.betas"=betas[,colnames(betas)%in%has.coef.name],"cov.mat"=cov.mat,"san.mse"=san.mse))
}
