#' Approximate DP Synthetic Data and Treatment Effect Table
#'
#' @inheritParams dp_model_summary
#' @inheritParams hybrid_synth
#' @param synth.epsilon privacy parameter for synthetic covariate and treatment data
#' @param synth.delta privacy parameter for synthetic covariate and treatment data. Default 0
#' @inheritParams synthdata_perturb_mvhist
#' @param treat.ppart.list list of p-value partitions for each treatment effect
#' @param return.treatment.pvalue indicator to return pvalues for treatment effects
#' @param return.treatment.CI indicator to return confidence intervals for treatment effects
#'      (last \code{alphas.list} value(s) will be used for confidence level)
#' @return synthetic data set (and possible treatment pvalue, CI and confidential summary table.)
#'
#' @details To become the example:
#' tempdf=mtcars[,1:6]
#' tempdf$am=as.factor(mtcars$am)
#' tempdf$cyl=as.factor(mtcars$cyl)
#' tform="mpg~am+cyl+hp+wt"
#' tcont=c("hp","wt")
#' dp_synthdata(tform,tempdf,synth.epsilon = 1,continuous.vars = tcont,
#' num.bin=3,epsilon.list=1,bd.mean.list=100,num.iters=300,use.san.residerror = T,
#' treat.ppart.list = 5)
#'
#' @importFrom parallel mclapply
#'
#'
#'
#' @export
dp_synthdata=function(formula,
                      confidential.data,
                      synth.data=NULL,
                      synth.epsilon=NA,
                      synth.delta=0,
                      continuous.vars=NULL,
                      num.bin=NULL,
                      bin.param=NA,
                      assign.type="simple",
                      blocks=NULL,
                      clusters=NULL,
                      within.blocks=TRUE,
                      epsilon.list,
                      delta.list=0,
                      bd.sd.list=c(2^(-15),2^(15)),
                      bd.mean.list,
                      alphas.list=0.05,
                      num.iters=1e4,
                      treatment.var=NULL,
                      rseed=NA,
                      return.time=FALSE,
                      return.confidential.table=FALSE,
                      return.san.summary=TRUE,
                      use.san.residerror=TRUE,
                      return.treatment.pvalue=TRUE,
                      return.treatment.CI=TRUE,
                      treat.ppart.list=NULL,
                      ...){
  start.time=proc.time()
  time.values=NULL
  #print("Inside dp_synthdata")

  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  if(is.null(synth.data)==TRUE){
    #if not synthetic data supplied, a positive privacy cost needs to be supplied
    stopifnot(is.numeric(synth.epsilon)&&synth.epsilon>0)
  }
  psf=parse_formula(formula=formula,confidential.data = confidential.data,
                    treatment.var = treatment.var)
  formula=psf[[1]]
  model.vars=psf[[2]]
  response.var=psf[[3]]
  synth.vars=psf[[4]]
  treatment.var=psf[[5]]

  synth.vars=c(model.vars[-1],blocks,clusters)
  synth.vars=synth.vars[!duplicated(synth.vars)]


  confidential.data=confidential.data[,colnames(confidential.data)%in%c(response.var,synth.vars)]

  model.fit.start=proc.time()
  #### Step 0: In paper
  #fit the confidential data to the model and get estimated coefficient and standard error
  conf.model<-stats::glm(formula=formula,family="gaussian",data=confidential.data,...)
  ####
  if(return.time==TRUE){
    time.values=c(time.values,
                  list("conf.model.fit"=(proc.time()-model.fit.start)[[3]]))
  }

  # if synthetic data is supplied, check it is data.frame w/ all the model.vars
  if((is.null(synth.data)==FALSE)&&
     (is.data.frame(synth.data)==TRUE)&&
     (sum(c(model.vars[model.vars!=response.var])%in%colnames(synth.data))==length(model.vars)-1)){
    synth.data=synth.data
  }else{ #if synth.data not supplied, check it closely
    synth.data<-handle_synthetic_data(confidential.data=confidential.data,
                                      synth.data=synth.data,
                                      treatment.var=treatment.var,
                                      synth.vars=synth.vars,
                                      model.vars=model.vars,
                                      epsilon=synth.epsilon,
                                      delta=synth.delta,
                                      continuous.vars=continuous.vars,
                                      num.bin=num.bin,
                                      bin.param=bin.param,
                                      add.cont.variation=TRUE,
                                      assign.type=assign.type,
                                      blocks=blocks,
                                      clusters=clusters,
                                      within.blocks=within.blocks,
                                      return.time=TRUE)
    time.values=c(time.values,
                  list("synthetic.T.X.time"=attr(synth.data,"comp.time")))
  }
  num.coefs=sum(!is.na(conf.model$coefficients))
  #get model matrix
  predictor.formula=stats::as.formula(paste0("~",as.character(formula)[3]))
  pred.newdata=data.frame(synth.data[,model.vars[-1]])
  colnames(pred.newdata)=model.vars[-1]
  modMat=stats::model.matrix(predictor.formula, pred.newdata)

  #print("after make modMat")

  iter.proxy.fit.start=proc.time()
  if(use.san.residerror==TRUE){
    #print("insdie use.san.residerror==TRUE")
    ### pulls off last element of list to be used for mse
    get_sanmse_params=function(param.list){
      if(is.list(param.list)==FALSE){
        param.list=list(param.list)
      }
      mse.param=param.list[[length(param.list)]]
      if(length(param.list)>1){
        param.list=param.list[seq(1,length(param.list)-1)]
      }
      list(mse.param,param.list)
    }


    #get epsilon, delta, and bounds.sd parameters for sanitizing mse
    tempout=get_sanmse_params(epsilon.list)
    mse.epsilon=tempout[[1]]
    epsilon.list=tempout[[2]]
    tempout=get_sanmse_params(delta.list)
    mse.delta=tempout[[1]]
    if(length(mse.delta)>1){
      mse.delta=mse.delta[1]
    }
    delta.list=tempout[[2]]
    tempout=get_sanmse_params(bd.sd.list)
    mse.bd.sd=tempout[[1]]
    bd.sd.list=tempout[[2]]

    warning("before dp_iter_hybrid")
    list2env(dp_iter_hybrid(conf.model=conf.model,
                            confidential.data=confidential.data,
                            synth.data=synth.data,
                            num.iters=num.iters,
                            model.vars=model.vars,
                            response.var=response.var,
                            mse.epsilon=mse.epsilon,
                            mse.delta=mse.delta,
                            mse.bd.sd=mse.bd.sd),envir=globalenv())


    coef.names=colnames(cov.mat)
  }else{#don't use residerror
    #print("inside use_san_residerror==FALSE")
    iter.betas=iter_hybrid(conf.model=conf.model,
                         confidential.data=confidential.data,
                         num.iters=num.iters,
                         predictor.vars=model.vars[-1],
                         response.var=model.vars[1])
    cov.mat=NULL
    coef.names=colnames(iter.betas)
  }
  #print("after get iter.betas")


  if(return.time==TRUE){
  iter.proxy.fit.stop=proc.time()
  time.values=c(time.values,
                list("iter.proxy.fit.time"=(iter.proxy.fit.stop-iter.proxy.fit.start)[[3]]))
  }

  ppart.list=list(NA)
  includepvalue=rep(FALSE,num.coefs)
  includeCI=rep(FALSE,num.coefs)
  if(return.treatment.pvalue==TRUE|return.treatment.CI==TRUE){
    pattern.treat=ifelse(length(treatment.var)==1,treatment.var,
                         paste0(treatment.var,collpase="|"))
    num.treat.coef=sum(grepl(pattern.treat,
                             names(conf.model$coefficients)))
    if(return.treatment.pvalue==TRUE){
      ppart.list=c(list(NA),
                  check_length_list(treat.ppart.list,"pvalue.partition",num.treat.coef),
                  rep(list(NA),num.coefs-num.treat.coef-1))
      includepvalue[seq(2,num.treat.coef+1)]=TRUE
    }
    if(return.treatment.CI==TRUE){
      includeCI[seq(2,num.treat.coef+1)]=TRUE
    }
  }



  check.lists=check_lists_coeftable(eps.list=epsilon.list, del.list=delta.list,
                                    alph.list=alphas.list,
                                    bd.sd.list=bd.sd.list,
                                    bd.mean.list=bd.mean.list,
                                    ppart.list=ppart.list,
                                    num.coefs=length(coef.names))


  epsilon.list=check.lists[[1]]
  delta.list=check.lists[[2]]
  alphas.list=check.lists[[3]]
  bd.sd.list=check.lists[[4]]
  bd.mean.list=check.lists[[5]]

  num.coefs=ncol(iter.betas)
  if(is.null(cov.mat)==TRUE){ #if not using sanitized residuals
  #get estimate and standard error
    san.summary=matrix(nrow=num.coefs,ncol=9)
    warning("before dp_coef_stats in is.null cov.mat==TRUE")
    for(idx in seq(1,length(coef.names))){
      temp.stats=dp_coef_stats(iter.betas[,idx],
                               bounds.sd=bd.sd.list[[idx]],
                               bound.mean=bd.mean.list[[idx]],
                               alphas=alphas.list[[idx]],
                               pvalue.partitions=ppart.list[[idx]],
                               epsilon = epsilon.list[[idx]],
                               delta=delta.list[[idx]],
                               include.san.pvalue=includepvalue[idx],
                               include.ci=includeCI[idx])
      temp.stats=c(temp.stats[1:2],
                   "san.CI.lower"=ifelse(includeCI[idx]==TRUE,temp.stats[3],NA),
                   "san.CI.upper"=ifelse(includeCI[idx]==TRUE,temp.stats[4],NA),
                   "confidence.level"=ifelse(includeCI[idx]==TRUE,temp.stats[5],NA),
                   "san.z.value"=ifelse(includepvalue[idx]==TRUE,temp.stats[(3+(3*as.numeric(includeCI[idx])))],NA),
                   "san.pvalue"=ifelse(includepvalue[idx]==TRUE,temp.stats[(4+(3*as.numeric(includeCI[idx])))],NA),
                   "epsilon"=temp.stats[(3+(3*as.numeric(includeCI[idx]))+(2*as.numeric(includepvalue[idx])))],
                   "delta"=temp.stats[(4+(3*as.numeric(includeCI[idx]))+(2*as.numeric(includepvalue[idx])))])
      san.summary[idx,]=temp.stats
    }
    cov.mat=stats::vcov(conf.model)
  }else{ #using sanitized residual variance

    warning(paste("before dp_coef_stats in is.null cov.mat==FALSE. class of bd.mean.list=",class(bd.mean.list),"length is",length(bd.mean.list),
                  "num.coefs is",length(coef.names)))
    san.summary=matrix(nrow=num.coefs,ncol=9)
    for(idx in seq(1,length(coef.names))){
      temp.stats=dp_coef_stats(iter.betas[,idx],
                               bounds.sd=bd.sd.list[[idx]],
                               bound.mean=bd.mean.list[[idx]],
                               alphas=alphas.list[[idx]],
                               pvalue.partitions=ppart.list[[idx]],
                               epsilon = epsilon.list[[idx]],
                               delta=delta.list[[idx]],
                               include.san.pvalue=includepvalue[idx],
                               include.ci=includeCI[idx])
      temp.stats=c(temp.stats[1:2],
                   "san.CI.lower"=ifelse(includeCI[idx]==TRUE,temp.stats[3],NA),
                   "san.CI.upper"=ifelse(includeCI[idx]==TRUE,temp.stats[4],NA),
                   "confidence.level"=ifelse(includeCI[idx]==TRUE,temp.stats[5],NA),
                   "san.z.value"=ifelse(includepvalue[idx]==TRUE,temp.stats[(3+(3*as.numeric(includeCI[idx])))],NA),
                   "san.pvalue"=ifelse(includepvalue[idx]==TRUE,temp.stats[(4+(3*as.numeric(includeCI[idx])))],NA),
                   "epsilon"=temp.stats[(3+(3*as.numeric(includeCI[idx]))+(2*as.numeric(includepvalue[idx])))],
                   "delta"=temp.stats[(4+(3*as.numeric(includeCI[idx]))+(2*as.numeric(includepvalue[idx])))])
      san.summary[idx,]=temp.stats
    }
  }
  rownames(san.summary)=colnames(iter.betas)
  colnames(san.summary)=c("san.coef","san.sd","san.CI.lower","san.CI.upper",
                          "confidence.level","san.z.value","san.pvalue","epsilon","delta")
  if(return.time==TRUE){
    san.summ.stop=proc.time()
    time.values=c(time.values,
                  list("san.summary.time"=(san.summ.stop-iter.proxy.fit.stop)[[3]]))
  }

  san.coefs=san.summary[,1]
  has.coef.name=names(san.coefs)[!is.na(san.coefs)] #only columns where coefficient!=NA
  san.coef.std.err=san.summary[,2]

  warning("before model.matrix")

  #model predictor matrix
  modMat=stats::model.matrix(predictor.formula, pred.newdata)
  modMat.cnames=colnames(modMat)
  not.in.modMat=has.coef.name[!(has.coef.name %in% modMat.cnames)]
  if(length(not.in.modMat)>0){ #if some variables don't appear
    zeros.mat=matrix(0,ncol=length(not.in.modMat),nrow=nrow(modMat))
    modMat=cbind(modMat,zeros.mat)
    colnames(modMat)=c(modMat.cnames,not.in.modMat)
  }
  modMat=modMat[,has.coef.name]

  ############################## PROBLEM HERE B/C MISSING A LEVEL OF BLOCK.
  #coef.rv=mvtnorm::rmvnorm(1, san.coefs[names(san.coefs)%in%has.coef.name],
  #                         cov.mat[rownames(cov.mat)%in%has.coef.name,colnames(cov.mat)%in%has.coef.name])

  sim.res=mvtnorm::rmvnorm(1, rep(0,nrow(modMat)), (san.mse*diag(nrow=nrow(modMat))))
  sim.response=as.numeric(c(modMat%*%matrix(san.coefs[names(san.coefs)%in% has.coef.name],ncol=1))+c(sim.res))

  if(return.time==TRUE){
    san.y.stop=proc.time()
    time.values=c(time.values,
                  list("san.response.time"=(san.y.stop-san.summ.stop)[[3]]))
  }

  #add synthetic response and reorder columns to match confidential data
  curr.synth.colnames=colnames(synth.data)
  conf.colnames=colnames(confidential.data)
  synth.data$response=sim.response
  colnames(synth.data)=c(curr.synth.colnames,response.var)
  synth.data=synth.data[,conf.colnames[conf.colnames %in% colnames(synth.data)]]
  output=list("synth.data"=synth.data)

  if(return.confidential.table==TRUE){
    output=c(output,list("confidential.summary"=summary(conf.model)))
  }
  if(return.san.summary==TRUE){
    san.summary=san.summary[,apply(san.summary,2,function(col)!all(is.na(col)))]
    output=c(output,list("san.summary"=san.summary))
  }
  if(return.time==TRUE){
    stop.time=proc.time()
    output=c(output,
             list("comp.time"=c(list("total.dp_synthdata.time"=(stop.time-start.time)[[3]]),
                                time.values)))
  }
  return(output)
}
