
#' Approximate DP Model Summary Table
#' @description Based on Algorithm 3 from paper.
#' @inheritParams hybrid_synth
#' @inheritParams synthdata_perturb_mvhist
#' @param epsilon.list a list where each entry is the epsilon value(s) for the
#'    corresponding model coefficient
#' @param delta.list a list where each entry is the delta value(s) for the
#'    corresponding model coefficient
#' @param bd.sd.list a list where each entry is the bounds.sd tuple for the
#'    corresponding model coefficient (see \code{\link{dp_confidence_interval}})
#' @param bd.mean.list a list where each entry is the bound.mean value for the
#'    corresponding model coefficient (see \code{\link{dp_confidence_interval}})
#' @param alphas.list a list where each entry is the alpha value(s) for the
#'    corresponding model coefficient (see the first 1 (or 3) values are passed to
#'    \code{\link{dp_confidence_interval}} and the last value is passed to
#'    \code{\link{dp_pvalue}})
#' @param p.partitions.list a list where each entry is the pvalue.partitions
#'    value for the corresponding model coefficient. (see \code{\link{dp_pvalue}})
#' @param synth.epsilon is the privacy parameter for making synthetic treatment and
#'    covariate data. If no \code{synth.data} is supplied and \code{use.san.residerr==TRUE},
#'    then synth.epsilon must be a positive value.
#' @param synth.delta is the privacy parameter for making synthetic treatment and
#'    covariate data. Default is 0.
#' @param num.iters the number of proxy response variable to produce
#' @param just.treatment indicator if just the treatment variable sanitized
#'    statistics (rather than all predictors) should be returned
#' @param return.confidential.table indicator if the confidential model summary
#'    table should be returned
#' @param return.time indicator if computation time should be returned
#' @param use.san.residerror indicator if a sanitized residual error should be used.
#' @returns a table with sanitized point estimate, standard deviation estimate,
#'    lower and upper 95\% confidence interval bounds, z value and P(>|z|)
#'    for each model coefficient.
#' @details
#' If only one value is given for \code{epsilon.list} or \code{delta.list}, then
#' this is taken to be the budget for the whole table and is divided evenly across the
#' For the elements of \code{returntypes}:
#' If "synth.data" is included, the full synthetic dataset is returned.
#' If "san.treat.effect", the sanitized (or privacy-preserving) treatment effect from fitting a model with the synthetic data is returned.
#' If "san.model" is included, the \code{\link{glm}} object of the model fitted with the synthetic data is outputted.
#' If "confidential.model" is included, the \code{glm} object of the model fitted with the confidential data is outputted without privacy protections.
#' If "comp.time" is included, the computation time of the algorithm is outputted. If "privacy.cost" is included, the overall privacy cost is returned as a vector of c(composed epsilon,composed delta).
#' Additional parameters ... are supplied to the \code{glm} function.
#' To specify random treatment assignment further. Use \code{\link{synthdata_perturb_mvhist}} with \code{with.treatment=TRUE} to create your synthetic data first.
#' Users may bring their own synthetic data and use \code{\link{treatment_assign}} to further control the specifications of the random treatment assignment methodology.
#' Then the synthetic dataset can be supplied to this function as \code{synth.data} input.
#'
#' @importFrom stats glm as.formula coefficients
#' @export

dp_model_summary=function(formula,
                           confidential.data,
                           epsilon.list,
                           delta.list=0,
                           bd.sd.list=c(2^(-15),2^(15)),
                           bd.mean.list,
                           alphas.list=0.05,
                          p.partitions.list,
                           num.iters=1e4,
                           treatment.var=NULL,
                           rseed=NA,
                           return.time=FALSE,
                           just.treatment=FALSE,
                           return.confidential.table=FALSE,
                          use.san.residerror=TRUE,
                          synth.data=NULL,
                          synth.epsilon=NA,
                          synth.delta=0,
                          continuous.vars=NULL,
                          num.bin=NULL,
                          bin.param=NULL,
                          add.cont.variation=TRUE,
                          assign.type="simple",
                          blocks=NULL,
                          clusters=NULL,
                          within.blocks=TRUE,
                           ...){

  start.time=proc.time()
  time.values=NULL
  #check inputs
  stopifnot(is.data.frame(confidential.data))

  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }


  psf=parse_formula(formula=formula,
                    confidential.data = confidential.data,
                    treatment.var = treatment.var)
  formula=psf[[1]]
  model.vars=psf[[2]]
  response.var=psf[[3]]
  treatment.var=psf[[5]]

  conf.model.start=proc.time()
  #### Step 0: In paper
  #fit the confidential data to the model and get estimated coefficient and standard error
  conf.model<-stats::glm(formula=formula,
                         family="gaussian",
                         data=confidential.data,...)
  ####
  if(return.confidential.table==TRUE){
  #get confidential coefficients and summary table
  conf.summary=base::summary(conf.model)$coefficients
  colnames(conf.summary)=paste0("Confidential ",colnames(conf.summary))
  zhalfalpha=qnorm(1-(sapply(alphas.list,function(x)sum(x,na.rm = TRUE))/2))
  margin.of.error=conf.summary[,2]*zhalfalpha
  conf.summary$`Confidential CI Lower`=conf.summary[,1]-margin.of.error
  conf.summary$`Confidential CI Upper`=conf.summary[,1]+margin.of.error
  conf.summary=t(conf.summary[,c(1,2,5,6,3,4)])

  }
  conf.model.stop=proc.time()
  time.values=c(time.values,
                list("confidential.model.time"=(conf.model.stop-conf.model.start)[[3]]))
      if(just.treatment==TRUE){ #just treatment coefficients sanitized
        treat.coefs=unlist(
          sapply(treatment.var,
                 function(tr.var)paste0(tr.var,levels(confidential.data[,tr.var])[-1])))
        num.coefs=length(treat.coefs)
      }else{
        num.coefs=sum(!is.na(conf.summary[,1]))
      }


  if(use.san.residerror==TRUE){
    # if synthetic data is supplied, check it is data.frame w/ all the model.vars
    if((is.null(synth.data)==FALSE)&&
       (is.data.frame(synth.data)==TRUE)&&
       (sum(c(model.vars[model.vars!=response.var])%in%colnames(synth.data))==length(model.vars)-1)){
      synth.data<-synth.data
    }else{ #if synth.data not supplied, check it closely
      synth.data=handle_synthetic_data(confidential.data=confidential.data,
                                        synth.data=synth.data,
                                        treatment.var=treatment.var,
                                        synth.vars=model.vars[model.vars!=response.var],
                                        model.vars=model.vars,
                                        epsilon=synth.epsilon,
                                        delta=synth.delta,
                                        continuous.vars=continuous.vars,
                                        num.bin=num.bin,
                                        bin.param=bin.param,
                                        add.cont.variation=add.cont.variation,
                                        assign.type=assign.type,
                                        blocks=blocks,
                                        clusters=clusters,
                                        within.blocks=within.blocks,
                                        return.time=TRUE)
      if(return.time==TRUE){
      time.values=c(time.values,
                    list("synthetic.T.X.time"=attr(synth.data,"comp.time")))
      }
    }

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

    list2env(dp_iter_hybrid(conf.model=conf.model,
                           confidential.data=confidential.data,
                           synth.data=synth.data,
                           num.iters=num.iters,
                           model.vars=model.vars,
                           response.var=response.var,
                           mse.epsilon=mse.epsilon,
                           mse.delta=mse.delta,
                           mse.bd.sd=mse.bd.sd),envir=globalenv())


    coef.names=names(conf.model$coefficients)
    coef.names=coef.names[!is.na(conf.model$coefficients)]
    ############ HERERERERERERERERER #####################
  }else{ #use.san.residerr==FALSE
  #generate num.iters proxy response variables from confidential model
  #fit a new model and return the coefficients
  iter.betas=iter_hybrid(conf.model=conf.model,
                         confidential.data = confidential.data,
                         num.iters = num.iters,
                         predictor.vars = model.vars[-1],
                         response.var=response.var)
  cov.mat=matrix(NA,nrow=num.coefs,ncol=num.coefs)
  coef.names=names(conf.model$coefficients)
  coef.names=coef.names[!is.na(conf.model$coefficients)]
  colnames(cov.mat)=coef.names
  rownames(cov.mat)=coef.names
  }
  if(return.time==TRUE){
    iter.proxy.fit.stop=proc.time()
    time.values=c(time.values,
                  list("iter.proxy.fit.time"=(iter.proxy.fit.stop-conf.model.stop)[[3]]))
  }
  #warning("in dp_model_summary before check lists")
  san.summ.start=proc.time()
  #check parameter lists are correct length (correct them if not)
      check.lists=check_lists_coeftable(epsilon.list,
                                        delta.list,
                                        alphas.list,
                                        bd.sd.list,
                                        bd.mean.list,
                                        p.partitions.list,
                                        num.coefs=length(coef.names))
      if(just.treatment==TRUE){
        diff.coef.num=ncol(conf.summary)-num.coefs-1
        check.lists=lapply(seq(1,5),
                           function(idx)c(list(NA), #intercept term
                                          check.lists[[idx]], #treatment terms
                                          rep(list(NA),diff.coef.num))) #NA for others

      }
      epsilon.list=check.lists[[1]]
      delta.list=check.lists[[2]]
      alphas.list=check.lists[[3]]
      bd.sd.list=check.lists[[4]]
      bd.mean.list=check.lists[[5]]
      p.partitions.list=check.lists[[6]]


      #warning(paste("class of bd.mean.list",class(bd.mean.list),"length is",length(bd.mean.list)))
      san.summary=sapply(seq(1,ncol(conf.summary)),
                 function(idx)
                   dp_coef_stats(iter.betas[,idx],
                                 san.sd=cov.mat[idx,idx],
                                 bounds.sd=bd.sd.list[[idx]],
                                 bound.mean=bd.mean.list[[idx]],
                                 alphas=alphas.list[[idx]],
                                 epsilon = epsilon.list[[idx]],
                                 delta=delta.list[[idx]],
                                 pvalue.partitions=p.partitions.list[[idx]],
                                 ))
      if(return.time==TRUE){
        san.summ.stop=proc.time()
        time.values=c(time.values,
                      list("san.summary.time"=(san.summ.stop-san.summ.start)[[3]]))
      }
      row.names(san.summary)=c("Sanitized Estimate","Sanitized Std. Error",
                       "Sanitized CI Lower","Sanitized CI Upper",
                       "Confidence Level","Sanitized z value","Sanitized P(>|z|)",
                       "Coefficient Epsilon","Coefficient Delta")
      if(return.confidential.table==TRUE){
        if(just.treatment==FALSE){ #just output treatment statistics
          output=rbind(san.summary,conf.summary)
          attr(output,"priv.cost")=c("epsilon"=sum(output["Coefficient Epsilon",],na.rm=T),
                                     "delta"=sum(output["Coefficient Delta",],na.rm=T))
        }else{
          attr(san.summary,"priv.cost")=c("epsilon"=sum(san.summary["Coefficient Epsilon",],na.rm=T),
                                     "delta"=sum(san.summary["Coefficient Delta",],na.rm=T))

          output=list("san.summary"=san.summary,"confidential.summary"=conf.summary)
        }
      }else{ #if return.confidential.table==FALSE
        output=san.summary
        attr(output,"priv.cost")=c("epsilon"=sum(output["Coefficient Epsilon",],na.rm=T),
                                   "delta"=sum(output["Coefficient Delta",],na.rm=T))

      }
      if(return.time==TRUE){
        stop.time=proc.time()
        return(list("summary.tables"=output,
                    "comp.time"=c(list("total.dp_model_summary.time"=(stop.time-start.time)[[3]]),
                                  time.values)))
      }else{
        return(output)
      }
    }




