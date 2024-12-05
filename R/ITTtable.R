#' Get ITT table values for a given regression model.
#' @param data is the data to use to fit the model
#' @param reg.model is the regression formula to fit
#' @param family is the family for the glm to use to fit the model. Default is "gaussian".
#' @param treat.vars a vector of variable names, indices, or logicals for the
#'    columns of the data that are treatment variables
#' @param bonferroni.npvals is a numeric value for the number of p-values to use
#'    a Bonferroni correction to adjust for multiple testing. The default is
#'    the number of treatment variables and level coefficients in the model
#'@param control.var the string for the name of the control variable in the data.
#'    If \code{control.var==NULL} (default), then the control mean will not be included in the table
#'@param add.pval.stars if \code{TRUE}, then standard errors will be formatted with
#'    3 decimal places in square brackets followed by one "*" if unadjusted p-value is <0.1,
#'    two if <0.05, and three if <0.01. For example, \code{"[0.110]**"}. If
#'    \code{add.pval.stars==FALSE}, then standard error will remain as numeric elements.
#'@param stderr.func (optional) is a function that takes the \link{glm} model output
#'    and returns the standard error of the estimated treatment coefficients. If NULL
#'    then the standard errors from the summary table will be used.
#'@param ... additional inputs will be passed to \code{\link{glm}}
#'@return data.frame with columns for the estimates, standard error (with stars
#'    for unadjusted p-value if \code{add.pval.stars==TRUE}),  Bonferroni corrected p-values,
#'    treatment variable/level names (or control), and response variable name.
#'    With rows for each treatment variable/level and level and the control variable
#'    (if \code{control.var!=NULL}).
#'
#'@importFrom stats glm as.formula var
#'
#'@export
ITTtable_oneresponse=function(data,reg.model,family="gaussian",
                              treat.vars,bonferroni.npvals=NULL,control.var=NULL,
                              add.pval.stars=TRUE,stderr.func=NULL,incl.df=FALSE,...){

  #fit model
  mod.fit=stats::glm(reg.model,family=family,data=data,...)
  mod.summary=summary(mod.fit)$coefficients
  mod.coefs=rownames(mod.summary)

  if(is.character(treat.vars)==FALSE){ #force treat.vars to be variable names
    treat.vars=colnames(data)[treat.vars]
  }
  #get factor levels if treat.vars is a factor
  treat.vars.factors=sapply(data[,c(treat.vars)],
                            function(x)base::is.factor(x)|base::is.character(x))
  if(sum(treat.vars.factors)==0){ #if no treat.vars are factors
    # if(length(base::setdiff(treat.vars,mod.coefs))>0){ #check all treat.vars in formula
    #   stop("Specified treatment variables are not in formula.")
    # }
    get.coefs=treat.vars
  }else{ #if there are/is treat.vars that are/is factor
    get.coefs=treat.vars[!treat.vars.factors] #not factor treat.vars
    if(length(base::setdiff(get.coefs,mod.coefs))>0){ #check not factor treat.vars in formula
      stop("Specified treatment variables are not in formula.")
    }
    remain.treat.vars=treat.vars[treat.vars.factors] #treat.vars that are factors

    treat.vars.levels=unlist( #list of levels for each factor treat.vars
      lapply(treat.vars[treat.vars.factors],
             function(x)paste0(x,t(base::unique(data[,x])))))
    #get coefficient rowname for all treat.vars
    get.coefs=c(get.coefs,treat.vars.levels[treat.vars.levels%in%mod.coefs])
  }

  sub.summary=mod.summary[mod.coefs%in%c(get.coefs),,drop=F] #regression summary info for treat.vars
  #if no bonferroni.npvals supplied, use number of coefficients for treat.vars
  if(is.null(bonferroni.npvals)==TRUE){
    bonferroni.npvals=nrow(sub.summary)
  }

  response.var=base::as.character(stats::as.formula(reg.model))[2]
  estimates=base::unname(sub.summary[,1])

  if(is.null(stderr.func)==TRUE){
    se.num=sub.summary[,2]
    pval=base::unname(sub.summary[,4])
  }else{
    se.num=stderr.func(mod.fit)
    se.num=se.num[names(se.num)%in%c(get.coefs)]
    z=abs(estimates)/se.num
    pval=pnorm(z,lower.tail=FALSE)
  }

  ### extracting ITT table information
  #get stars for the pvalues of each coefficient, format se with stars
  if(add.pval.stars==TRUE){
    se.stars=ifelse(pval<0.1,"*","")
    se.stars[pval<0.05]="**"
    se.stars[pval<0.01]="***"
    se.str=paste0("[",as.character(base::format(round(se.num,3),nsmall=3)),"]",se.stars)
  }else{
    se.str=se.num
  }


  adj.p=base::pmin(base::unname(pval)*bonferroni.npvals,1)


  warning(paste("get estimates",paste0(estimates,collapse=", "),
                "get stderr",paste0(se.str,collapse=", "),
                "get pvals",paste0(adj.p,collapse=", "),
                "get treatment",paste0(rownames(sub.summary),collapse=", ")))

  tab.out=data.frame("ITT"=estimates,
                     "StdErr"=se.str,
                     "AdjPvalue"=adj.p,
                     "nObs"=length(mod.fit$residuals),
                     "Treatment"=rownames(sub.summary),
                     "Response"=rep(response.var,length(adj.p)))
  if(incl.df==TRUE){
    tab.out$deg.freedom=mod.fit$df.residual
  }



  if(is.null(control.var)==FALSE){ #if we are adding mean of treatment control
    #control.var should take value 0,1 or FALSE,TRUE.
    indic.control=(data[,control.var]==1|data[,control.var]==TRUE)
    control.responses=unlist(data[indic.control,response.var])
    se.control=base::sqrt(stats::var(control.responses,na.rm=T)/sum(!is.na(control.responses)))
    se.control=ifelse(add.pval.stars==TRUE,
                      paste0("[",as.character(base::format(round(se.control,3),nsmall=3)),"]"),
                      se.control)
    #add row for control with "Estimate" being the mean response for control treatment,
    # "StdErr" is standard error of mean, no AdjPvalue included
    controldf=data.frame("ITT"=mean(control.responses,na.rm=T),
                         "StdErr"=se.control,
                         "AdjPvalue"=NA,
                         "nObs"=length(mod.fit$residuals),
                         "Treatment"="Control",
                         "Response"=response.var)
    tab.out=rbind(tab.out,controldf)
  }

  return(tab.out)
}

#' ITT table for multiple models
#'
#' @inheritParams ITTtable_oneresponse
#' @param reg.models vector of strings or formulas for the regression models. If
#'    \code{reg.models==NULL} (default), then \code{treat.vars} and \code{covariate.vars},
#'    cannot be \code{NULL} and are used to make a linear regression model for each
#'    \code{response.vars}. If \code{reg.models} is only one formula or string, then
#'    the predictors of that formula are used with each \code{response.vars} to make
#'    a vector of models.
#'@param response.vars a vector of names, indicators, or indices of the response
#'    variables within \code{data}.
#'@param families a vector of the family for each \code{reg.models} \link{\code{glm}}. The
#'    \code{length(families)} must be 1 or equal to \code{length(reg.models)} (or \code{response.vars}).
#'    If the length is 1, then the same distribution family is used for all models. Default is
#'    "gaussian".
#'@param covariate.vars a vector of names, indicators, or indices of the covariate
#'    variables within \code{data} to create the regression model formulas if \code{reg.models==NULL}.
#'@param mult.test.correct a vector with elements "treatments" or "responses". If
#'    \code{bonferroni.npvals==NULL}, then \code{mult.test.correct} is used to determine
#'    the bonferroni correction. If \code{mult.test.correct=="treatments"}, then
#'    \code{bonferroni.npvals} is the number of treatment variables/levels. If
#'    \code{mult.test.correct=="responses"}, then \code{bonferroni.npvals} is the
#'    length of \code{reg.models} or \code{response.vars}. If \code{mult.test.correct}
#'    is the vector of "treatments" and "responses", then \code{bonferroni.npvals} is
#'    the number of treatment variables/levels plus the number of treatment variables/levels.
#'@param pivot.xtreat if \code{TRUE}, then the resulting data.frame is pivoted wider
#'    (see \link{\code{pivot_wider}}) to have "Estimate", "StdErr", "AdjPvalue" columns for
#'    each treatment variable/level (control.var if applicable). If \code{control.var!=NULL}
#'    and \code{only.control.mean==TRUE}, then the control mean (not the StdErr or AdjPvalue)
#'    is also a column. If \code{pivot.xtreat==FALSE},then data frame is left in long form
#'    with a "Treatment" column.
#'@param only.control.mean is indicator if the wider pivot data.frame only includes
#'    the StdErr and AdjPvalue columns.
#'@param model.names (optional) is a vector of names for each regression model.
#'@param stderr.func (optional) is a function that takes the \link{glm} model output
#'    and returns the standard error of the estimated treatment coefficients. If NULL
#'    then the standard errors from the summary table will be used.
#'@return a data.frame with the coefficient estimates, standard errors, bonferroni
#'  adjusted p-values for each treatment variable/level and each regression
#'  model (or response variables).
#'
#'@importFrom dplyr bind_rows
#'@importFrom tidyr pivot_wider
#'@importFrom stats as.formula
#'
#'@export
ITTtable=function(data,reg.models=NULL,
                  response.vars=NULL,
                  families="gaussian",
                  treat.vars,
                  covariate.vars=NULL,
                  control.var=NULL,
                  bonferroni.npvals=NULL,
                  mult.test.correct=c("treatments","responses"),
                  pivot.xtreat=TRUE,
                  only.control.mean=TRUE,
                  add.pval.stars=TRUE,
                  model.names=NULL,
                  stderr.func=NULL,
                  include.df=FALSE,
                  ...){

  if(is.character(treat.vars)==FALSE){ #make treat.vars to a vector of names
    treat.vars=colnames(data)[treat.vars]
  }


  #### deal with reg.models and corresponding inputs ###
  # if reg.models!=NULL and has length>1 and != length(response.vars),
  #       then use reg.models as is
  # if reg.models!=NULL and has length==1 and length(response.vars)>1,
  #       then use predictors from reg.models with response.vars
  #if reg.models==NULL,
  #       then use response.vars, treat.vars, and covariate vars to make reg.models
  stopifnot(is.null(response.vars)==FALSE|is.null(reg.models)==FALSE)
  if(is.null(reg.models)==FALSE){ #reg.models supplied
    if(is.null(covariate.vars)==FALSE){ #if covariates!=NULL, warn that reg.model is used
      warning("reg.models and covariate.vars !=NULL. reg.models are used. covariate.vars are not.")
    }
    #if reg.models and response.vars supplied
    if(is.null(response.vars)==FALSE){
      if(length(reg.models)!=length(response.vars)){ #uneven length
      #if one reg.model given and multiple response.vars, use predictors from reg.model with response.vars
      if((length(response.vars)>1)&(length(reg.models)==1)){
        ##get just the predictors of the reg.models
      reg.models=as.character(reg.models)
      if(sum(grepl("~",reg.models))>0){ #if input reg.models is formula,it contains "~"
        reg.models=as.character(stats::as.formula(reg.models))[3] #reg.models predictors
      }
      #if input reg.models was a string without "~", then assume it is just the predictors
      #combine response.vars with regression predictor formula string
      reg.models=paste(response.vars,reg.models,sep="~")
    }else if(length(reg.models)>1){ #end length(reg.models)==1 & length(response.vars)>1
      warning(paste0("Length of reg.models is >1 but != length(response.vars).\n",
              "Default will only use reg.models and ignore response.vars"))
    }#end length(reg.models)>1
    }else{ #end uneven length => if length(reg.models)==length(response.vars)
      if(length(reg.models)==1){
        y.from.reg.models=as.character(stats::as.formula(reg.models))[2]
      }else{ #end length=1, more than one regression model
      y.from.reg.models=sapply(reg.models,function(x)as.character(stats::as.formula(x))[2])
      if(base::setequal(y.from.reg.models,response.vars)==FALSE){
        warning(paste0("Inputted response.vars do not match the responses variable of reg.models.\n",
                       "Default reg.models will be used and response.vars will be ignored."))
      }
      }#end length>1
    } #end even length
} #end response.vars!=NULL
    response.vars=sapply(reg.models,function(x)as.character(stats::as.formula(x))[2])
  }else{ #if no reg.models supplied
    #force response.vars, and covariate.vars to be (vectors of) strings
    if(is.character(response.vars)==FALSE){
      response.vars=colnames(data)[response.vars]
    }
    if(is.character(covariate.vars)==FALSE){
      covariate.vars=colnames(data)[covariate.vars]
    }
    #if covariate.vars is a vector of variables force it to be an additive formula
    if(length(covariate.vars)>1){
      covariate.vars=paste0(covariate.vars,collapse="+")
    }
    #make vector of reg.models
    reg.models=paste0(response.vars,
                      paste0(paste0(c(treat.vars),collapse="+"),"+",covariate.vars), #predictors
                      sep="~") #response.vars seperated from predictors with "~"
  }
  if(length(response.vars)==1){
    data[,response.vars]=as.numeric(as.character(data[,response.vars]))
  }else{
  data[,response.vars]=apply(data[,response.vars],2,function(x)(as.numeric(as.character(x))))
  }
  ####

  ### deal with families input ###
  if(length(families)!=length(reg.models)){ #if families and reg.models have unequal length
    if(length(families)==1){ #if length(families)==1, repeat it for all reg.models
      warning("Only one family value supplied. It is used for all models.")
      families=rep(families,length(reg.models))
    }else{ #Error for unequal families and reg.models length
      stop(paste0("Length of families, ",length(families),
                  ", not equal to length of reg.models, ",length(reg.models),"."))
    }
  }

  #### deal with mult.test.correct and bonferroni.npvals ###
  if(is.null(bonferroni.npvals)==TRUE){ #if ... then use mult.test.correct
    if(is.null(mult.test.correct)==TRUE){
      bonferroni.npvals=1
      pvalname="Pvalue"
    }else{
      pvalname="AdjPvalue"
    #check that mult.test.correct element(s) are either "response(s)" or "treatment(s)".
    if(sum(grepl("response|treatment",c(mult.test.correct)))==0){
      stop("When bonferroni.npvals==NULL, mult.test.correct must contain or be equal to 'treatments' or 'responses'.")
    }
    #if mult.test.correct contains "response(s)", correction for number of reg.models
    bonferroni.npvals=ifelse(sum(grepl("response",c(mult.test.correct)))>=1,length(reg.models),1)
    if(sum(grepl("treatment",c(mult.test.correct)))>=1){ #if mult.test.correct contains "treatment(s)"
      #get which treat.vars are factors
      treat.vars.factors=sapply(data[,c(treat.vars)],
                                function(x)base::is.factor(x)|base::is.character(x))
      n.treat.vars.levels=length(treat.vars[!treat.vars.factors]) #number of numeric treat.vars
      if(sum(treat.vars.factors)>1){ #if there are factor treat.vars

        #check if reg.models used fixed intercepts (not fitted)
        #if so one treatment factor will have coefficients fitted for all levels
        predictor.model=sapply(reg.models,
                               function(x)base::trimws(as.character(as.formula(x))[3]))
        has.fixed.intercept=grepl("^[0-9]",predictor.model)
        if(sum(has.fixed.intercept)>0){
          warning(paste0("Some reg.models have fixed intercepts.\n",
                         "Thus all levels of the first treatment variable will have a estimated coefficient."))
        }
        #get number of levels -1 for each factor treat.vars
        n.tr.factor.levels=sapply(data[,c(treat.vars[treat.vars.factors])],
                                  function(x)length(base::unique(x))-1)
        #number of treatment variable/level coefficients fitted
        n.treat.vars.levels=n.treat.vars.levels+sum(n.tr.factor.levels)+has.fixed.intercept
      } #end treat.vars include factor variables.
      bonferroni.npvals=bonferroni.npvals+n.treat.vars.levels
    } #end mult.test.correct contains "treatment(s)"
  }#end mult.test.correct!=NULL
}else{#end input bonferroni.pvals==NULL
  pvalname=ifelse(bonferroni.npvals==1,"Pvalue","AdjPvalue")
}

  ###

  #if model.names provided, then it must be a vector of length equal to reg.models or response.vars
  if(is.null(model.names)==FALSE){
    stopifnot(length(model.names)==length(reg.models))
  }

  if(length(reg.models)==1){
    table.out=ITTtable_oneresponse(data=data,
                                reg.model=reg.models,
                                family=families,
                                treat.vars=treat.vars,
                                bonferroni.npvals =bonferroni.npvals,
                                control.var=control.var,
                                add.pval.stars=add.pval.stars,
                                stderr.func=stderr.func,...)

    if(is.null(model.names)==FALSE){
      table.out$ModelName=model.names
      }
  }else{
  ### internal function which calls ITTtable_oneresponse, and adds model.names if
  ### model.names!=NULL. Inputs to ITTtable_oneresponse are determined by and index value.
  one_response_func=function(idx,
                             df=data,
                             reg.mods=reg.models,
                             fams=families,
                             tr.vars=treat.vars,
                             bonf.npv=bonferroni.npvals,
                             ctrl.var=control.var,
                             plus.stars=add.pval.stars,
                             mod.names=model.names,
                             stderr.f=stderr.func,
                             incl.df=include.df,
                             ...){
    y.var=as.character(stats::as.formula(reg.mods[idx]))[2]
    tab.df=ITTtable_oneresponse(data=df,
                                reg.model=reg.mods[idx],
                                family=fams[idx],
                                treat.vars=tr.vars,
                                bonferroni.npvals =bonf.npv,
                                control.var=ctrl.var,
                                add.pval.stars=plus.stars,
                                stderr.func=stderr.f,...)

    if(is.null(model.names)==FALSE){
      tab.df=cbind(tab.df,mod.names[idx])
      colnames(tab.df)=c(colnames(tab.df)[seq(1,ncol(tab.df)-1)],"ModelName")
    }
    return(tab.df)
  }
  ####

  #combine ITTtable_oneresponse for each reg.models into one data.frame
  table.out=dplyr::bind_rows(lapply(seq(1,length(reg.models)),
                                    function(idx)
                                      one_response_func(idx,...)))
  }
  cnames=colnames(table.out)
  cnames[grepl("Pvalue",cnames,ignore.case = FALSE)]=pvalname
  colnames(table.out)=cnames

  #if pivot.xtreat==TRUE, pivot data wider with cols for each treatment variable/level coefficient
  if(pivot.xtreat==TRUE){
    table.out=tidyr::pivot_wider(table.out,names_from = "Treatment",
                                 values_from = c("ITT","StdErr",pvalname))
    #if only.control.mean==TRUE then remove StdErr_control and AdjPvalue_control columns
    if((is.null(control.var)==FALSE)&(only.control.mean==TRUE)){
      control.cols=grepl("Control",colnames(table.out))
      control.mean=table.out[,"ITT_Control"]
      table.out=table.out[,!control.cols]
      other.colnames=colnames(table.out)
      table.out$ControlMean=unlist(control.mean)
    }

    ### reorder columns to be keep columns for treatment values together
    colnms.estimates=colnames(table.out)[grepl("^ITT",colnames(table.out))]
    tr.var.suffix=sapply(base::strsplit(colnms.estimates,"_"),function(x)x[2])
    temp.colnames=base::expand.grid(c("ITT_","StdErr_",paste0(pvalname,"_")),tr.var.suffix)
    temp.colnames=paste0(temp.colnames[,1],temp.colnames[,2])
    other.colnames=colnames(table.out)[!(colnames(table.out)%in%temp.colnames)]
    table.out=table.out[,c(other.colnames,temp.colnames)]
  }#end pivot.xtreat==TRUE
  return(table.out)
}


