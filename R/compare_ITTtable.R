#' Dataframe to Compare ITT Across Synthetic and Confidential Data
#'
#' @inheritParams ITTtable
#' @param data.list a named list of data sets to compare the ITT tables of. If no names are
#'      given data sets will be named Data 1, Data 2, ...
#' @param model.names Names for each regression model as they should appear in the table.
#'      If NULL (default), then the variable names for the response of each
#'      \code{reg.models} wil be used instead.
#' @param include.ci.overlap indicator if the overlap of the ITT confidence intervals
#'      should be included in the outputted dataset.
#' @param ci.overlap.ref.name is the name of the dataset in the \code{data.list} that
#'      should be used as the reference for all the other datasets to be compared to.
#'      This must not be \code{NULL} if \code{include.ci.overlap==TRUE}.
#' @param ci.confidence is the confidence level (as a decimal) of the intervals. Default is 0.95.
#' @param ... Further arguments will be passed to the \link{kable} function.
#' @returns a data.frame with columns: name of dataset the model is fit on; ITT;
#'      Std. Err; and p-value (or Adj. p-value); number of observations; treatment name;
#'      response variable/regression model; ITT Confidence Interval Lower and Upper Bounds;
#'      and Confidence interval overlap with a reference dataset. There is a row for
#'      each treatment level and control for each response variable/regression model
#'      and dataset.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom stats qnorm
#' @importFrom dplyr bind_rows
#'
#'
#'
#' @export
compare_ITTdata=function(data.list,
                                reg.models,
                                response.vars=NULL,
                                families="gaussian",
                                treat.vars,
                                control.var="control",
                                bonferroni.npvals=NULL,
                                mult.test.correct=c("treatments","responses"),
                                model.names=NULL,
                                stderr.func=NULL,
                         only.control.mean=FALSE,
                                include.ci.overlap=FALSE,
                         ci.overlap.ref.name=NULL,
                         ci.confidence=0.95,
                         pivot.xtreat=FALSE,
                         col.var.names=NULL,
                         treat.names=NULL,
                                ...){

  #Make sure data.list is a list of named objects
  if(is.list(data.list)==FALSE){
    data.list=list("Data"=data.list)
  }
  if(is.null(names(data.list))==TRUE){ #if no names provided name them Data 1, Data 3,...
    names(data.list)=paste("Data",seq(1,length(data.list)))
  }

  if(include.ci.overlap==TRUE){
    add.stars=FALSE
    include.df=TRUE
  }else{
    add.stars=TRUE
    include.df=FALSE
  }

  #make the ITT table for each data set in the data.list and stack them
  table.combine=dplyr::bind_rows(
    lapply(data.list,
           function(df)
               ITTtable(data=df,families=families,
                      reg.models=reg.models,
                      response.vars=response.vars,
                      treat.vars=treat.vars,
                      control.var=control.var,
                      bonferroni.npvals=bonferroni.npvals,
                      mult.test.correct = mult.test.correct,
                      pivot.xtreat=FALSE,
                      add.pval.stars=add.stars,
                      model.names=model.names,
                      stderr.func=stderr.func,
                      only.control.mean=TRUE,
                      include.df=include.df)),.id="data")
  if(setequal(unique(table.combine$data),names(data.list))==FALSE){
    table.combine$data=as.factor(table.combine$data)
    levels(table.combine$data)=names(data.list)
  }
  if(include.ci.overlap==TRUE){ #if including confidence interval overlap...
    stopifnot((is.null(ci.overlap.ref.name)==FALSE))#&&(ci.overlap.ref.name%in% unique(table.combine$data)))
    stopifnot(ci.confidence>0 & ci.confidence<1)
    #z multiplier
    z.half.alpha=stats::qnorm(0.5*(1+ci.confidence))
    #get upper and lower bounds for each confidence interval.
    cis.all=sapply(seq(1,nrow(table.combine)),
                   function(i){
                     table.combine$ITT[i]+
                       (c(-z.half.alpha,z.half.alpha)*table.combine$StdErr[i])})
    table.combine$ci.lower=cis.all[1,]
    table.combine$ci.upper=cis.all[2,]
    #indicator of reference data set rows
    ref.indic=table.combine$data==ci.overlap.ref.name
    ref.tab=table.combine[ref.indic,]
    table.combine$ci.overlap=1 #initialize confidence interval overlap.
    not.ref=unique(as.character(table.combine$data))
    not.ref=not.ref[not.ref!=ci.overlap.ref.name]
    warning(paste("In compare_ITT. not.ref is",paste0(not.ref,collapse=", ")," and table.combine$data is",paste(unique(table.combine$data),collapes=", ")))
    for(dta in not.ref){ #for data which is not the reference data...
      dta.indic=table.combine$data==dta #indicator of the dataset
      table.combine$ci.overlap[dta.indic]= #get ci overlap
        CI_overlap(table.combine$ci.lower[dta.indic],table.combine$ci.upper[dta.indic],
                   ref.tab$ci.lower,ref.tab$ci.upper)
    }
  }

  if(pivot.xtreat==TRUE){
    #### deal with mult.test.correct and bonferroni.npvals ###
    pvalname="Pvalue"
    if((is.null(bonferroni.npvals)==FALSE)&&(bonferroni.npvals>1)){
        pvalname="AdjPvalue"
      }else if((is.null(bonferroni.npvals)==TRUE)&(is.null(mult.test.correct)==FALSE)){
        pvalname="AdjPvalue"
      }
    treat.vars=unique(table.combine$Treatment)
    pivot.values=c("ITT","StdErr",pvalname)
    if(include.ci.overlap==TRUE){
      pivot.values=c(pivot.values,"ci.lower","ci.upper","ci.overlap")
    }
    table.combine=tidyr::pivot_wider(table.combine,names_from = "Treatment",
                                 values_from = pivot.values)

    #if only.control.mean==TRUE then remove StdErr_control and AdjPvalue_control columns
    if((is.null(control.var)==FALSE)&(only.control.mean==TRUE)){
      control.cols=grepl("Control",colnames(table.combine))
      control.mean=table.combine[,"ITT_Control"]
      table.combine=table.combine[,!control.cols]
      other.colnames=colnames(table.combine)
      table.combine$ControlMean=unlist(control.mean)
      treat.vars=treat.vars[treat.vars!=control.var]
    }

    ### reorder columns to be keep columns for treatment values together
    colnms.estimates=colnames(table.combine)[grepl("^ITT",colnames(table.combine))]
    tr.var.suffix=sapply(base::strsplit(colnms.estimates,"_"),function(x)x[2])
    temp.colnames=base::expand.grid(pivot.values,tr.var.suffix)
    temp.colnames=paste(temp.colnames[,1],temp.colnames[,2],sep="_")
    other.colnames=colnames(table.combine)[!(colnames(table.combine)%in%temp.colnames)]
    table.combine=table.combine[,c(other.colnames,temp.colnames)]

    #get the desired columns in the desired order
    cnames=colnames(table.combine)
    col.order=c(ifelse("ModelName"%in%cnames,"ModelName","Response"),"data")
    if(only.control.mean==TRUE){col.order=c(col.order,"ControlMean")}
    for(tr.var in treat.vars){
      col.order=c(col.order,paste(pivot.values,tr.var,sep="_"))
    }
    col.order=c(col.order,"nObs")
    table.combine=table.combine[,col.order]
  }
  return(table.combine)
}

