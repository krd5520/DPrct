#' Print Table to Compare ITT Across Synthetic and Confidential Data
#'
#' @inheritParams ITTtable
#' @param data.list a named list of data sets to compare the ITT tables of. If no names are
#'      given data sets will be named Data 1, Data 2, ...
#' @param model.names Names for each regression model as they should appear in the table.
#'      If NULL (default), then the variable names for the response of each
#'      \code{reg.models} wil be used instead.
#' @param digits a value for how many digits to round the values of the table to
#' @param label (Optional) a label for the table
#' @param caption (Optional) a caption for the table
#' @param col.var.names (Optional) the variable names of the columns in the order they
#'      should appear in the table.
#' @param treat.names (Optional) the names of the treatments as they should appear on
#'      on the table. If NULL (default) then the variable names will be used.
#' @param ... Further arguments will be passed to the \link{kable} function.
#' @returns a kableExtra table with columns: Outcome; Control mean (if \code{control.var!=NULL});
#'      ITT, Std. Err, and p-value (or Adj. p-value) for each treatment; and
#'      number of observations. The columns are grouped by the treatment level.
#'      The rows are grouped by data used for the regression with a row for each
#'      response variable/regression model for each data set.
#'
#' @importFrom knitr kable
#' @importFrom kableExtra add_header_above pack_rows
#'
#' @export
print_compare_ITTtable=function(data.list,
                                reg.models,
                                response.vars=NULL,
                                families="gaussian",
                                treat.vars,
                                control.var="control",
                                bonferroni.npvals=NULL,
                                mult.test.correct=c("treatments","responses"),
                                add.pval.stars=TRUE,
                                model.names=NULL,
                                stderr.func=NULL,
                                digits=3,
                                label=NULL,
                                caption=NULL,
                                col.var.names=NULL,
                                treat.names=NULL,
                                return.tabdf=FALSE,
                                ...){

  #Make sure data.list is a list of named objects
  if(is.list(data.list)==FALSE){
    data.list=list("Data"=data.list)
  }
  if(is.null(names(data.list))==TRUE){ #if no names provided name them Data 1, Data 3,...
    names(data.list)=paste("Data",seq(1,length(data.list)))
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
                 mult.test.correct=mult.test.correct,
                 pivot.xtreat=TRUE,
                 add.pval.stars=TRUE,
                 model.names=model.names,
                 stderr.func=stderr.func,
                 only.control.mean=TRUE)))


  #get the desired columns in the desired order
  cnames=colnames(table.combine)
  colnum.start.treatment=2
  firstcol=c(ifelse("ModelName"%in%cnames,"ModelName","Response"),"ControlMean")
  if(is.null(col.var.names)==FALSE){
    table.combine=table.combine[,col.var.names]
  }
  if(is.null(treat.names)==TRUE){ #if no treat.names supplied, use the variable names
      treat.names=treat.vars
  }


  # name columns to match Table 2b in Blattman et al.
  tabcnames=c("Outcome","Control mean",
              rep(c("ITT","Std. Err.",
                    ifelse(sum(grepl("AdjPvalue",cnames))>0,"Adj. p-value","p-value")),
                  length(treat.vars)),
              "N")

  #make kable table
  table.kab=knitr::kable(table.combine,
                    digits=digits,caption=caption,label=label,
                    col.names=tabcnames,
                    ...)

  #add column header above for each treatment
  treat.header=rep(3,length(treat.vars))
  names(treat.header)=treat.names
  table.kab=kableExtra::add_header_above(table.kab,
                                    header=c(" "=2,
                                      treat.header," "=1),align="c")

  #add overall header above
  table.kab=kableExtra::add_header_above(table.kab,
                                    header=c(" "=2,
                                      "ITT regression"=(3*length(treat.vars))+1),align="c")

  #if more than one dataset pack rows with the label from the data.list names
  if(length(data.list)>1){
    starti=1
  for(df.name in names(data.list)){
    table.kab=kableExtra::pack_rows(table.kab,df.name,starti,starti+max(length(reg.models),length(response.vars))-1)
    starti=starti+max(length(reg.models),length(response.vars))
  }
  }
  if(return.tabdf==TRUE){
    return(list("table.kable"=table.kab,"table.df"=table.combine))
  }
  return(table.kab)
}

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
    stopifnot((is.null(ci.overlap.ref.name)==FALSE)&&(ci.overlap.ref.name%in% unique(table.combine$data)))
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
  if(include.ci.overlap==TRUE){ #if including confidence interval overlap...
    stopifnot((is.null(ci.overlap.ref.name)==FALSE)&&(ci.overlap.ref.name%in% unique(table.combine$data)))
    stopifnot(ci.confidence>0 & ci.confidence<1)
    #z multiplier
    z.half.alpha=stats::qnorm(0.5*(1+ci.confidence))
    #get upper and lower bounds for each confidence interval.
    cis.all=sapply(seq(1,nrow(table.combine)),
                   function(i){
                     table.combine$ITT[i]+
                       (c(-z.half.alpha,z.half.alpha)*table.combine$StdErr[i])})
    table.combine$lower.ci=cis.all[1,]
    table.combine$upper.ci=cis.all[2,]
    #indicator of reference data set rows
    ref.indic=table.combine$data==ci.overlap.ref.name
    ref.tab=table.combine[ref.indic,]
    table.combine$ci.overlap=1 #initialize confidence interval overlap.
    not.ref=unique(table.combine$data)
    not.ref=not.ref[not.ref!=ci.overlap.ref.name]
    for(dta in not.ref){ #for data which is not the reference data...
      dta.indic=table.combine$data==dta #indicator of the dataset
      table.combine$ci.overlap[dta.indic]= #get ci overlap
        CI_overlap(table.combine$lower.ci[dta.indic],table.combine$upper.ci[dta.indic],
                   ref.tab$lower.ci,ref.tab$upper.ci)
    }
  }

  if(pivot.xtreat==TRUE){
    #### deal with mult.test.correct and bonferroni.npvals ###
    pvalname="Pvalue"
    if(is.null(mult.test.correct)==FALSE){
        pvalname="AdjPvalue"
      }else if((is.null(bonferroni.npvals)==FALSE)&&(bonferroni.npvals>1)){
        pvalname="AdjPvalue"
      }
    treat.vars=unique(table.combine$Treatment)
    pivot.values=c(c("ITT","StdErr",pvalname),
                   ifelse(include.ci.overlap==TRUE,c("ci.lower","ci.upper","ci.overlap"),NULL))
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
    temp.colnames=paste0(temp.colnames[,1],temp.colnames[,2])
    other.colnames=colnames(table.combine)[!(colnames(table.combine)%in%temp.colnames)]
    table.combine=table.combine[,c(other.colnames,temp.colnames)]

    #get the desired columns in the desired order
    cnames=colnames(table.combine)
    colnum.start.treatment=2
    firstcol=c(ifelse("ModelName"%in%cnames,"ModelName","Response"),"ControlMean")
    col.order=c(firstcol)
    for(tr.var in treat.vars){
      col.order=c(col.order,paste(pivot.values,tr.var,sep="_"))
    }
    col.order=c(col.order,"nObs")
    table.combine=table.combine[,col.order]
  }
  return(table.combine)
}

