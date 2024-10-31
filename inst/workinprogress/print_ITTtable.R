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
