
#' Internal-use within handle_synthetic_data to separate treatment columns
#' @param treat.col is a vector of strings which are the treatment column name
#'    followed by the level for each treatment variable separated by an "_".
#' @param treat.var is a vector of strings which are the treatment column names.
#' @returns a dataframe of the treatment values.
#' @keywords internal
#'
#' @noRd
seperate_multiple_treatment=function(treatment.col,treat.var){
  #split each treatment string at the "_" and place the pieces in separate columns
  treat.mat=matrix(
    unlist(
      strsplit(treatment.col,"_")),nrow=length(treatment.col),byrow=T)
  #remove column name from each string piece
  treat.mat=sapply(seq(1,ncol(treat.mat)),
                         function(x)gsub(treat.var[x],"",treat.mat[,x]))
  #make a data.frame and name the columns after the treatment.var
  treat.df=as.data.frame(treat.mat)
  colnames(treat.df)=treat.var
  return(treat.df)
}



#' Internal-function to parse formula for algorithms to extract variable names
#' @param formula either string or \code{\link{as.formula}} object with response variable as the first variable in the formula. If \code{treatment.var} is not supplied, the second variable in the formula is taken as the treatment variable.
#' @param confidential.data a data.frame object with the confidential data to be used.
#' @inheritParams synthdata_perturb_mvhist
#' @inheritParams treatment_assign
#' @returns list with unified strings as var names and public data.frame
#'
#' @keywords internal
#' @importFrom stats as.formula
#' @noRd

parse_formula<-function(formula,confidential.data,treatment.var=NULL){
  if(is.character(formula)==TRUE){ #if needed, formula changed from string to formula
    formula=stats::as.formula(formula)
  }
  model.vars=all.vars(formula)

  #### Dissect model formula for response, treatment and covariates ####
  #get names of response, treatment, and covariate variables from model
  response.var=model.vars[1]
  conf.colnames=colnames(confidential.data)
  if("."%in%model.vars){ #if formula is response~. use all columns
    model.vars=conf.colnames
  }
  predictor.vars=model.vars[model.vars!=response.var]
  n.model.vars=length(model.vars)

  #check that variables in the formula are in the confidential data
  if(sum(model.vars %in% colnames(confidential.data))!=n.model.vars){
    stop("Variables in formula are not in the confidential data.")
  }

  treatment.var=parse_treatment_var(treatment.var=treatment.var,
                                   confidential.data = confidential.data,
                                   predictor.vars = predictor.vars)




  #covariates to be made synthetic is necessary
  synth.vars<-predictor.vars[!(predictor.vars%in%c(treatment.var))]
  treat.in.form.idx=which(predictor.vars%in%treatment.var)
  if(!identical(treat.in.form.idx,seq(1,max(treat.in.form.idx)))){
    predictor.form=gsub("-","+",as.character(formula)[3])
    split.pred.form=strsplit(predictor.form,"+")
    treat.terms.idx=grepl(paste(treatment.var,collpase="|"),split.pred.form)
    new.form=paste(paste0(as.character(formula)[1:2],collapse=""),
                         paste(split.pred.form[treat.terms.idx],collapse="+"),"+",
                         paste(split.pred.form[!treat.terms.idx],collpase="+"))
    warning(
      paste("Treatment variable(s) should be the first predictor variable(s) in the formula. The formula is reorganized to be:",
                  new.form))
    formula=stats::as.formula(new.form)
  }





  #variables to be made synthetic are any that are not, treatment, response, or public
  synth.vars<-conf.colnames[(!(conf.colnames%in%response.var))&
                              (!(conf.colnames%in%treatment.var))]

  return(list("formula"=formula,
                    "model.vars"=model.vars,
                    "response.var"=response.var,
                    "synth.vars"=synth.vars,
                    "treatment.var"=treatment.var))

}


#' Internal-use to extract names of treatment variable(s)
#' @inheritParams parse_formula
#' @param check.pred.vars is a logical to indicate if the function should check
#'    if the treatment variable is in the predictor variables.
#' @returns string or vector of strings that is the names of the treatment columns
#'
#' @keywords internal
#' @noRd
parse_treatment_var=function(treatment.var,confidential.data,predictor.vars=NULL,check.pred.vars=TRUE){
  if(is.null(treatment.var)==T){ #if treatment.var not supplied, first predictor used
    treatment.var<-predictor.vars[1]
    message(paste0("No treatment.var input. Using '",treatment.var,"' as treatment variable."))
    treatment.col<-confidential.data[,treatment.var]
    #check this first predictor is factor or string and has no NA values
    #if(((is.character(treatment.col)+is.factor(treatment.col))==0)){
    #  stop(paste("No specified treatment variable.",
    #                   "First predictor in formula cannot be used as treatment.",
    #                   "It is not a column of strings or factor type."))
    #}
  }else{ #treatment variable is specified
    if(is.character(treatment.var)==F){
      treatment.var=colnames(confidential.data)[treatment.var]
    }

    if((check.pred.vars==TRUE)&&
       sum(treatment.var %in% predictor.vars)!=length(treatment.var)){
      stop(paste0("Specified treatment variable(s) [ ",
                        paste(treatment.var,collapse=", "),
                        "] is not in the provided formula."))
    }
  }
  return(treatment.var)
}


