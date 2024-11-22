

#' Hybrid-DP creates a synthetic dataset for RCTs using DP-based methods
#' @description (Algorithm 1 in paper) which uses estimates of coefficients,
#'    treatment effect, and standard error from the confidential data to create
#'    a synthetic dataset and estimated treatment effect based on differential
#'    privacy methods.
#'
#' @param formula either string or \code{\link{as.formula}} object
#'    with response variable as the first variable in the formula.
#'    If \code{treatment.var} is not supplied, the second variable in the formula
#'    is taken as the treatment variable.
#' @param family a description of the error distribution and link function to be
#'    supplied to \code{\link{glm}}.
#' @param confidential.data a data.frame object with the confidential data to be used.
#' @inheritDotParams stats::glm
#' @inheritParams synthdata_perturb_mvhist
#' @inheritParams treatment_assign
#' @param rseed a value to set the randomization seed. If \code{rseed=NA} (default),
#'    then no random seed it set.
#' @param treatment.var a vector of column indices, logical, or strings of column
#'    names to denote what columns of \code{confidential.data} are treatment columns
#' @param returntypes a vector of strings where the elements can take values
#'    \code{"synth.data"}, \code{"san.treat.effect"}, \code{"san.model"},
#'    and/or \code{"confidential.model"}. These determine what values are returned
#'    from the function.
#' @returns a list or element determined by \code{returnvalues}. The elements are
#'    returned in a list in the same order as the vector requests them. (See details)
#' @details
#' For the elements of \code{returntypes}:
#' If "synth.data" is included, the full synthetic dataset is returned.
#' If "san.treat.effect", the sanitized (or privacy-preserving) treatment effect
#' from fitting a model with the synthetic data is returned.
#' If "san.model" is included, the \code{\link{glm}} object
#' of the model fitted with the synthetic data is outputted.
#' If "confidential.model" is included, the \code{\link{glm}} object
#' of the model fitted with the confidential data is outputted without privacy protections.
#' If "comp.time" is included, the computation time of the algorithm is outputted.
#'
#' Additional parameters are passed to the \code{\link{glm}} function.
#' To specify random treatment assignment further. Use \code{synthdata_perturb_mvhist()}
#' with \code{with.treatment=TRUE} to create your synthetic data first.
#' Users may bring their own synthetic data and use \code{treatment_assign} to
#' further control the specifications of the random treatment assignment methodology.
#' Then the synthetic dataset can be supplied to this function as \code{synth.data} input.
#'
#' @importFrom stats glm as.formula coefficients
#'
#' @export
#'
hybrid_synth<-function(formula,
                    family="gaussian",
                    confidential.data,
                    epsilon=NA,
                    delta=0,
                    treatment.var=NULL,
                    synth.data=NULL,
                    continuous.vars=NULL,
                    num.bin=NULL,
                    bin.param=NA,
                    add.cont.variation=FALSE,
                    assign.type="simple",
                    blocks=NULL,
                    within.blocks=TRUE,
                    clusters=NULL,
                    returntypes=c("synth.data","san.model"),
                    rseed=NA,...){
  start.time=proc.time()
  return.time="comp.time" %in% returntypes
  time.values=NULL
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  if(is.null(synth.data)==TRUE){
    #if not synthetic data supplied, a positive privacy cost needs to be supplied
    stopifnot(is.numeric(epsilon)&&epsilon>0)
  }

  psf=parse_formula(formula=formula,confidential.data = confidential.data,
                    treatment.var = treatment.var)
  formula=psf[[1]]
  model.vars=psf[[2]]
  response.var=psf[[3]]
  synth.vars=psf[[4]]
  treatment.var=psf[[5]]

  model.fit.start=proc.time()
  #### Step 0: In paper
  #fit the confidential data to the model and get estimated coefficient and standard error
  conf.model<-stats::glm(formula=formula,family=family,data=confidential.data,...)
  ####
  model.fit.stop=proc.time()


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
                                      epsilon=epsilon,
                                      delta=delta,
                                      continuous.vars=continuous.vars,
                                      num.bin=num.bin,
                                      bin.param=bin.param,
                                      add.cont.variation=add.cont.variation,
                                      assign.type=assign.type,
                                      blocks=blocks,
                                      within.blocks=within.blocks,
                                      clusters=clusters,
                                      return.time=return.time)
  if(return.time==TRUE){
    time.values=c(time.values,
                  list("synthetic.T.X.time"=attr(synth.data,"comp.time")))
  }
  }

  warning("right before simulate_response_glm")
  san.y.start=proc.time()
  #simulate response from synthetic data and confidential model
  synth.model.data=synth.data[,colnames(synth.data)%in%model.vars]
  synth.data$response=simulate_response_glm(conf.model,newdata=synth.model.data)

  #if response variable is factor, replace values with level names
  if(is.factor(confidential.data[,response.var])==TRUE){
    synth.data$response=as.factor(
      levels(confidential.data[,response.var])[1+synth.data$response])
  }
  if(return.time==TRUE){
    san.y.stop=proc.time()
    time.values=c(time.values,list("fit.model.time"=(model.fit.stop-model.fit.start)[[3]],
                                "san.response.time"=(san.y.stop-san.y.start)[[3]]))
  }


  warning("right after simulate_response_glm")
  #make response name match confidential data
  n.synth.cols=length(colnames(synth.data))
  colnames(synth.data)<-c(colnames(synth.data)[seq(1,n.synth.cols-1)],
                                response.var)

  #make column order match that of confidential data
  conf.colnames=colnames(confidential.data)
  col.order=conf.colnames%in% colnames(synth.data)
  synth.data=synth.data[,conf.colnames[col.order]]

  rownames(synth.data)=NULL

  san.model=NULL
  if(sum(c("san.model","san.treatment.effect")%in%returntypes)>0){
    san.model=stats::glm(formula=formula,family=family,data=synth.data,...)
  }

  ## PARSE returnvalues
  output=output_returntypes(returntypes = returntypes,
                            synth.data = synth.data,
                            san.model = san.model,
                            treatment.var = treatment.var,
                            confidential.model=conf.model,
                            start.time=start.time,
                            time.values=time.values,
                            synth.priv.cost = c("epsilon"=epsilon,"delta"=delta))

  return(output)
}



