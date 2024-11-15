#' Internal-use function to check if synthetic data is inputted is of the expected format. If necessary, the synethtic data is created with multivariate perturbed histogram and/or the treatment is randomly assigned.
#' @inheritParams synthdata_perturb_mvhist
#' @returns data.frame of the synthetic data with the randomly assigned treatments.
#' @keywords internal
#'
#' @noRd
handle_synthetic_data<-function(confidential.data,
                                synth.data=NULL,
                                treatment.var,
                                synth.vars,
                                model.vars,
                                epsilon=NA,
                                delta=0,
                                continuous.vars=NULL,
                                num.bin=NULL,
                                bin.param=NA,
                                num.arms=NULL,
                                add.cont.variation=FALSE,
                                assign.type="simple",
                                within.blocks=TRUE,
                                blocks=NULL,clusters=NULL,return.time=TRUE){
  ##### Handling synthetic data ####
  start.time=proc.time()
  time.values=NULL
  if(is.null(synth.data)==TRUE){ #if no synthetic data provided, make synthetic data
    delta=ifelse(is.na(delta)==TRUE,0,delta)
    stopifnot(is.numeric(epsilon)==TRUE&&epsilon>0&delta>=0)

    subset.conf.data<-data.frame(confidential.data[,synth.vars])

    if(ncol(subset.conf.data)==1){
      if(is.character(synth.vars)){
        colnames(subset.conf.data)=synth.vars
      }else{
        colnames(subset.conf.data)=colnames(confidential.data)[synth.vars]
      }
    }


    #if there is more than one treatment column
    treat.parsed=parse_treatment_var(treatment.var,confidential.data,
                                     predictor.vars=synth.vars,check.pred.vars = FALSE)
    treatment.var=treat.parsed[1]

    if(length(treatment.var)>1){
      treat.combos=expand.grid( #get all combinations of treatment levels
        sapply(
          treatment.var,
          function(x){
            paste0(x,levels(as.factor(confidential.data[,x])))
          })
      )
      #name the conditions
      treat.conditions=apply(treat.combos,1,paste,collapse="_")
      treatment.colname="treatment"
    }else{ #one treatment column
      treat.conditions=levels(as.factor(confidential.data[,treatment.var]))
      treatment.colname=treatment.var
    }
    synth.data<-synthdata_perturb_mvhist(data=subset.conf.data,
                                         epsilon=epsilon,
                                         delta=delta,
                                         continuous.vars=continuous.vars,
                                         num.bin=num.bin,
                                         bin.param=bin.param,
                                         add.cont.variation=add.cont.variation,
                                         with.treatment=TRUE,
                                         assign.type=assign.type,
                                         treatment.colname=treatment.colname,
                                         blocks=blocks,
                                         clusters=clusters,
                                         within.blocks=within.blocks,
                                         conditions=treat.conditions,
                                         return.time=return.time)
    time.values=c(time.values,attr(synth.data,"comp.time"))
  }else{ #if synthetic data provided
    check.synth.start=proc.time()
    if(sum(!(treatment.var%in%colnames(synth.data)))!=0){#treatment.var not in synth data
      #if there is more than one treatment column
      treatment.var=colnames(confidential.data[,treatment.var])
      if(length(treatment.var)>1){
        treat.combos=expand.grid( #get all combinations of treatment levels
          sapply(
            treatment.var,
            function(x){
              paste0(x,levels(as.factor(confidential.data[,x])))
            })
        )
        #name the conditions
        treat.conditions=apply(treat.combos,1,paste,collapse="_")
        treatment.colname="treatment"
      }else{
        treat.conditions=levels(as.factor(confidential.data[,treatment.var]))
        treatment.colname=treatment.var
      }

      #add treatment variable
      synth.data<-treatment_assign(synth.data=synth.data,
                                   assign.type=assign.type,
                                   treatment.colname=treatment.colname,
                                   blocks=blocks,
                                   clusters=clusters,
                                   conditions=treat.conditions)
      assign.treat.stop=proc.time()
      time.values=c(time.values,c("assign.treat.time"=(check.synth.start-assign.treat.stop)[[3]]))

    }
  }

  if(length(treatment.var)>1){
    treatment.df=seperate_multiple_treatment(synth.data$treatment,treat.var = treatment.var)
    synth.data=cbind(synth.data[,-which(colnames(synth.data)=="treatment")],treatment.df)
  }

  #if the covariates and treatment is not in the synthetic data sent error
  if(sum(c(model.vars) %in% colnames(synth.data))!=length(model.vars)-1){
    stop(paste0("Treatment and covariate variables are not in the synthetic data.",
                "Synthetic data variables include:",
                paste(colnames(synth.data),collapse=", ")))
  }
  if(sum(!(colnames(synth.data)%in%continuous.vars))>0){
  synth.data[,!(colnames(synth.data)%in%continuous.vars)]=apply(synth.data[,!(colnames(synth.data)%in%continuous.vars)],2,as.factor)
  }
  if(return.time==TRUE){
    synthstop=proc.time()
    time.values=c(time.values,c("make.or.check.synthdata.time"=(synthstop-start.time)[[3]]))
    attr(synth.data,"comp.time")=time.values
  }
  return(synth.data)
  #### End Handling Synthetic Data ####
}
