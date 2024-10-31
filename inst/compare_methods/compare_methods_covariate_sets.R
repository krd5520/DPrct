####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

## Considering 3 sets of covariates
sub.covariates=c('age_b', 'asbhostil_b', 'drugssellever_b', 'drinkboozeself_b',
                 'druggrassself_b', 'harddrugsever_b', 'steals_b', 'wealth_indexstd_b',
                 'school_b', 'cognitive_score_b')
red.sub.covariates=c('age_b', 'asbhostil_b', 'drugssellever_b', 'drinkboozeself_b',
                     'druggrassself_b', 'harddrugsever_b', 'steals_b')

baseline.full=c('age_b', 'livepartner_b', 'mpartners_b', 'hhunder15_b', 'famseeoften_b',
                'muslim_b', 'school_b', 'schoolbasin_b', 'literacy_b', 'mathscore_b', 'health_resc_b',
                'disabled_b', 'depression_b', 'distress_b', 'rel_commanders_b', 'faction_b', 'warexper_b',
                'profitsump99avg7d_b', 'wealth_indexstd_b', 'homeless_b', 'slphungry7dx_b', 'savstockp99_b',
                'loan50_b', 'loan300_b', 'illicit7da_zero_b', 'agricul7da_zero_b', 'nonagwage7da_zero_b',
                'allbiz7da_zero_b', 'nonaghigh7da_zero_b', 'agriculeveramt_b', 'nonagbizeveramt_b',
                'nonaghigheveramt_b', 'drugssellever_b', 'drinkboozeself_b', 'druggrassself_b',
                'grassdailyuser_b', 'harddrugsever_b', 'harddrugsdailyuser_b', 'steals_b',
                'stealnb_nonviol_b', 'stealnb_felony_b', 'disputes_all_b', 'asbhostil_b',
                'conscientious_b', 'neurotic_b', 'grit_b', 'rewardresp_b', 'locuscontr_b', 'impulsive_b',
                'selfesteem_b', 'patient_game_real_b', 'inconsistent_game_resc_b', 'risk_game_resc_b',
                'timedecl_b', 'riskdecl_b', 'cognitive_score_b', 'ef_score_b')
########

#list all possible response variables for round==5
response.vars.all = c('fam_asb_lt',
                      paste0(c('drugssellever','stealnb','disputes_all_z','carryweapon',
                               'arrested','asbhostilstd','domabuse_z'),"_ltav"))
treat.vars=c("cashassonly","tpassonly","tpcashass") #treatment variables
block.vars=c("tp_strata_alt","cg_strata") #block variables

#other global paramaters
mv.bins=80 #bins for multivariate histograms
synthdata.budget=c(1,0) #overall budget epsilon, delta
use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables

## read in data
load(paste0(basepath,"/data/SubsetLiberiaWithFactors.Rda"))
list2env(attr(liberia.sub,"column.classifications"),envir=globalenv())
cont.vars.all=c(cont.as.cont,cont.as.cat,response.vars.all)
cont.vars.all=cont.vars.all[!duplicated(cont.vars.all)]
liberia.sub[,colnames(liberia.sub)!="treatment"]=apply(liberia.sub[,colnames(liberia.sub)!="treatment"],2,
                                                       function(x)as.numeric(as.character(x)))
liberia.sub[,block.vars]=apply(liberia.sub[,block.vars],2,
                               function(x)as.factor(as.character(x)))
#########


#### parameters for Full DP Model- Based
# privacy budget variables #
mv.prop.budget=0.02
## within each response, the residual variance gets 10% of the budget
win.y.resvar=0.1
##      the treatment coefficients get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.3
##
# other params #
bound.means=100
bound.sds=c(2^(-15),2^(15))
n.iters=300
alpha=0.05
########

mod.names=c("Anti-social Behaviors, z-score",
            "Usually Sells Drugs",
            "# of thefts/robberies in past 2 weeks",
            "Disputes and fights in past 2 weeks, z-score",
            "Carries a weapon on body",
            "Arrested in past 2 weeks",
            "Aggressive behaviors, z-score",
            "Verbal/physical abuse of partner, z-score")
treat.names=c("Therapy Only","Cash Only","Both")


liberia_compare_methods=function(liberia.sub,
                                 synthdata.budget,mv.bins,use.continuous.noise,
                                 covariate.vars,block.vars,treat.vars,response.vars.all,
                                 cont.as.cont,
                                 mv.prop.budget=NA,win.y.resvar=NA,win.y.trcoef=NA,
                                 bound.means=NULL,bound.sds=NULL,n.iters=NA,range.alpha=0.05,
                                 mod.names,treat.names,
                                 mult.test.correct=c("treatments","responses"),bonferroni.npvals=NA,
                                 ci.confidence=0.95,
                                 rseed=NA,
                                 sim=NA,
                                 include.full.model.based=TRUE){

  start.time=proc.time()
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  if(include.full.model.based==TRUE){
    stopifnot((is.na(mv.prop.budget)==FALSE)&&(mv.prop.budget>0 & mv.prop.budget<1))
    stopifnot((is.na(win.y.resvar)==FALSE)&(is.na(win.y.trcoef)==FALSE))
    stopifnot((win.y.resvar>0)&(win.y.trcoef>0)&((win.y.resvar+win.y.trcoef)<1))
    stopifnot(is.null(bound.means)==FALSE&is.null(bound.sds)==FALSE&is.na(n.iters)==FALSE)

  }

  liberia.sub=liberia.sub[,colnames(liberia.sub)%in%c(treat.vars,block.vars,response.vars.all,covariate.vars,"treatment","control")]

  if(length(synthdata.budget)==1){
    synthdata.budget.eps=synthdata.budget
    synthdata.budget.del=0
  }else{
    synthdata.budget.eps=synthdata.budget[1]
    synthdata.budget.del=synthdata.budget[2]
  }

  ### Constructing regression model formulas for each response variable
  families=rep("gaussian",length(response.vars.all))

  covariate.formula=paste(c(covariate.vars,block.vars),collapse="+")
  treat.expand.formula=paste(treat.vars,collapse="+")

  reg.formulas=paste(response.vars.all,
                     paste(treat.expand.formula,covariate.formula,sep="+"),
                     sep="~")


  #variables
  covariate.data=liberia.sub[,colnames(liberia.sub)%in%c(covariate.vars,block.vars)]
  n.response=length(response.vars.all)
  #### intialize time keeping ###
  comp.times.perturb.hist=data.frame("SyntheticMethod"=c("Perturbed MV Histogram","Not DP Hybrid","Full DP Model-Based"),
                                   "nHistVars"=c(ncol(liberia.sub[,colnames(liberia.sub)!="treatment"]),
                                                 ncol(covariate.data),ncol(covariate.data)),
                                   "SynthTreatmentAssigned"=c(F,T,T),
                                   "FullDP"=c(T,F,T),
                                   "MakeMVHist.Time"=rep(NA,3),
                                   "SanHistProps.Time"=rep(NA,3),
                                   "SampleHist.Time"=rep(NA,3),
                                   "Total.Time"=rep(NA,3))
  comp.times.by.response=data.frame("SyntheticMethod"=c(rep("Not DP Hybrid",2),rep("Full DP Model-Based",4)),
                                    "nPredictorVars"=rep(ncol(covariate.data)+length(block.vars)+length(treat.vars),6),
                                    "Timed.Step"=c("total.time","san.response.time",
                                                   "total.time","iter.proxy.fit","san.summary","san.response"))
  y.times.df=as.data.frame(matrix(rep(NA,n.response*6),ncol=n.response))
  colnames(y.times.df)=response.vars.all
  comp.times.by.response=cbind(comp.times.by.response,y.times.df)

  if(is.na(sim)==FALSE){
    comp.times.perturb.hist$Sim=rep(sim,3)
    comp.times.by.response$Sim=rep(sim,6)
  }
  if(include.full.model.based==FALSE){
    comp.times.perturb.hist=comp.times.perturb.hist[1:2,]
    comp.times.by.response=comp.times.by.response[1:2,]
  }
  get.synths.start=proc.time()
  comp.times.summary=list("Set Up Time"=(get.synths.start-start.time)[[3]])
  print("Starting Synthetic Data Methods")

  ############# Multivariate Histogram Fully Synthetic epsilon=1 ############################
  #get multivariate histogram synthetic data
  mv.out=synthdata_perturb_mvhist(data=liberia.sub[,colnames(liberia.sub)!="treatment"],
                                  epsilon=synthdata.budget.eps,
                                  delta=synthdata.budget.del,
                                  continuous.vars=cont.as.cont[cont.as.cont%in%colnames(covariate.data)],
                                  num.bin=mv.bins,
                                  add.cont.variation=use.continuous.noise,
                                  treatment.colname=treat.vars,
                                  assign.type="block",
                                  return.time=TRUE,
                                  with.treatment=FALSE, #treatment is kept in the mv histogram
                                  within.blocks=FALSE,
                                  blocks=block.vars,
                                  block.sizes=NULL)

  mv.synth.comp.time=attr(mv.out,"comp.time") #computation time
  comp.times.perturb.hist[1,5:8]=mv.synth.comp.time
  comp.times.summary=c(comp.times.summary,
                       list("Perturb MV Histogram"=(proc.time()-get.synths.start)[[3]]))

  print("Done with Perturb MV Histogram")
  ################ Hybrid DP epsilon=1 #######################
  all.y.hybrid.start=proc.time()


  mv.covariates=synthdata_perturb_mvhist(data=covariate.data,
                                         epsilon=synthdata.budget.eps,
                                         delta=synthdata.budget.del,
                                         continuous.vars=cont.as.cont[cont.as.cont%in%colnames(covariate.data)],
                                         num.bin=mv.bins,
                                         add.cont.variation=use.continuous.noise,
                                         treatment.colname="treatment",
                                         assign.type="block",
                                         return.time=TRUE,
                                         with.treatment=TRUE,
                                         within.blocks=TRUE,
                                         blocks=block.vars,
                                         block.sizes=NULL,
                                         conditions=c(treat.vars,"control"))
  for(tr.var in c(treat.vars,"control")){
    mv.covariates[,tr.var]=ifelse(mv.covariates$treatment==tr.var,1,0)
  }
  #liberia.sub[,treat.vars]=apply(liberia.sub[,treat.vars],2,as.factor)
  mv.covariates=mv.covariates[,!(colnames(mv.covariates)%in%c("control"))]
  mv.covariate.comp.time=attr(mv.covariates,"comp.time")





  get_synth_response_hybriddp=function(idx,synth.data=mv.covariates[,!(colnames(mv.covariates)=="treatment")],
                                       formulas=reg.formulas,
                                       y.vars=response.vars.all,
                                       fams=families,
                                       confidential.data=liberia.sub,
                                       tr.vars=treat.vars){
    temp.out=hybrid_synth(formulas[idx],
                          family=fams[idx],
                          confidential.data=confidential.data,
                          epsilon=NA,
                          delta=0,
                          treatment.var=tr.vars,
                          synth.data=synth.data,
                          continuous.vars=NULL,
                          num.bin=NULL,
                          bin.param=NA,
                          add.cont.variation=FALSE,
                          assign.type="simple",
                          blocks=NULL,
                          within.blocks=TRUE,
                          clusters=NULL,
                          returntypes=c("synth.data","comp.time"),
                          rseed=NA)
    response.var=y.vars[idx]
    temp.synth=temp.out[[1]]
    list(temp.synth[,response.var],temp.out[[2]])
  }

  synth.response=parallel::mclapply(seq(1,n.response),get_synth_response_hybriddp)
  hybrid.dp.comp.times=sapply(synth.response,"[[",2)
  colnames(hybrid.dp.comp.times)=paste0(response.vars.all,".comp.time")
  synth.ys=dplyr::bind_cols(lapply(synth.response,"[[",1))
  hybrid.synthdata=cbind(synth.ys,mv.covariates)
  hybrid.synthdata$control=ifelse(hybrid.synthdata$treatment=="control",1,0)
  all.y.hybrid.comp.time=(proc.time()-all.y.hybrid.start)[[3]]

  comp.times.perturb.hist[2,5:8]=mv.covariate.comp.time
  comp.times.summary=c(comp.times.summary,
                       list("Not DP Hybrid"=all.y.hybrid.comp.time))
  comp.times.by.response[1:2,seq(4,3+n.response)]=hybrid.dp.comp.times

  hybrid.dp.comp.times=list("total.hybrid.time"=all.y.hybrid.comp.time,
                            "hybrid.response.times"=hybrid.dp.comp.times,
                            "mv.T.X.time"=mv.covariate.comp.time)

  print("Done with Not-DP Hybrid Method")

  ####################

  if(include.full.model.based==TRUE){
    ################ My Full DP Synthetic epsilon=1 #######################
    synthdata.start=proc.time()


    ################### Set Parameters ##############################
    ##each of the responses get equal proportion of the the overall budget after the sythetic covariate cost
    per.y.prop=(1-mv.prop.budget)/n.response
    ## within each response, the residual variance gets 10% of the budget
    ##      the other coefficients (and the intercept) get the remaining portion of budget within each response
    win.y.x.coefs=(1-win.y.resvar-win.y.trcoef)
    ############

    ##############################
    n.x.coefs=length(covariate.vars)+1+length(c(t(unique(liberia.sub[,block.vars[1]])),
                                                t(unique(liberia.sub[,block.vars[2]]))))
    n.treat.coefs=length(treat.vars)
    win.y.budget.props=c(rep(win.y.trcoef/n.treat.coefs,n.treat.coefs),
                         rep(win.y.x.coefs/n.x.coefs,n.x.coefs),win.y.resvar)
    #get synthetic covariate and treatment data
    #covariate.data from Hybrid DP section
    mv.covariates2=synthdata_perturb_mvhist(data=covariate.data,
                                            epsilon=synthdata.budget.eps,
                                            delta=synthdata.budget.del,
                                            continuous.vars=cont.as.cont[cont.as.cont%in%colnames(covariate.data)],
                                            num.bin=mv.bins,
                                            add.cont.variation=use.continuous.noise,
                                            treatment.colname="treatment",
                                            assign.type="block",
                                            return.time=TRUE,
                                            rseed=3,
                                            with.treatment=TRUE,
                                            within.blocks=TRUE,
                                            blocks=block.vars,
                                            block.sizes=NULL,
                                            conditions=c(treat.vars,"control"))
    for(tr.var in c(treat.vars,"control")){
      mv.covariates2[,tr.var]=ifelse(mv.covariates2$treatment==tr.var,1,0)
    }
    mv.covariates2=mv.covariates2[,!(colnames(mv.covariates2)%in%c("control"))]
    mv.covariate.comp.time2=attr(mv.covariates2,"comp.time")

    get_synth_response_synthdatadp=function(idx,
                                            synth.data=mv.covariates2,
                                            formulas=reg.formulas,
                                            y.vars=response.vars.all,
                                            fams=families,
                                            confidential.data=liberia.sub,
                                            tr.vars=treat.vars,
                                            per.response.eps=per.y.prop*synthdata.budget.eps,
                                            per.response.del=per.y.prop*synthdata.budget.del,
                                            bd.sd=bound.sds,
                                            bd.mean=bound.means,
                                            nits=n.iters,
                                            budget.prop.list=win.y.budget.props,
                                            alph=range.alpha,...

    ){
      temp.out=dp_synthdata(formula=formulas[idx],
                            confidential.data=confidential.data,
                            synth.data=synth.data,
                            synth.epsilon=NULL,
                            synth.delta=NULL,
                            continuous.vars=NULL,
                            num.bin=NULL,
                            bin.param=NA,
                            assign.type="simple",
                            blocks=NULL,
                            clusters=NULL,
                            within.blocks=TRUE,
                            epsilon.list=as.list(per.response.eps*budget.prop.list),
                            delta.list=as.list(per.response.del*budget.prop.list),
                            bd.sd.list=bd.sd,
                            bd.mean.list=bd.mean,
                            alphas.list=alph,
                            num.iters=nits,
                            treatment.var=tr.vars,
                            rseed=NA,
                            return.time=TRUE,
                            return.confidential.table=FALSE,
                            return.san.summary=FALSE,
                            use.san.residerror=TRUE,
                            return.treatment.pvalue=FALSE,
                            return.treatment.CI=FALSE,
                            treat.ppart.list=NULL,
                            ...)

      response.var=y.vars[idx]
      temp.synth=temp.out[[1]]
      list(temp.synth[,response.var],temp.out[[2]])
    }

    synth.response=parallel::mclapply(seq(1,n.response),get_synth_response_synthdatadp)
    full.dp.synth=dplyr::bind_cols(mv.covariates2,dplyr::bind_cols(sapply(synth.response,"[[",1)))
    full.dp.synth$control=ifelse(full.dp.synth$treatment=="control",1,0)
    by.y.comp.times=sapply(synth.response,"[[",2)
    full.dp.comp.times=list("total.dp.time"=(proc.time()-synthdata.start)[[3]],
                            "dp.response.time"=by.y.comp.times,
                            "mv.T.X.time"=mv.covariate.comp.time2)
    comp.times.perturb.hist[3,5:8]=mv.covariate.comp.time2
    comp.times.summary=c(comp.times.summary,
                         list("Full DP Model-Based"=(proc.time()-synthdata.start)[[3]]))
    for(nc in seq(1,n.response)){
      comp.times.by.response[3:6,3+nc]=unlist(by.y.comp.times[,nc])
    }
    print("Done with Full DP Model-Based")
  }

  ##############
  comp.times.summary=c(comp.times.summary,
                       list("All Synthetic Data Methods"=(proc.time()-get.synths.start)[[3]]))

  table.start=proc.time()
  sterr.type="HC1"
  stderr.func=function(mod,std.type=sterr.type){
    sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
  }


  table.data.list=list("Confidential Data"=liberia.sub,
                       "DP Multivariate Histogram Synthetic Data"=mv.out,
                       "Not-DP Hybrid Synthetic Data"=hybrid.synthdata)
  if(include.full.model.based==TRUE){
    table.data.list=c(table.data.list,
                      list( "DP Model-based Synthetic Data Using Sanitized Residual Error"=full.dp.synth))
  }

  out.df=compare_ITTdata(table.data.list,
                         reg.formulas,
                         response.vars=response.vars.all,
                         families=rep("gaussian",length(reg.formulas)),
                         treat.vars=treat.vars,
                         control.var="control",
                         bonferroni.npvals=bonferroni.npvals,
                         mult.test.correct=mult.test.correct,
                         add.pval.stars=TRUE,
                         model.names=mod.names,
                         stderr.func=stderr.func,
                         treat.names=treat.names,
                         include.ci.overlap = TRUE,
                         ci.overlap.ref.name = "Confidential Data",
                         ci.confidence=ci.confidence,
                         pivot.xtreat=TRUE)
  if(is.na(sim)==FALSE){
    out.df$Sim=sim
  }

  end.time=proc.time()
  comp.times.summary=c(comp.times.summary,list("Make Table"=(end.time-table.start)[[3]],
                                               "Total Time"=(end.time-main.start)[[3]])
  )
  attr(out.df,"all.comp.times")=list("Summary"=comp.times.summary,
                                     "Perturb.Hist.Steps"=comp.times.perturb.hist,
                                     "By.Response.Steps"=comp.times.by.response)

  return(list(out.df,comp.times.summary,comp.times.perturb.hist,comp.times.by.response))
}
#save(out.df,comp.times,file=paste0(basepath,"/inst/tables/ITT_compare_methods_data.Rda"))

out1=liberia_compare_methods(liberia.sub=liberia.sub,
                             synthdata.budget=synthdata.budget,
                             mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                             covariate.vars=baseline.full,
                             block.vars=block.vars,treat.vars=treat.vars,response.vars.all=response.vars.all,
                             cont.as.cont=cont.as.cont,
                             mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                             bound.means=bound.means,bound.sds=bound.sds,n.iters=100,range.alpha=0.05,
                             mod.names=mod.names,treat.names=treat.names,
                             mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                             ci.confidence=0.95,
                             rseed=1,
                             sim="Full Baseline Covariates",
                             include.full.model.based=TRUE)

out2=liberia_compare_methods(liberia.sub=liberia.sub,
                             synthdata.budget=synthdata.budget,
                             mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                             covariate.vars=sub.covariates,
                             block.vars=block.vars,treat.vars=treat.vars,response.vars.all=response.vars.all,
                             cont.as.cont=cont.as.cont,
                             mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                             bound.means=bound.means,bound.sds=bound.sds,n.iters=100,range.alpha=0.05,
                             mod.names=mod.names,treat.names=treat.names,
                             mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                             ci.confidence=0.95,
                             rseed=1,
                             sim="Subset of Covariates",
                             include.full.model.based=TRUE)

out3=liberia_compare_methods(liberia.sub=liberia.sub,
                             synthdata.budget=synthdata.budget,
                             mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                             covariate.vars=red.sub.covariates,
                             block.vars=block.vars,treat.vars=treat.vars,response.vars.all=response.vars.all,
                             cont.as.cont=cont.as.cont,
                             mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                             bound.means=bound.means,bound.sds=bound.sds,n.iters=100,range.alpha=0.05,
                             mod.names=mod.names,treat.names=treat.names,
                             mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                             ci.confidence=0.95,
                             rseed=1,
                             sim="Reduced Subset of Covariates",
                             include.full.model.based=FALSE)

full.comparison=dplyr::bind_rows(out1[[1]],out2[[1]],out3[[1]])
summary.times=dplyr::bind_rows(unlist(out1[[2]]),unlist(out2[[2]]),unlist(out3[[2]]))
summary.times$CovariateSet=c("Full Baseline","Subset","Reduced Subset")
summary.times$nPredictors=c(length(block.vars)+length(treat.vars)+length(baseline.full),
                            length(block.vars)+length(treat.vars)+length(sub.covariates),
                            length(block.vars)+length(treat.vars)+length(red.sub.covariates))
perturb.hist.times=dplyr::bind_rows(out1[[3]],out2[[3]],out3[[3]])
by.response.times=dplyr::bind_rows(out1[[4]],out2[[4]],out3[[4]])

save(full.comparison,summary.times,perturb.hist.times,by.response.times,
     file=paste0(basepath,"/inst/tables/ITT_compare_methods_data.Rda"))


