####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()


#library(DPrtc)
webshot::install_phantomjs() #for saving kable as image
#treat everything with gaussian and continuous
###################### Global for Synthetic Data Variables ###################
#
# #number of response variables (and how many are binary responses)
# n.response=8 #fam_asb_lt, drugssellever_ltav, stealnb_ltav, disputes_all_z_ltav,
# #             #carryweapon_ltav, arrested_ltav, asbhostilstd_ltav, domabuse_z_ltav
# n.bin.response=0 #"drugssellever_ltav"   "carryweapon_ltav"     "arrested_ltav"

#number of covariates and number of levels for categorical covariates
########## Global for Model Variables ############
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

#list all possible response variables for round==5
response.vars.all = c('fam_asb_lt',
                      paste0(c('drugssellever','stealnb','disputes_all_z','carryweapon',
                               'arrested','asbhostilstd','domabuse_z'),"_ltav"))
treat.vars=c("cashassonly","tpassonly","tpcashass") #treatment variables
block.vars=c("tp_strata_alt","cg_strata") #block variables

#other global paramaters
mv.bins=80 #bins for multivariate histograms
synthdata.budget.eps=1 #overall budget epsilon
synthdata.budget.del=0 #overall budget delta
use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables

#robust standard errors
sterr.type="HC1"
sterr.func=function(mod,std.type=sterr.type){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}
########


## read in data
load(paste0(basepath,"/data/SubsetLiberiaWithFactors.Rda"))
list2env(attr(liberia.sub,"column.classifications"),envir=globalenv())
cont.vars.all=c(cont.as.cont,cont.as.cat,response.vars.all)
cont.vars.all=cont.vars.all[!duplicated(cont.vars.all)]
liberia.sub[,block.vars]=apply(liberia.sub[,block.vars],2,
                               function(x)as.factor(as.character(x)))
response.vars.all=response.vars.all[1:7]



# ######### Model Coefficients & Residual Error
# ## to generate continuous covariate variables
# # we use a normal distibution with a normal prior distribution for the mean
# # and a uniform prior distribution for the standard deviation
# cont.cov.prior.params=list("mu"=10,"sigma"=5,"lw.bd"=0,"up.bd"=30)
# ##
#
# ## to generate categorical covariate variables
# # we use probabilities with a dirichlet distribution prior whose shape parameter
# # is sampled with replacement from the number of levels* a multiplier
# cat.cov.prior.multiplier=3
#
# # we generate covariate coefficients with a normal prior distribution.
# # The residual error is generated with a uniform prior distribution
# response.prior.params=list("mu"=1,"sigma"=1.5,"lw.bd"=10,"up.bd"=39)
#

sim.seed=1 #random seed to create simulated data

##########################


############ Global Parameters for Synthetic Data Generation ###############
#other global paramaters
mv.bins=80 #bins for multivariate histograms
synthdata.budget.eps=1 #overall budget epsilon
synthdata.budget.del=0 #overall budget delta
use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables
########################

#robust standard errors
sterr.type="HC1"
sterr.func=function(mod,std.type=sterr.type){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}
########


### Constructing regression model formulas for each response variable
families=rep("gaussian",length(response.vars.all))

base.formula=paste(c("treatment",baseline.full,block.vars),collapse="+")
reg.formulas=paste(response.vars.all,
                   base.formula,
                   sep="~")

#get matrix with an estimated coefficient values (coef.mat) and residual error variance (res.err.var)
# for each response variable
#coef.mat: rows= estimated coefficients, cols=response variables
get_coefs_reserrvar=function(conf.data=liberia.sub,
                             response.vars=response.vars.all,
                             base.models=base.formula){
  coef.mat=NULL
  res.err.var=NULL
  last.col=NULL
  for(i in seq(1,length(response.vars))){
    mod.summary=summary(lm(paste0(response.vars[i],"~",base.models),data=conf.data))
    est.coefs=mod.summary$coefficients[,1]
    if(is.null(last.col)){
      last.col=names(est.coefs)
    }
    print(length(est.coefs))
    print(sum(is.na(est.coefs)))
    print(dim(coef.mat))
    coef.mat=dplyr::bind_cols(coef.mat,est.coefs))
    res.err.var=c(res.err.var,mod.summary$sigma^2)

  }
  coef.mat=t(coef.mat)
  colnames(coef.mat)=response.vars
 return(list("coefs.mat"=tcoef.mat,"res.err.var"=res.err.var))
}

get_coefs_reserrvar(response.vars=response.vars.all[1:7])
########## Creating simulated data ##################
#coef.mat needs to have col for each response and row for each coefficient
simdata_origpredictor=function(predictor.data,
                                sim.seed=NA,
                        n.response=NA,
                        response.colnames=NULL,
                        response.prior.params=NULL,
                        coefs.mat=NULL,
                        res.err.var=NULL,
                        return.time=TRUE,
                        pred.colnames=NULL,
                        pred.cat=NULL,
                        treat.colname="treatment"){

  if(return.time==TRUE){ #if returning computation time, start timer
    start.time=proc.time()
  }
  stopifnot((is.na(n.response)==FALSE)|(is.null(response.colnames)==FALSE))
  stopifnot((is.null(response.prior.params)==FALSE)|((is.null(coefs.mat)==FALSE)&(is.null(res.err.var)==FALSE)))

   if(is.na(sim.seed)==FALSE){ #is random seed supplied
    set.seed(sim.seed)
  }else{
    message("sim.seed==NA. No random seed set within this function.
            Random seeds help replicability.")
  }
  if(is.data.frame(predictor.data)==FALSE){ #make sure predictor data is data.frame
  predictor.data=try(as.data.frame(predictor.data))
  }

  if(is.null(pred.colnames)==FALSE){ #if predictor colnames are supplied.
    #if there is a predictor column name (in pred.colnames) for each column, set column names
    #if the column names (in pred.colnames) are a subset of colnames(predictor.data), select those columns
    if(length(pred.colnames)==ncol(predictor.data)){
      colnames(predictor.data)=pred.colnames
    }else if(sum(pred.colnames %in% colnames(predictor.data))==length(pred.colnames)){
      predictor.data=predictor.data[,pred.colnames,drop=FALSE]
    }else{
      warning(paste("pred.colnames is not a vector of same length as ncol(predictor.data)
              and not a vector of column names already present in predictor.data. Instead
              the existing (or default) column names of predictor.data is used:",
                    paste0(colnames(predictor.data)[seq(1,min(3,ncol(predictor.data)))],collapse=", "),
                    "..."))
    }
  }
  if(is.null(pred.cat)==FALSE){
    predictor.data[,pred.cat]=apply(predictor.data[,pred.cat],2,as.factor)
  }

if(is.na(n.response)==TRUE){
  n.response=length(response.colnames)
  predictor.data=predictor.data[,!(colnames(sim.data)%in%response.colnames)]
}else{
  response.colnames=paste0("y",seq(1,n.response))
}
sim.data=predictor.data

if(length(treat.colname)>1){ #if more than 1 treatment column supplied (combine them)
  tr.df=sim.data[,treat.colname,drop=FALSE]
  sim.data$treatment=rep(colnames(tr.df),nrow(tr.df))[as.logical(t(as.matrix(tr.df)))]
  sim.data$treatment[is.na(sim.data$treatment)]=="control"
  pred.data=sim.data[,!(colnames(sim.data)%in%treat.colname)]
  treat.colname="treatment"
}else if(sum(unique(na.omit(sim.data[,treat.colname]))%in%colnames(sim.data))>1){
  #if values of the treatment column are also column names, remove them
  pred.data=sim.data[,!(colnames(sim.data)%in%unique(na.omit(sim.data[,treat.colname])))]
}else{
  pred.data=sim.data
}


  ## generate response variables
  predictor.model=paste0("~",paste0(colnames(pred.data),collapse="+"))

  #predictor matrix for model
  modMat=stats::model.matrix(as.formula(predictor.model),data=pred.data)
  n.coefs=ncol(modMat)
  out.coefs.mat=matrix(NA,nrow=n.coefs+1,ncol=n.response,
                           dimnames = list(c(colnames(modMat),"residual.error.var"),response.vars.all))

  #for each response variable to simulate
  #generate model coefficients and residual error from response.prior.params
  for(i in seq(1,n.response)){
    if(is.null(coefs)==TRUE){
      coefs=stats::rnorm(n.coefs,mean=response.prior.params[[1]],sd=response.prior.params[[2]])
    }else{
      coefs=coefs.mat[,i]
    }
    if(is.null(res.err.var)==TRUE){
      res.err=stats::runif(1,min=response.prior.params[[3]],max=response.prior.params[[4]])
    }else{
      res.err=sqrt(res.err.var[i])
    }
      fit.response=tcrossprod(modMat,t(coefs))
      model.sim=stats::rnorm(n.obs,mean=fit.response,sd=res.err)
      sim.data[,response.var.all[i]]=model.sim
      out.coefs.mat[,i]=c(coefs,res.err^2)
    }

  out.list=list("sim.data"=sim.data,
                "predictor.model"=predictor.model,
                "coefficient.matrix"=out.coefs.mat)
  if(return.time==TRUE){
    out.list=c(out.list,list("sim.data.time"=(start.time-proc.time())[[3]]))
  }
  return(out.list)
}

temp.out=get_coefs_reserrvar()
list2env(simdata_origpredictors(predictor.data=liberia.sub[,c("treatment",baseline.full,block.vars)],
                        sim.seed=sim.seed,
                        response.colnames=response.vars.all,
                        coefs.mat=temp.out$coefs.mat,
                        res.err.var=temp.out$res.err.var,
                        return.time=TRUE,
                        pred.cat=cat.vars[!(cat.vars%in%response.vars.all)],
                        return.time=TRUE),globalenv())



save(sim.data,file=paste0(basepath,"/data/simulated_data_v2.rda"))


### Constructing regression model formulas for each response variable
reg.formulas=paste(response.vars.all,predictor.model)

cnames=colnames(sim.data)
cat.vars=c(cnames[sapply(cnames,function(x)is.factor(sim.data[,x]))],
           cnames[grepl("Indicator",cnames)])
cont.as.cont=cnames[!(cnames %in% cat.vars)]
treat.vars=cnames[(grepl("treatment",cnames))&(cnames!="treatment")]
block.vars=cnames[grepl("block",cnames)]

ITTtable(sim.data,reg.formulas,response.vars.all,families=families,treat.vars = treat.vars)

print("Start Synthetic Data Generation")

get.synths.start=proc.time()
############# Multivariate Histogram Fully Synthetic epsilon=1 ############################
#get multivariate histogram synthetic data
mv.out=synthdata_perturb_mvhist(data=sim.data[,colnames(sim.data)!="treatment"],
                                epsilon=synthdata.budget.eps,
                                delta=synthdata.budget.del,
                                continuous.vars=cont.as.cont,
                                num.bin=mv.bins,
                                add.cont.variation=use.continuous.noise,
                                treatment.colname=treat.vars,
                                assign.type="block",
                                return.time=TRUE,
                                rseed=1,
                                with.treatment=FALSE, #treatment is kept in the mv histogram
                                within.blocks=FALSE,
                                blocks=block.vars,
                                block.sizes=NULL)

mv.synth.comp.time=attr(mv.out,"comp.time") #computation time

print("Multivariate Histogram Done")
################ Hybrid DP epsilon=1 #######################
all.y.hybrid.start=proc.time()

#variables
n.response=length(response.vars.all)

covariate.data=sim.data[,!(colnames(sim.data)%in%c(response.vars.all,treat.vars,"control","treatment"))]
mv.covariates=synthdata_perturb_mvhist(data=covariate.data,
                                       epsilon=synthdata.budget.eps,
                                       delta=synthdata.budget.del,
                                       continuous.vars=cont.as.cont[cont.as.cont%in%colnames(covariate.data)],
                                       num.bin=mv.bins,
                                       add.cont.variation=use.continuous.noise,
                                       treatment.colname="treatment",
                                       assign.type="block",
                                       return.time=TRUE,
                                       rseed=2,
                                       with.treatment=TRUE,
                                       within.blocks=TRUE,
                                       blocks=block.vars,
                                       block.sizes=NULL,
                                       conditions=c(treat.vars,"control"))
for(tr.var in c(treat.vars,"control")){
  mv.covariates[,tr.var]=ifelse(mv.covariates$treatment==tr.var,1,0)
}
#sim.data[,treat.vars]=apply(sim.data[,treat.vars],2,as.factor)
mv.covariates=mv.covariates[,!(colnames(mv.covariates)%in%c("control"))]
mv.covariate.comp.time=attr(mv.covariates,"comp.time")





get_synth_response_hybriddp=function(idx,synth.data=mv.covariates[,!(colnames(mv.covariates)=="treatment")],
                                     formulas=reg.formulas,
                                     y.vars=response.vars.all,
                                     fams=families,
                                     confidential.data=sim.data,
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
colnames(synth.ys)=response.vars.all
hybrid.synthdata=cbind(synth.ys,mv.covariates)
hybrid.synthdata$control=ifelse(hybrid.synthdata$treatment=="control",1,0)
all.y.hybrid.comp.time=(proc.time()-all.y.hybrid.start)[[3]]
hybrid.dp.comp.times=list("total.hybrid.time"=all.y.hybrid.comp.time,
                          "hybrid.response.times"=hybrid.dp.comp.times,
                          "mv.T.X.time"=mv.covariate.comp.time)

print("Hybrid Synthetic Done")
####################

################ My Full DP Synthetic epsilon=1 #######################
synthdata.start=proc.time()


################### Set Parameters ##############################
#### privacy budget variables ####
mv.prop.budget=0.02
##each of the 8 responses get 12.25% of the overall budget (i.e. 0.98/8)
per.y.prop=(1-mv.prop.budget)/n.response
## within each response, the residual variance gets 10% of the budget
win.y.resvar=0.1
##      the treatment coefficients get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.3
##      the other coefficients (and the intercept) get the remaining 60%
win.y.x.coefs=(1-win.y.resvar-win.y.trcoef)
############

### other params #####
bound.means=100
bound.sds=c(2^(-15),2^(15))
n.iters=300
alpha=0.05
########

##############################
n.cat.cov=length(n.cat.covariate.levels)
n.x.coefs=n.covariates-(2*n.cat.cov)+sum(n.cat.covariate.levels)+1+
  length(c(t(unique(sim.data[,block.vars[1]])),
           t(unique(sim.data[,block.vars[2]]))))
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
                                        confidential.data=sim.data,
                                        tr.vars=treat.vars,
                                        per.response.eps=per.y.prop*synthdata.budget.eps,
                                        per.response.del=per.y.prop*synthdata.budget.del,
                                        bd.sd=bound.sds,
                                        bd.mean=bound.means,
                                        nits=n.iters,
                                        budget.prop.list=win.y.budget.props,
                                        alph=alpha,...

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
synth.ys=dplyr::bind_cols(sapply(synth.response,"[[",1))
colnames(synth.ys)=response.vars.all
full.dp.synth=dplyr::bind_cols(mv.covariates2,synth.ys)
full.dp.synth$control=ifelse(full.dp.synth$treatment=="control",1,0)
full.dp.comp.times=list("total.dp.time"=(proc.time()-synthdata.start)[[3]],
                        "dp.response.time"=sapply(synth.response,"[[",2),
                        "mv.T.X.time"=mv.covariate.comp.time2)


print("DP Model-Based Synthetic Done")
# ################ Thier Full DP Synthetic epsilon=1 #######################
# synthdata.start2=proc.time()
#
#
# ################### Set Parameters ##############################
# #### privacy budget variables ####
# mv.prop.budget=0.02
# ##each of the 8 responses get 12.25% of the overall budget (i.e. 0.98/8)
# per.y.prop=(1-mv.prop.budget)/n.response
# ##      the treatment coefficients get 30% of the budget (i.e. 10% for each)
# win.y.trcoef2=0.3
# ##      the other coefficients (and the intercept) get the remaining 60%
# win.y.x.coefs2=(1-win.y.trcoef)
# ############
#
# ### other params #####
# bound.means=100
# bound.sds=c(2^(-15),2^(15))
# n.iters=300
# alpha=0.05
# ########
#
# ##############################
# n.x.coefs=length(baseline.full)+1+length(c(t(unique(sim.data[,block.vars[1]])),
#                                            t(unique(sim.data[,block.vars[2]]))))
# n.treat.coefs=length(treat.vars)
# win.y.budget.props2=c(rep(win.y.trcoef2/n.treat.coefs,n.treat.coefs),
#                      rep(win.y.x.coefs2/n.x.coefs,n.x.coefs))
#
# #get synthetic covariate and treatment data
# #covariate.data from Hybrid DP section
# mv.covariates3=synthdata_perturb_mvhist(data=covariate.data,
#                                         epsilon=synthdata.budget.eps,
#                                         delta=synthdata.budget.del,
#                                         continuous.vars=cont.as.cont[cont.as.cont%in%colnames(covariate.data)],
#                                         num.bin=mv.bins,
#                                         add.cont.variation=use.continuous.noise,
#                                         treatment.colname="treatment",
#                                         assign.type="block",
#                                         return.time=TRUE,
#                                         rseed=4,
#                                         with.treatment=TRUE,
#                                         within.blocks=TRUE,
#                                         blocks=block.vars,
#                                         block.sizes=NULL,
#                                         conditions=c(treat.vars,"control"))
# for(tr.var in c(treat.vars,"control")){
#   mv.covariates3[,tr.var]=ifelse(mv.covariates3$treatment==tr.var,1,0)
# }
# mv.covariates3=mv.covariates3[,!(colnames(mv.covariates3)%in%c("control"))]
# mv.covariate.comp.time3=attr(mv.covariates3,"comp.time")
#
# get_synth_response_synthdatadp_noresiderr=function(idx,
#                                         synth.data=mv.covariates3,
#                                         formulas=reg.formulas,
#                                         y.vars=response.vars.all,
#                                         fams=families,
#                                         confidential.data=sim.data,
#                                         tr.vars=treat.vars,
#                                         per.response.eps=per.y.prop*synthdata.budget.eps,
#                                         per.response.del=per.y.prop*synthdata.budget.del,
#                                         bd.sd=bound.sds,
#                                         bd.mean=bound.means,
#                                         nits=n.iters,
#                                         budget.prop.list=win.y.budget.props2,
#                                         alph=alpha,...
#
# ){
#   print(idx)
#   temp.out=dp_synthdata(formula=formulas[idx],
#                         confidential.data=confidential.data,
#                         synth.data=synth.data,
#                         synth.epsilon=NULL,
#                         synth.delta=NULL,
#                         continuous.vars=NULL,
#                         num.bin=NULL,
#                         bin.param=NA,
#                         assign.type="simple",
#                         blocks=NULL,
#                         clusters=NULL,
#                         within.blocks=TRUE,
#                         epsilon.list=as.list(per.response.eps*budget.prop.list),
#                         delta.list=as.list(per.response.del*budget.prop.list),
#                         bd.sd.list=bd.sd,
#                         bd.mean.list=bd.mean,
#                         alphas.list=alph,
#                         num.iters=nits,
#                         treatment.var=tr.vars,
#                         rseed=NA,
#                         return.time=TRUE,
#                         return.confidential.table=FALSE,
#                         return.san.summary=FALSE,
#                         use.san.residerror=FALSE,
#                         return.treatment.pvalue=FALSE,
#                         return.treatment.CI=FALSE,
#                         treat.ppart.list=NULL,
#                         ...)
#
#   response.var=y.vars[idx]
#   temp.synth=temp.out[[1]]
#   list(temp.synth[,response.var],temp.out[[2]])
# }
# get_synth_response_synthdatadp_noresiderr(1)
#
# synth.response=lapply(seq(1,n.response),get_synth_response_synthdatadp_noresiderr)
# full.dp.synth2=dplyr::bind_cols(mv.covariates3,dplyr::bind_cols(sapply(synth.response,"[[",1)))
# full.dp.synth2$control=ifelse(full.dp.synth2$treatment=="control",1,0)
# full.dp.comp.times2=list("total.dp.time"=(proc.time()-synthdata.start2)[[3]],
#                         "dp.response.time"=sapply(synth.response,"[[",2),
#                         "mv.T.X.time"=mv.covariate.comp.time2)

##############
comp.times=list("MV Histogram"=mv.synth.comp.time,
                "Hybrid"=hybrid.dp.comp.times,
                "Full Model-Based"=full.dp.comp.times)#,
#"Full Model-Based Not Using ResidError"=full.dp.comp.times2)


print("Starting Table")

table.start=proc.time()
sterr.type="HC1"
stderr.func=function(mod,std.type=sterr.type){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}
mod.names=response.vars.all
col.order=c("ModelName","ControlMean",
            paste0(c("ITT_","StdErr_","AdjPvalue_"),treat.vars[1]),
            paste0(c("ITT_","StdErr_","AdjPvalue_"),treat.vars[2]),
            paste0(c("ITT_","StdErr_","AdjPvalue_"),treat.vars[3]),"nObs")

#when this code is run in a Rmarkdown with yaml header containing:
#output: pdf_document:
#  keep_tex: yes
#in a chunk {r,results="asis"} it formats the table properly
out.table=print_compare_ITTtable(list("Confidential Data"=sim.data,
                                      "DP Multivariate Histogram Synthetic Data"=mv.out,
                                      "Not-DP Hybrid Synthetic Data"=hybrid.synthdata,
                                      "DP Model-based Synthetic Data Using Sanitized Residual Error"=full.dp.synth),
                                 #"DP Model-based Synthetic Data Using Sanitized Coefficient Standard Error"=full.dp.synth2),
                                 reg.formulas,
                                 response.vars=response.vars.all,
                                 families=families,
                                 treat.vars=treat.vars,
                                 control.var="control",
                                 bonferroni.npvals=NULL,
                                 mult.test.correct=c("treatments","responses"),
                                 add.pval.stars=TRUE,
                                 model.names=mod.names,
                                 stderr.func=stderr.func,
                                 digits=3,
                                 label=NULL,
                                 caption=NULL,
                                 col.var.names = col.order,
                                 treat.names=c("Therapy Only","Cash Only","Both"),
                                 booktabs=TRUE,
                                 align="lcclcclcclcc")

out.table=kableExtra::kable_styling(out.table,bootstrap_options = "striped")
out.table=kableExtra::column_spec(out.table,1,width="1000em")
for(i in seq(2,2+(length(treat.vars)*3))){
  out.table=kableExtra::column_spec(out.table,i,width="20em")
}
out.table=kableExtra::add_footnote(out.table,"Notes: The table reports intent to treat estimates of each treatment are after 12-13 weeks. The standard errors reports are heteroskedastic robust standard errors with the unajusted p-values reported as follows: ***p<0.01, **p<0.05$, *p<0.1. The adjusted p-values are adjusted for the 8 outcomes and 3 treatment arms using a Bonferroni adjustment. Blattman et al. use a Westfall-Young adjustment method.")
kableExtra::save_kable(out.table,paste0(basepath,"/inst/tables/sim_ITT_table_compare.png"),keep_tex=F)

end.time=proc.time()
comp.times=c(list("total.time"=(end.time-main.start)[[3]],
                  "table.time"=(end.time-table.start)[[3]]),
             comp.times)
comp.times



