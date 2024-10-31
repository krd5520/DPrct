################ Full DP Synthetic epsilon=1 #######################
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
win.y.x.coefs=1-win.y.resvar-win.y.trcoef
############

### other params #####
bound.means=100
bound.sds=c(2^(-15),2^(15))
n.iters=300
########

##############################
n.x.coefs=length(baseline.full)+1
n.treat.coefs=length(treat.vars)
win.y.budget.props=c(win.y.resvar,rep(win.y.trcoef,n.treat.coefs),
                     rep(win.y.x.coefs,n.x.coefs))
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
                                        keep.block.vars=TRUE,
                                        block.vars=block.vars,
                                        block.sizes=NULL,
                                        conditions=c(treat.vars,"control"))
for(tr.var in c(treat.vars,"control")){
  mv.covariates2[,tr.var]=as.factor(ifelse(mv.covariates2$treatment==tr.var,1,0))
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
                                        budget.prop.list=win.y.budget.props

){
  temp.out=dp_synthdata(formula=formulas[idx],
                        confidential.data=liberia.sub,
                        synth.data=NULL,
                        synth.epsilon=NA,
                        synth.delta=0,
                        continuous.vars=NULL,
                        num.bin=NULL,
                        bin.param=NA,
                        assign.type="simple",
                        blocks=NULL,
                        clusters=NULL,
                        within.blocks=TRUE,
                        epsilon.list,
                        delta.list=0,
                        bd.sd.list=c(2^(-15),2^(15)),
                        bd.mean.list,
                        alphas.list=0.05,
                        num.iters=1e4,
                        treatment.var=NULL,
                        rseed=NA,
                        return.time=FALSE,
                        return.confidential.table=FALSE,
                        return.san.summary=TRUE,
                        use.san.residerror=TRUE,
                        return.treatment.pvalue=TRUE,
                        return.treatment.CI=TRUE,
                        treat.ppart.list=NULL,
                        ...){

    response.var=y.vars[idx]
    temp.synth=temp.out[[1]]
    list(temp.synth[,response.var],temp.out[[2]])
  }










  ##############
