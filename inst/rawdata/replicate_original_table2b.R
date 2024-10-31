####### Processing  Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

########## Global for Model Variables ############
## To be used to subset replication data

#Replicating Table 2b in Original Paper using AER.do file as reference

#base variable in AER.do
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


#dvs_t2 in AER.do
response.stem=c('fam_asb_lt','drugssellever','stealnb','disputes_all_z','carryweapon',
                'arrested','asbhostilstd','domabuse_z')

#strata in AER.do
block.vars=c("tp_strata_alt","cg_strata") #block variables

#focusing on round==5 (term=lt in AER.do file)
#read in the liberia data
liberia.full=haven::read_stata(paste0(basepath,"/inst/rawdata/STYL_Final.dta"))
#Data from:
##Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. Replication
##data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral
##Therapy in Liberia. Nashville, TN: American Economic Association [publisher], 2017.
##Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor],
##2019-10-12. https://doi.org/10.3886/E113056V1
liberia.full=labelled::remove_attributes(liberia.full,"format.stata")


y.sub.data=liberia.full[,grepl(paste0(response.stem,collapse="|"),colnames(liberia.full))]
# Instructions.pdf say "_st" and "_b" are short term and pretreatment values
sapply(y.sub.data[,grepl("_lt",colnames(y.sub.data))],function(x)attr(x,"label"))
#fam_asb_lt is the correct choice. The other repsonse variables have "_ltav"

#list all possible response variables for round==5
response.vars.allltav = c('fam_asb_lt',paste0(response.stem[-1],"_ltav"))
treat.vars=c("cashassonly","tpassonly","tpcashass") #treatment variables
########

### If not investigating the choices of what variables are continuous and which
### are categorical, you can start with "LiberiaRound5.Rda" instead of "STYL_Final.dta"
### However, the following code shows how LiberiaRound5.Rda is created by subsetting
### and processing STYL_Final.dta.

col.vars=c(response.vars.allltav,
           treat.vars,block.vars,baseline.full,c("control"))

#line 334 in AER.do shows the regression model to restruct to unfound_wave==0 and round==5 (for long term)
row.conditions=((liberia.full$round==5)&(liberia.full$unfound_wave==0))
#subset data
liberia.sub=liberia.full[row.conditions,col.vars[!duplicated(col.vars)]]

## Optional: save liberia.sub
#save(liberia.sub,file=paste0(basepath,"/data/LiberiaRound5.Rda")) #save smaller data

## If starting with "LiberiaRound5.Rda" uncomment the following line:
#load(paste(basepath,"LiberiaRound5.Rda",sep="/"))

covariate.formula=paste(c(baseline.full,block.vars),collapse="+")

#liberia.sub[,colnames(liberia.sub)]=
#  apply(liberia.sub[,colnames(liberia.sub)],2,
#        function(x)as.numeric(as.character(x)))

#line 33 of AER.do sets block variables to factors
liberia.sub[,block.vars]=
  apply(liberia.sub[,block.vars],2,
        function(x)as.factor(as.character(x)))
#on line 334 of AER.do file the ",r" indicates that robust (sandwich estimator) of variance
reg.formula=paste(response.vars.allltav,
                  paste(c(treat.vars,baseline.full,block.vars),collapse="+"),sep="~")
sterr.func=function(mod){
  sqrt(diag(sandwich::vcovHC(mod)))
}
ITTtable(data=liberia.sub,
         families=rep("gaussian",length(reg.formula)),
         reg.models=reg.formula,
         treat.vars=treat.vars,
         control.var="control",
         stderr.func=sterr.func,
         add.pval.stars = FALSE,
         mult.test.correct=NULL)
#The robust sandwhich method does not work. However there are other types of it.

sterr.type="HC0"
sterr.func=function(mod,std.type=sterr.type){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}


type.list=c("HC0","HC1","HC2","HC3","HC4","HC4m","HC5","const")
for(sterr.type in type.list){
    tab=ITTtable(data=liberia.sub,
                   families=rep("gaussian",length(reg.formula)),
                 reg.models=reg.formula,
                 treat.vars=treat.vars,
                 control.var = "control",
                 stderr.func = sterr.func,
                 add.pval.stars = FALSE,
                 mult.test.correct=NULL)
    if((sum(is.nan(tab$StdErr_cashassonly))==0)&&
      (sum(abs(round(tab$StdErr_cashassonly,3)-
             c(0.097,0.030,0.388,0.090,0.035,0.025,0.107,0.113)))<10^(-6))){
      print(sterr.type)
        tempdf=tab[,!grepl("Pvalue",colnames(tab))]
        print(tempdf[,colnames(tempdf)!="nObs"])
      }

    }
#### Anaylsis is done with additive block indicators, _ltav response variables,
## and robust HC type="HC1" errors

### Function to Check Regression Assumptions ###
reg_assumptions=function(reg.model,data=liberia.sub,family="gaussian",response.var=NULL,pt.sz=2.5,ln.sz=1.5,bs.sz=7,...){
  mod.fit=glm(reg.model,data=data,family=family,... )
  if(is.null(response.var)==TRUE){
    response.var=names(attr(tempmod$terms,"dataClasses"))[1]
  }
  data.temp=data.frame("fitted"=mod.fit$fitted,"residuals"=mod.fit$residuals,"zeros"=0)
  plot.fit.resid=ggplot2::ggplot(data.temp,ggplot2::aes(x=fitted,y=residuals))+
    ggplot2::geom_point(size=pt.sz)+
    ggplot2::geom_smooth(se=F,col="red",linewidth=ln.sz)+
    ggplot2::geom_line(ggplot2::aes(x=fitted,y=zeros),col="black",linetype=2,size=ln.sz)+
    ggplot2::labs(x=paste("Fitted",response.var),y="Residuals")+
    ggplot2::ggtitle(paste0(response.var,": Residual vs. Fitted"))+
    ggplot2::theme_minimal(base_size=bs.sz)
  data.temp$norm=stats::rnorm(nrow(data.temp))
  plot.qqnorm=ggplot2::ggplot(data=data.temp,ggplot2::aes(sample=residuals))+
    ggplot2::geom_qq(size=pt.sz,color="black")+
    ggplot2::geom_qq_line(color="red",linewidth=ln.sz)+
    ggplot2::labs(x="Theoretical Quantiles",y="Residuals Quantiles")+
    ggplot2::ggtitle(paste0(response.var,": Normal Q-Q Plot"))+
    ggplot2::theme_minimal(base_size=bs.sz)



  #variance.test=car::ncvTest.glm(mod.fit)
  return(cowplot::plot_grid(plot.fit.resid,plot.qqnorm,ncol=2))
}
tempout= lapply(seq(1,length(response.vars.allltav)),
                function(idx)reg_assumptions(reg.formula[idx],
                                             response.var=response.vars.allltav[idx],
                                             pt.sz=3,
                                             ln.sz=2,
                                             bs.sz=20))
png(paste0(basepath,"/inst/rawdata/diagnosticFitPlots.png"),width=3000,height=1600)
print(cowplot::plot_grid(plotlist=tempout,ncol=2))
dev.off()

