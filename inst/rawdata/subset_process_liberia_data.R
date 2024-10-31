####### Processing  Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

########## Global for Model Variables ############
## To be used to subset replication data

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
########

### If not investigating the choices of what variables are continuous and which
### are categorical, you can start with "LiberiaRound5.Rda" instead of "STYL_Final.dta"
### However, the following code shows how LiberiaRound5.Rda is created by subsetting
### and processing STYL_Final.dta.

#read in the liberia data
## If starting with "LiberiaRound5.Rda" comment the following line.
liberia.full=haven::read_stata(paste0(basepath,"/inst/rawdata/STYL_Final.dta"))
#Data from:
##Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. Replication
##data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral
##Therapy in Liberia. Nashville, TN: American Economic Association [publisher], 2017.
##Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor],
##2019-10-12. https://doi.org/10.3886/E113056V1

######## Subset and Process Replication Data #######
## If starting with "LiberiaRound5.Rda" comment the following 4 lines:
liberia.full=labelled::remove_attributes(liberia.full,"format.stata")
col.vars=c(response.vars.all,treat.vars,block.vars,baseline.full,c("control"))
row.conditions=((liberia.full$round==5)&(liberia.full$unfound_wave==0))
liberia.sub=liberia.full[row.conditions,col.vars[!duplicated(col.vars)]]

## Optional: save liberia.sub
#save(liberia.sub,file=paste0(basepath,"/LiberiaRound5.Rda")) #save smaller data

## If starting with "LiberiaRound5.Rda" uncomment the following line:
#load(paste(basepath,"LiberiaRound5.Rda",sep="/"))


#### Further processing of liberia.sub data

#combine treatment variables into 1 column
# (this is useful for the dp_synthdata and hybrid_dp function)
liberia.sub$treatment=ifelse(liberia.sub$cashassonly==1,"cashassonly","tpassonly")
liberia.sub$treatment[liberia.sub$tpcashass==1]="tpcashass"
liberia.sub$treatment[liberia.sub$control==1]="control"

# Investigate what are continuous and what are categorical variables.
# These decisions where checked manually using the labels of liberia.sub
cov.tab.vals=sapply(colnames(liberia.sub),
                    function(rvar)length(table(liberia.sub[,rvar],useNA="ifany")))

#<=3 is often 0,1 indicators with potential NA rows.
#Also includes "literacy_b" variable which is a literacy index taking value 0,1,2
cat.vars=names(cov.tab.vals)[cov.tab.vals<=3]
## Uncomment following line to investigate the labels for these variables
#str(liberia.sub[,cat.vars])


low.value.counts=sapply(names(cov.tab.vals)[cov.tab.vals<10&cov.tab.vals>3],
                        function(x)table(liberia.sub[,x]))
low.value.counts
#something strange is happening with mpartners_b (1 person has 0.10000005?)
t(liberia.sub[round(liberia.sub$mpartners_b,3)==0.1,])
liberia.sub.wodd=liberia.sub #save dataframe with the strange entry for later comparison
liberia.sub=liberia.sub[!(round(liberia.sub$mpartners_b,3)==0.1),] #remove odd entry

#looking at the labels for the low value entries
sapply(liberia.sub[,names(low.value.counts)],function(x)attr(x,"label"))
### Most of these are discrete counts, or indexes. Indexes could be categorical or numeric.
## When treating them as numeric we assume the change in response when there is an increase of 1 unit
## 1 unit in the index score (holding all else constant) is constant no matter what the value
## of the base index score was.
## The original replication material appear to treat these variables as numeric/continuous, so
## we will also consider them numeric/continuous when fitting our models.


attr(liberia.sub$famseeoften_b,"labels")
##we see "famseeoften_b" is categorical
##   with values 1,2,3,4 indicating "everytime", "somethime", "one one time" and "never"
##   and values 97 and 98 indicate "don't know" and "refuse to answer".

cat.vars=c(cat.vars,block.vars,"famseeoften_b") #block variables are also categorical

## for the purpose of the multivariate histograms, variables taking a limited number of values
## can be treated as categorical. However, they will be reverted to continuous for model fitting.
cont.as.cat=names(cov.tab.vals)[(cov.tab.vals<80)&(!(names(cov.tab.vals)%in%cat.vars))]
mid.value.counts=sapply(cont.as.cat[!(cont.as.cat %in% names(low.value.counts))],
                        function(x)table(liberia.sub[,x]))
mid.value.counts
mid.labs=sapply(liberia.sub[,names(mid.value.counts)],function(x)attr(x,"label"))
mid.labs[!(grepl("index",mid.labs,ignore.case = T))]
#These labels are years, hours (or average hours), and discrete counts so are continuous/numeric


cont.as.cont=c(names(cov.tab.vals)[cov.tab.vals>=80])
cov.tab.vals[names(cov.tab.vals)%in% cont.as.cont]
high.labs=sapply(liberia.sub[,cont.as.cont],function(x)attr(x,"label"))
high.labs #all appear to be continuous numeric

#continuous variables
cont.vars.all=c(cont.as.cat,cont.as.cont)
#continuous variables to treat as continuous
liberia.sub[,cat.vars]=apply(liberia.sub[,cat.vars],2,as.factor)

liberia.sub=labelled::remove_labels(liberia.sub)

clean.data.time=(proc.time()-main.start)[[3]]
attr(liberia.sub,"column.classifications")=
  list("cat.vars"=cat.vars,
       "cont.as.cont"=cont.as.cont,
       "cont.as.cat"=cont.as.cat)
attr(liberia.sub,"processing.time")=clean.data.time

file.location=paste0(basepath,"/data/SubsetLiberiaWithFactors.Rda")

save(liberia.sub,file=file.location)
