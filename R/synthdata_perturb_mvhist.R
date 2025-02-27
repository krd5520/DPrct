#' Make Approximate DP synthetic data from mixed data types
#'
#' @description The synthetic data is created using perturbed multivariate
#'      histogram methods  that satisfy (\code{epsilon}, \code{delta})-differential privacy. The
#'      methods are introduced by Dwork et al. (2006) and Wasserman and Zhou (2008).
#'
#' @param data a data.frame (or vector of data) to generate (\code{epsilon}, \code{delta})-DP version of
#' @param epsilon a positive value that defines the privacy cost
#' @param delta a nonnegative value that defines the privacy cost. Default is \code{delta=0}.
#' @param continuous.vars vector of strings, indices, or logical values to indicate
#'    the columns of \code{data} with continuous variables. If \code{continuous.vars=NULL},
#'    all variables are assumed to be categorical.
#' @param num.bin a numeric value or vector of numeric values for the number of
#'    bins for each continuous variable. If \code{num.bin=NA}, the \code{bin.param}
#'    is used to determine the number of bins.
#' @param bin.param a numeric value in (0,1) to determine number of bins for
#'    continuous variables to be used if \code{num.bins=NA}.
#' @param add.cont.variation a logical if \code{TRUE}, add uniform noise to midpoint
#'    to get a random value within the range of bin. Default is \code{FALSE}.
#' @param with.treatment a logical if \code{TRUE}, randomly assignment treatments
#'    are added using \code{treatment_assign} function. Default is \code{TRUE}.
#' @param assign.type a string taking values "simple","complete","block",
#'    "cluster", or "block_and_cluster" to specify what kind of random assignment
#'    should be used. Default is "simple".
#' @param rseed a value to set the randomization seed. If \code{rseed=NA} (default),
#'    then no random seed it set.
#' @param return.time a logical to indicate if the computation times should be
#'    returned with the synthetic data. Default is \code{FALSE}.
#' @param treatment.colname a string to rename the treatment column. Default is "treatment".
#' @param within.blocks a logical to indicate if the data should be sampled within the
#'      block levels in order to maintain the same block sizes as the original data.
#' @param blocks column name(s) for the variables for block treatment assignment. This
#'      is only needed if \code{assign.type=="block"} or \code{=="block_and_cluster"}.
#' @param block.sizes sample size of each block level. If NULL, then the size of the
#'      blocks in the \code{data} will be used.
#' @param perturb is an indicator as to whether the histogram should be perturbed to
#'      satisfy DP. Otherwise it will just be a resampled dataset.
#' @return a synthetic data.frame for confidential data (inputted as \code{data})
#'    that satisfies (\code{epsilon},\code{delta})-DP using the perturbed
#'    multivariate histogram method.
#' @param ... additional arguments are such as \code{blocks} or \code{clusters}
#'    are supplied to the \code{treatment_assign} function.
#' @importFrom stats setNames
#'
#'
#' @examples
#' synthdata_perturb_mvhist(mtcars[,2:6],1,continuous.vars=c(FALSE,TRUE,TRUE,TRUE,TRUE),num.bin=3)
#'
#' @references
#' Dwork, C., McSherry, F., Nissim, K., and Smith, A. (2006). Calibrating noise to
#'     sensitivity in private data analysis. volume Vol. 3876, pages 265–284.
#'Wasserman, L. and Zhou, S. (2008). A statistical framework for differential privacy.
#'
#'
#' @family syntheticData
#' @export
#'

synthdata_perturb_mvhist<-function(data,
                                   epsilon,
                                   delta=0,
                                   continuous.vars=NULL,
                                   num.bin=NULL,
                                   bin.param=NA,
                                   add.cont.variation=FALSE,
                                    treatment.colname="treatment",
                                   with.treatment=TRUE,
                                   assign.type="simple",
                                   return.time=FALSE,
                                   rseed=NA,
                                   within.blocks=TRUE,
                                   blocks=NULL,
                                   clusters=NULL,
                                   block.sizes=NULL,
                                   factorial=FALSE,
                                   conditions=NULL,
                                   perturb=T,
                                   ...){

  start.time=proc.time()
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }

  #if there is only 1 column to make synthetic (aka data is a vector)
  if((is.data.frame(data)==FALSE)&&(length(data)>0)){
    data=data.frame("covariate"=data)
  }
  ### Check Inputs ###
  stopifnot(is.numeric(epsilon)&&epsilon>0)

  which.cat=sapply(seq(1,ncol(data)),function(x)is.factor(data[,x])|is.character(data[,x]))
  #get multivariate histogram for mixed data types.
  mv.hist.out=multivariate_histogram(data=data,continuous.vars = continuous.vars,
                                     num.bin = num.bin,bin.param=bin.param,
                                     which.cont.out=TRUE,check.cont=T)
  if(is.list(mv.hist.out)==TRUE){
    freq.df<-mv.hist.out[[1]]
    which.cont=mv.hist.out[[2]]
  }else{
    freq.df=mv.hist.out
    which.cont= rep(FALSE,ncol(data)) #else all variables categorical
  }
  mv.hist.stop=proc.time()
  mv.hist.time=(mv.hist.stop-start.time)[[3]]
  if(perturb==TRUE){
    freq.df=dp_perturbed_hist(hist.df=freq.df,epsilon = epsilon,delta=delta)
  }else{
    freq.df$san.prop=freq.df$Freq/sum(freq.df$Freq)
  }

  san.hist.stop=proc.time()
  san.hist.time=(san.hist.stop-mv.hist.stop)[[3]]

  #if treatment within blocks with set block size, must maintain that block size
  if((within.blocks==TRUE)&(is.null(blocks)==FALSE)){
    freq.df$san.prop.in.block=NA #initialize a new prop to be normalized to block size
    row.sample=NULL #initialize row.sample vector of indices to be sampled from freq.df
    full.idx=seq(1,nrow(freq.df)) #index rows of freq.df
    if(length(blocks)==1){ #if only one block variable
      #get sum of sanitized probabilities per block level
      within.blocks.san.prop=tapply(freq.df$san.prop,freq.df[,blocks],sum)
      if(is.null(blocks)==TRUE){ #if no block size is given, use block size from data
        block.sizes=table(as.factor(freq.df[,blocks]))
      }else if(is.null(names(block.sizes))==TRUE){ #if no level names given for
        if(length(block.sizes)==1){ #if one block size given, assume it is constant across levels
          block.sizes=rep(block.sizes,length(names(within.blocks.san.prop)))
        }
        names(block.sizes)=names(within.blocks.san.prop) #give names of levels
      }
      for(block.level in names(within.blocks.san.prop)){ #for each block level
        block.level.indic=(freq.df[,blocks]==block.level) #logical for block level
        #renormalize for block level
        freq.df$san.prop.in.block[block.level.indic]=
          freq.df$san.prop[block.level.indic]/within.blocks.san.prop[block.level]
        #sample block.size of the bin with this block.level with replacement
        block.level.sample=sample(sum(block.level.indic),block.sizes[block.level],
                                  replace = TRUE,
                                  prob = freq.df$san.prop.in.block[block.level.indic])
        sub.idx=full.idx[block.level.indic] #get the indices corresponding to this level
        row.sample=c(row.sample,sub.idx[block.level.sample]) #get the indices for this sample
      }
    }else{ #more than 1 block variable
        within.blocks.san.prop=sapply(blocks,
                      function(bl.var)tapply(freq.df$san.prop,
                                             freq.df[,bl.var],sum))
      if(is.null(block.sizes)==TRUE){
        block.sizes=lapply(blocks,
                           function(bl.var)table(freq.df[,bl.var]))
      }else if(length(block.sizes)==1){
        block.sizes=lapply(blocks,
                           function(bl.v)stats::setNames(rep(block.sizes,nlevels(as.factor(freq.df[,bl.v])),
                                                  levels(as.factor(freq.df[,bl.v])))))
      }
      for(i in seq(1,length(blocks))){
        bl.sizes=block.sizes[[i]]
        within.levels.san.prop=within.blocks.san.prop[[i]]
        if(is.null(names(bl.sizes))==TRUE){
        if(length(bl.sizes)==1){
          bl.sizes=rep(bl.sizes,length(levels(as.factor(freq.df[,blocks[i]]))))
        }
          names(bl.sizes)=names(within.levels.san.prop)
        }
        for(l in names(within.levels.san.prop)){
          block.level.indic=(freq.df[,blocks[i]]==l) #logical for block level
          #renormalize for block level
          freq.df$san.prop.in.block[block.level.indic]=
            freq.df$san.prop[block.level.indic]/within.levels.san.prop[l]
          #sample block.size of the bin with this block.level with replacement
          block.level.sample=sample(sum(block.level.indic),bl.sizes[l],
                                    replace = TRUE,
                                    prob = freq.df$san.prop.in.block[block.level.indic])
          sub.idx=full.idx[block.level.indic] #get the indices corresponding to this level
          row.sample=c(row.sample,sub.idx[block.level.sample]) #get the indices for this sample
        }

      }
    }

    # #if using delta>0, and number of bins is high enough
    # if((delta>0)&(num.bins.level>(2/delta))){ #use Bun et al. 2016
    #   threshold=((2*base::log(2/delta))/(epsilon/n))+(1/n)
    # }else{ #pure-DP and low number of bins don't use a threshold
    #   threshold=0
    # }

  }
  #sample rows of mv hist with probabilities equal to norm.san with replacement
  # get sample of size equal to number of rows of data.
  row.sample<-sample(1:nrow(freq.df),nrow(data),
                           replace=TRUE,prob=freq.df$san.prop)

  #synthetic data has values from the sample of the histogram
  # (remove frequency, san.prop columns)
  synth.data<-freq.df[row.sample,colnames(freq.df)%in%colnames(data)]

  #force columns to be factors or numeric (add variation if needed)
  synth.data[,which.cat]=lapply(synth.data[,which.cat],
                                function(col)as.factor(as.character(col)))
  if(add.cont.variation==TRUE){# if adding uniform variation for continuous values
    cont.v=colnames(synth.data[,which.cont,drop=F])
    if(length(cont.v)==0){
      synth.data=synth.data
    }else if(length(cont.v)==1){
      synth.data[,cont.v]=synth_continuous_variation(synth.data[,cont.v])
    }else{
    count.nas=sapply(cont.v,function(x)sum(is.na(synth.data[,cont.v])))
    count.vals=sapply(cont.v,function(x)length(unique(as.character(synth.data[,cont.v])))<=1)
    if(sum(count.nas,na.rm=T)>0){
      warning(paste(cont.v[count.nas>0],"has",count.nas," NA values",collapse="\n"))
    }
    if((sum(count.vals,na.rm=T)>0)|(sum(is.na(count.vals))>0)){
      warning(paste("some continuous variable only has one value?:",paste0(cont.v[count.vals],collapse=", "),
                    " or something weird with NAs ",sum(is.na(count.vals))))
    }
    synth.data[,which.cont]=lapply(synth.data[,which.cont],synth_continuous_variation)
    }
  }else{
      synth.data[,which.cont]=lapply(synth.data[,which.cont],
                                           function(col)as.numeric(as.character(col)))
  }
  synth.data[,!which.cat]=lapply(synth.data[,!which.cat],
                                function(col)as.numeric(as.character(col)))



  samp.hist.stop=proc.time()
  samp.hist.time=(samp.hist.stop-san.hist.stop)[[3]]
  if(with.treatment==TRUE){
    if(factorial==TRUE){
      #warning("in factorial==TRUE")
      if(is.null(conditions)==TRUE){
        if(length(treatment.colname)>1){
          #warning(paste("treatment.colname is a vector.",paste0(treatment.colname,collapse=", ")))
          conditions=treatment.colname
        }else{
          conditions=c("T1","T2","Both")
        }
      }
      conditions=conditions[conditions!="control"]
      nblocks=length(blocks)
      if((nblocks>2)|(length(conditions)>3)){
        stop("Only 2X2 factorial design supported")
      }
        if(nblocks==0){
          blocks.ls=list(NULL,NULL)
        }else if(nblocks==1){
          blocks.ls=list(blocks,blocks)
        }else if(nblocks==2){
         blocks.ls=list(blocks[1],blocks[2])
        }
      if(length(treatment.colname)!=1){
        warn.mess=paste("column named 'treatment' will summarize the treatment variables.",
                        ifelse(length(treatment.colname)==0," ",
                               "Inputted treatment.colname will be the column names of the indicator columns."))
        warning(warn.mess)
        treatment.colname="treatment"
      }
      for(i in c(1,2)){
        #warning(paste("inside for loop i=",i))
#        if((length(treatment.colname)==length(conditions))&(length(treatment.colname)>1)){
#          warning("treatment.colname has same length as conditions")
          synth.data=treatment_assign(synth.data=synth.data,
                                      assign.type=assign.type,
                                      treatment.colname=conditions[i],
                                      blocks=unlist(blocks.ls[[i]]))#,conditions=c("1","0"))#,#clusters=clusters, ...)
          # warning(paste("done with treatment_assign.",
          # "dim of treatment is",paste0(dim(synth.data[,ncol(synth.data),drop=F]),collapse=", "),
          # "synthhead of treat col",paste0(head(c(synth.data[,ncol(synth.data),drop=T])),collapse=", ")))
          tr.col=as.numeric(c(unlist(synth.data[,ncol(synth.data)]))) #newly added treatment column

          synth.data[,ncol(synth.data)]=tr.col
          #colnames(synth.data)=c(colnames(synth.data)[-1],conditions[i])
          }
      treff.rowsums=rowSums(synth.data[,colnames(synth.data)%in%conditions])
      synth.data$control=ifelse(treff.rowsums==0,1,0)
      synth.data[,conditions[3]]=ifelse(treff.rowsums==2,1,0)
      # }
      #
     #  if(length(treatment.colname)==length(conditions)){
     #    synth.data$control=ifelse(rowSums(synth.data[,colnames(synth.data)%in%treatment.colname])==0,1,0)
     #
     #  }
     #  treateffs=synth.data[,colnames(synth.data)%in%conditions]
     #  treatcol=sapply(seq(1,nrow(treateffs)),function(x)paste0(treateffs[x,],collapse=""))
     #  #warning(paste("length of treat col is",length(treatcol),"length of data is",nrow(synth.data)))
     #  synth.data[,treatment.colname]=treatcol
     #  warning(paste("unique treatcol is",paste0(unique(unlist(synth.data[,treatment.colname]))),collapse=", "))
     #  synth.data[,treatment.colname]=base::trimws(synth.data$treatment)
     #  warning(paste("control count is",sum(synth.data[,treatment.colname]==" ")))
     #  synth.data[synth.data[,treatment.colname]==" ",treatment.colname]="control"
     #
     #  cond.idx=length(blocks)+1
     #  #warning("before i,j")
     #  for(i in seq(1,length(blocks))){
     #    for(j in seq(i+1, length(blocks))){
     #      synth.data[synth.data$treatment==paste0(conditions[i],conditions[j]),treatment.colname]=conditions[cond.idx]
     #      cond.idx=cond.idx+1
     #    }
     #  }
     # # warning("outside i,j ")
     #  for(cond in c(conditions,"control")){
     #    synth.data[,cond]=ifelse(synth.data$treatment==cond,1,0)
     #  }
    }else{
      #warning("factorial==FALSE")
    synth.data=treatment_assign(synth.data=synth.data,
                                 assign.type=assign.type,
                                 treatment.colname=treatment.colname,
                                 blocks=blocks,conditions=conditions,#clusters=clusters,
                                 ...)
    #### NOTE: double check this handles multiple treatment variables as well?
    }
    #warning("end within.treatment==TRUE")
  }

  #warning(paste0("control group has ",sum(synth.data$control==1)))

  #warning(paste("colnames are ",paste0(colnames(synth.data),collapse=", ")))
  attr(synth.data,"priv.cost")=c("epsilon"=epsilon,"delta"=delta)

  if(return.time==TRUE){
    #warning("inside return.time==TRUE")
    end.time=proc.time()
    attr(synth.data,"comp.time")=c("mv.hist.time"=mv.hist.time,"san.hist.time"=san.hist.time,
                                   "samp.hist.time"=samp.hist.time,
                                   "total.synthdata_perturb_mvhist.time"=(end.time-start.time)[[3]])
  }
  return(synth.data)
}

#' Internal function: add uniform variation to continuous values within the interval of thier category
#'
#' @param cont.var is vector of the discretized interval midpoints of a continuous random variable as a factor
#' @return a numeric vector of synthetic continuous values
#' @noRd
#' @keywords internal
#'
#' @importFrom stats runif
synth_continuous_variation<-function(cat.var){
  n.rw=length(cat.var)
  midpoints<-base::sort(as.numeric(as.character(cat.var))) #get midpoint values
  if(sum(is.na(midpoints))>0){
    warning(paste("There are NA values in midpoints:",sum(is.na(midpoints)),"NA values"))
  }
  n.levels=length(midpoints)
  #half the interval length is 1/2 difference between midpoints
  half.widths=abs(midpoints[2:n.levels]-midpoints[1:n.levels-1])/2

  if(sum(!is.na(half.widths))==0){
    stop(paste("all half widths are NA. cat.var head is:",paste(head(cat.var),sep=", "),
               "dim is:",paste(try(dim(cat.var)),collapse=", ")))
  }
  if((sum(half.widths,na.rm=T)/sum(!is.na(half.widths)))==half.widths[1]){ #if all widths are equal
    variation<-stats::runif(n.rw,-half.widths[1],half.widths[1]) #uniform rv
    cont.var=as.numeric(as.character(cat.var))+variation
  }else{ #if all widths are NOT equal
    cont.var=rep(NA,n.rw) #initialize column
    for(hw in unique(half.widths)){ #for each unique width value
      level.idx=which(half.widths==hw) #find corresponding levels
      cat.idx=which(cat.var %in% unique(as.character(cat.var))[level.idx]) #find rows with those levels
      variation=stats::runif(length(cat.idx),-hw,hw)
      cont.var[cat.idx]=as.numeric(as.character(cat.var[cat.idx]))+variation
    }
  }
    return(as.numeric(cont.var))
}

