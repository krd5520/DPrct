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
#' @param within.blocks (antiquated)
#' @param blocks column name(s) for the variables for block treatment assignment. This
#'      is only needed if \code{assign.type=="block"} or \code{=="block_and_cluster"}.
#' @param block.sizes (antiquated)
#' @param continuous.limits a list of upper and lower bounds for each continuous variable.
#'    If there are continuous variables, then a limit must be supplied for each.
#' @param perturb is an indicator as to whether the histogram should be perturbed to
#'      satisfy DP. Otherwise it will just be a resampled dataset.
#' @param standardized.cont a vector of column names to transform to a standardized
#'    form $(x_i-\bar{x})/std.dev.(x)$. If \code{NULL}, then no variable is transformed.
#' @param std.limits a numeric value for the number of standard deviations away from the mean
#'    to set the transformed standardized column limits to. Default is 15.
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
                                   within.blocks=FALSE,
                                   blocks=NULL,
                                   clusters=NULL,
                                   block.sizes=NULL,
                                   factorial=FALSE,
                                   conditions=NULL,
                                   perturb=T,
                                   standardize.cont=NULL,
                                   std.limits=15,
                                   continuous.limits=NULL,
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
                                     continuous.limits=continuous.limits,
                                     num.bin = num.bin,bin.param=bin.param,
                                     which.cont.out=TRUE,levels.out=T,std.limits=std.limits,standardize.cont=standardize.cont)
  #if(is.list(mv.hist.out)==TRUE){
  freq.df<-mv.hist.out[["mv.histogram"]]
  which.cont=mv.hist.out[["which.cont"]]
  levels.list=(mv.hist.out[["levels.list"]])[[1]]
  #}else{
  #  freq.df=mv.hist.out
  #  which.cont= rep(FALSE,ncol(data)) #else all variables categorical
  #}
  mv.hist.stop=proc.time()
  mv.hist.time=(mv.hist.stop-start.time)[[3]]

  san.prop.zero.to.add=NULL
  if(perturb==TRUE){
    possible.combos=prod(sapply(levels.list,length))
    freq.df=dp_perturbed_hist(hist.df=freq.df,epsilon = epsilon,delta=delta,possible.combos=possible.combos)
    if(is.list(freq.df)==TRUE){
      san.prop.zero.to.add=freq.df[[2]]
      freq.df=freq.df[[1]]
    }
  }else{
    warning("No privacy noise added since perturb=FALSE was inputted. This is to produce a epsilon=infinity case that is simply a resampling algorithm, not a DP algorithm.")
    freq.df$san.prop=freq.df$Freq/sum(freq.df$Freq)
  }

  san.hist.stop=proc.time()
  san.hist.time=(san.hist.stop-mv.hist.stop)[[3]]


  #if treatment within blocks with set block size, must maintain that block size


  # #if using delta>0, and number of bins is high enough
  # if((delta>0)&(num.bins.level>(2/delta))){ #use Bun et al. 2016
  #   threshold=((2*base::log(2/delta))/(epsilon/n))+(1/n)
  # }else{ #pure-DP and low number of bins don't use a threshold
  #   threshold=0
  # }

  idx.to.sample=seq(1,nrow(freq.df)+length(san.prop.zero.to.add))
  #sample rows of mv hist with probabilities equal to norm.san with replacement
  # get sample of size equal to number of rows of data.
  row.sample<-sample(idx.to.sample,size=nrow(data),
                     replace=TRUE,prob=c(freq.df$san.prop,san.prop.zero.to.add))
  #synthetic data has values from the sample of the histogram
  # (remove frequency, san.prop columns)
  synth.data<-freq.df[row.sample[row.sample<=nrow(freq.df)],colnames(freq.df)%in%colnames(data)]
  zero.to.sample=row.sample[row.sample>nrow(freq.df)]
  if(length(zero.to.sample)>0){
    zero.to.sample=as.numeric(as.factor(zero.to.sample)) #reindex 1,...,nsample
    nsample=length(unique(zero.to.sample))
    create.zero.rows=unrealized_sampler(realized.df=freq.df[,colnames(freq.df)%in%colnames(data)],
                                        levels.list = levels.list,n.realized=nrow(freq.df),
                                        nsample=nsample,orig.nreal = nrow(freq.df))
    unrealized.rows=create.zero.rows[zero.to.sample,]
    synth.data= dplyr::bind_rows(synth.data,unrealized.rows)

    rownames(synth.data)=NULL
    stopifnot(nrow(synth.data)==nrow(data))
  }

  #force columns to be factors or numeric (add variation if needed)
  synth.data[,which.cat]=lapply(synth.data[,which.cat],
                                function(col)as.factor(as.character(col)))
  if(add.cont.variation==TRUE){# if adding uniform variation for continuous values
    cont.v=colnames(synth.data[,which.cont,drop=F])
    if(length(cont.v)==0){
      synth.data=synth.data
    }else if(length(cont.v)==1){
      synth.data[,cont.v]=synth_continuous_variation(synth.data[,cont.v],levels.list[[cont.v]])
    }else{
      # count.nas=sapply(cont.v,function(x)sum(is.na(synth.data[,cont.v])))
      # count.vals=sapply(cont.v,function(x)length(unique(as.character(synth.data[,cont.v])))<=1)
      # if(sum(count.nas,na.rm=T)>0){
      #   warning(paste(cont.v[count.nas>0],"has",count.nas," NA values",collapse="\n"))
      # }
      # if((sum(count.vals,na.rm=T)>0)|(sum(is.na(count.vals))>0)){
      #   warning(paste("some continuous variable only has one value?:",paste0(cont.v[count.vals],collapse=", "),
      #                 " or something weird with NAs ",sum(is.na(count.vals))))
      # }
      synth.data[,which.cont]=lapply(cont.v,function(cv)synth_continuous_variation(synth.data[,cv],levels.list[[cv]]))
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
        message(warn.mess)
        treatment.colname="treatment"
      }
      for(i in c(1,2)){
        #warning(paste("inside for loop i=",i))
        #        if((length(treatment.colname)==length(conditions))&(length(treatment.colname)>1)){
        #          warning("treatment.colname has same length as conditions")

        synth.data=treatment_assign(synth.data=synth.data,
                                    assign.type=assign.type,
                                    treatment.colname=conditions[i],
                                    blocks=unlist(blocks.ls[[i]]))
        #,conditions=c("1","0"))#,#clusters=clusters, ...)
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

    }else{
      #warning("factorial==FALSE")
      synth.data=treatment_assign(synth.data=synth.data,
                                  assign.type=assign.type,
                                  treatment.colname=treatment.colname,
                                  blocks=blocks,conditions=conditions)#,clusters=clusters,...)
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
#' @param lvls is vector of all possible midpoints (levels of the factor)
#' @return a numeric vector of synthetic continuous values
#' @noRd
#' @keywords internal
#'
#' @importFrom stats runif
synth_continuous_variation<-function(cat.var,lvls){
  n.rw=length(cat.var)
  midpoints=as.numeric(as.character(lvls))
  n.levels=length(midpoints)
  #half the interval length is 1/2 difference between midpoints
  half.width=(midpoints[2]-midpoints[1])/(2)
  variation=stats::runif(n.rw,-half.width,half.width) #uniform rv
  cont.var=as.numeric(as.character(cat.var))+variation
  return(as.numeric(cont.var))
}


#' Sample from unrealized combinations of variable values
#'
#' @description When sampling from a perturbed histogram, some rows will be sampled that
#'      did not appear in the original data set. This sampler samples without enumerating
#'      all possible combinations of column values.
#'
#' @param realized.df dataset of unique rows observed in the original data
#' @param levels.list a named list where the name is each column of the original data
#'    and the values are the possible levels for that column
#' @param nsample number of unique unrealized rows to be generated
#' @param n.realized default=NA (used for recursion)
#' @param current.iter default=0 (used for recursion)
#' @param orig.nreal default=NA (used for recursion)
#' @param max.iter default=1000, maximum number of iterations to try to sample the unrealized rows
#' @return a data.frame with nsample unique unrealized rows
#' @importFrom dplyr distinct bind_rows
#'
#'
#' @family syntheticData
#' @export
#'
unrealized_sampler=function(realized.df,levels.list,nsample,n.realized=NA,current.iter=0,orig.nreal=NA,max.iter=1000){
  current.iter=current.iter+1
  #rownames(realized.df)=NULL
  if(is.na(n.realized)==TRUE){
    n.realized=nrow(realized.df)
  }
  if(is.data.frame(realized.df)==FALSE){
    realized.df=as.data.frame(realized.df)
  }
  if(is.na(orig.nreal)==TRUE){
    orig.nreal=n.realized
  }
  if(current.iter>max.iter+1){
    stop(paste("Iterations in unrealized sampler exceed the max iterations: ",max.iter))
  }
  if(nrow(realized.df)>=n.realized+nsample){
    base::message(paste("Iterations to sample unrealized combinations of variables:",current.iter))
    return(realized.df[-seq(1,orig.nreal),])
  }else{
    if(nsample==1){
      possible=data.frame(matrix(nrow=1,ncol=ncol(realized.df)))
      colnames(possible)=colnames(realized.df)
      possible[1,]=unname(sapply(levels.list,function(x)sample(x,nsample,replace=T)))
    }else{
      possible=as.data.frame(sapply(levels.list,function(x)sample(x,nsample,replace=T)))
    }
    #rownames(possible)=NULL
    #colnames(possible)=colnames(realized.df)
    combo=dplyr::distinct(dplyr::bind_rows(realized.df,possible))
    new.n.realized=nrow(combo)
    new.nsample=nsample-(new.n.realized-n.realized)
    return(unrealized_sampler(realized.df = combo,levels.list=levels.list,
                              nsample=new.nsample,n.realized=new.n.realized,
                              orig.nreal=orig.nreal,current.iter=current.iter,max.iter=max.iter))
  }
}

