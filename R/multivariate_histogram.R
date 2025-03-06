#'Generate a frequency data.frame for multivariate histogram of mixed variable types
#'
#' @param data a data.frame to generate (epsilon,delta) differentially private version of
#' @param continuous.vars list of strings for the names, vector of indexes, or vector
#'    of logical for columns of data with continuous variables
#' @param num.bin a numeric value or vector of numeric values for the number of bins
#'    for each continuous variable. If =NA use bin.param to determine the number of bins
#' @param bin.param a numeric value in (0,1) to determine number of bins for continuous
#'    variables to be used if num.bins==NA.
#' @param continuous.limits a list of upper and lower bounds for each continuous variable.
#'    If there are continuous variables, then a limit must be supplied for each.
#' @param standardized.cont a vector of column names to transform to a standardized
#'    form $(x_i-\bar{x})/std.dev.(x)$. If \code{NULL}, then no variable is transformed.
#' @param std.limits a numeric value for the number of standard deviations away from the mean
#'    to set the transformed standardized column limits to. Default is 15.
#' @param which.cont.out a logical if TRUE, return which.count and histogram data.frame.
#'    Default FALSE only returns histogram data.frame
#' @param return.time a logical to indicate if the computation times should be
#'    returned with the synthetic data. Default is \code{FALSE}.
#' @return if which.count.out=FALSE, returns a data frame with frequencies of each combination of categorical and discretized continuous variables. If which.cont.out=TRUE and there are continuous variables, returns a list with the first element is the mv histogram and the second is the which.cont variable.
#'
#' @importFrom dplyr bind_rows mutate_all
#' @importFrom parallel mclapply
#'@examples
#'multivariate_histogram(mtcars[,2:4],continuous.vars=c("disp","hp"),num.bin=3)
#'
#'@family syntheticData
#'@export
multivariate_histogram<-function(data,continuous.vars=NULL,
                                 num.bin=NULL,bin.param=NA,continuous.limits=NULL,
                                 which.cont.out=FALSE,levels.out=FALSE,
                                 return.time=TRUE,std.limits=15,standardize.cont=NULL){#,check.cont=F){
  ### Check Inputs ###
  stopifnot(base::is.data.frame(data)) #check data input

  ### Multivariate Number of Bins Vector ####
  ## each entry is the number of bins for the corresponding column in data

  #get number of continuous variables
  if(base::is.null(continuous.vars)){ #if continuous.vars=NULL
    num.continuous=0 #number continuous variables=0, no need num.bin or bin.param
    #mx.n.combos=prod(sapply(colnames(data),function(x)length(table(data[,x],useNA = "ifany"))))
    #return(base::data.frame(base::table(data)))
  }else{ #if continuous.vars is supplied
    cont.data=data[,continuous.vars,drop=F] #subset data to be continuous variables
    if(base::ncol(cont.data)==0){
      num.continuous=0
      message(paste0("No continuous variables detected. Dimension of cont.data is...",dim(cont.data),"... and length of continuous.vars is...",length(continuous.vars)))

    }

    #number of continuous variables
    num.continuous=base::ncol(cont.data)
    if(is.null(standardize.cont)==FALSE){
      std.not.in.cont=NULL
      if(sum(colnames(cont.data)%in% c(standardize.cont))!=length(standardize.cont)){
        std.not.in.cont=standardize.cont[!(standardize.cont%iN%colnames(cont.data))]
        warning(paste("Some columns in standardize.cont are not in continuous.vars and/or the dataset. Combining them:",
                      paste0(std.not.in.cont,collapse=", ")))
        cont.data=rbind(cont.data,data[,std.not.in.cont,drop=F])
      }
      new.num.continuous=ncol(cont.data)
      std.idx=standardize.cont%in% colnames(cont.data)
      cont.data[,standardize.cont[std.idx]]=
        lapply(cont.data[,standardize.cont[std.idx]],
               function(x)(x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T)))
      ncontlim=base::length(continuous.limits)
      if((ncontlim<2)&(num.continuous>1)){
        message("Only one continuous limit supplied. It will be used for all the continuous variables. That are not standardized.")
        continuous.limits=c(base::rep(continuous.limits,num.continuous),rep(list(-std.limits,std.limits),sum(std.idx)))
      }else if(ncontlim==num.continuous){
        continuous.limits=c(continuous.limits,rep(list(-std.limits,std.limits),new.num.continuous-num.continuous))
      }else if(ncontlim==new.num.continuous){
        message("countinus.limits provides a bound for each continuous column. std.limits input is ignored.")
      }else{
        stop("continuous.limits length with continuous.vars and standardize.cont inputs do not lead to a clear continuous.limits value for each continuous variable.")
      }
      if(is.null(names(continuous.limits))==TRUE){ #if no names, name them after the continuous columns
        names(continuous.limits)=colnames(cont.data)
      }else{ #otherwise put continuous.limits in the same order as cont.data
        names(continuous.limits)[names(continuous.limits)==""]=standardize.cont[!(standardize.cont%in%names(continuous.limits))]
        continuous.limits=continuous.limits[colnames(cont.data)]
      }
      num.continuous=new.num.continuous
    }else{ #no standardization
    if((base::length(continuous.limits)<2)&(num.continuous>1)){
      message("Only one continuous limit supplied. It will be used for all the continuous variables.")
      continuous.limits=base::rep(continuous.limits,num.continuous)
    }
    stopifnot(length(continuous.limits)==num.continuous)
    if(is.null(names(continuous.limits))==TRUE){ #if no names, name them after the continuous columns
      names(continuous.limits)=colnames(cont.data)
    }else{ #otherwise put continuous.limits in the same order as cont.data
      continuous.limits=continuous.limits[colnames(cont.data)]
    }
    }

  if(num.continuous>0){ #if there are continuous variables
    if(base::sum(base::sapply(cont.data,base::is.numeric))!=num.continuous){ #check columns are numeric
      stop(
        base::paste("Continuous variables selected are not numeric values. Columns Selected:",
                    base::paste(base::colnames(cont.data),collapse=", ")))
    } #end if continuous variables are not numeric

    #if there are continuous variables, either num.bin or bin.param must be numeric
    stopifnot(base::is.numeric(num.bin)|base::is.numeric(bin.param))
    if(base::is.numeric(num.bin)==T){ #num.bin must be a single value or a value for each column
      stopifnot(length(num.bin)==1|base::length(num.bin)==num.continuous)

      if(base::length(num.bin)==1){ #if num.bin is single value
        num.bin=base::rep(num.bin,num.continuous)
      } #end if num.bin single value

    }else{ #if num.bin not supplied bin.param must be a single value between 0 and 1
      stopifnot(bin.param>0&bin.param<1)
      num.bin=base::rep(base::ceiling((base::nrow(cont.data)^bin.param)),num.continuous)
    }
    #if(check.cont==T){
    #cont.gr.bin=sapply(seq(1,num.continuous),function(i)length(unique(cont.data[,i]))>num.bin[i])
    #if(sum(!cont.gr.bin)>0){
    #  cont.data=cont.data[,cont.gr.bin,drop=F]
    #  num.continuous=ncol(cont.data)
    #}
    }
    if(num.continuous>0){
    which.cont=base::colnames(data)%in%base::colnames(cont.data) #logical if continuous variable
    data[,which.cont]=
      base::lapply(base::seq(1,ncol(cont.data)),
                                   function(i)continuous_bins(cont.data[,i],num.bins=num.bin[i],cont.limit=continuous.limits[[i]]))

    }else{
      which.cont=rep(F,ncol(data))
    }
  }# end if there are continuous variables
  data[,sapply(data,function(x)!is.factor(x))]=lapply(data[,sapply(data,function(x)!is.factor(x))],as.factor)
  mv.histogram=dplyr::count(data,data[seq(1,nrow(data)),],name="Freq")

  #if which.cont.out=TRUE return mv.histogram and which.cont,
  # if levels.out=TRUE include levels.out in output
  # else only mv.histogram returned
  if(which.cont.out==TRUE){
    if(levels.out==TRUE){
      levels.list=list(lapply(data,levels))
      out=list(mv.histogram,which.cont,levels.list)
      names(out)=c("mv.histogram","which.cont","levels.list")
    }else{
      out=list(mv.histogram,which.cont)
      names(out)=c("mv.histogram","which.cont")
    }
  }else{
    if(levels.out==TRUE){
      levels.list=list(lapply(data,levels))
      out=list(mv.histogram,levels.list)
      names(out)=c("mv.histogram","levels.list")
    }else{
      out=mv.histogram
    }
  }

  return(out)
}

