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
                                 return.time=TRUE,std.limits=5,standardize.cont=NULL){#,check.cont=F){
  ### Check Inputs ###
  stopifnot(base::is.data.frame(data)) #check data input

  if((is.null(continuous.limits)==F)&&(is.list(continuous.limits)==FALSE)){
    continuous.limits=as.list(continuous.limits)
  }

  stopifnot((is.null(num.bin)==FALSE)|(is.null(bin.param)==FALSE))
  if(is.null(num.bin)==TRUE){
    num.bin=nrow(data)^bin.param
  }


  ### Multivariate Number of Bins Vector ####
  ## each entry is the number of bins for the corresponding column in data

  orig.colnms=colnames(data)
  #get number of continuous variables
  if(base::is.null(continuous.vars)){ #if continuous.vars=NULL
    num.continuous=0 #number continuous variables=0, no need num.bin or bin.param
    cont.data=NULL
  }else{ #if continuous.vars is supplied
    continuous.vars=continuous.vars[continuous.vars%in%colnames(data)]
    cont.data=data[,continuous.vars,drop=F] #subset data to be continuous variables
    if(base::ncol(cont.data)==0){
      continuous.vars=NULL
      num.continuous=0
      message(paste0("No continuous variables detected. Dimension of cont.data is...",dim(cont.data),"... and length of continuous.vars is...",length(continuous.vars)))
    }

    #number of continuous variables
    orig.contdata.nm=colnames(cont.data)
    num.continuous=base::ncol(cont.data)
    if(is.null(standardize.cont)==FALSE){ #variables will be standardized
      std.not.in.cont=NULL
      standardize.cont=standardize.cont[standardize.cont%in%colnames(data)]

      #if there are variables in standardize.cont that are not in cont.data. add them and adjust num.continuous
      if(sum(colnames(cont.data)%in% c(standardize.cont))!=length(standardize.cont)){
        std.not.in.cont=standardize.cont[!(standardize.cont%in%colnames(cont.data))]
        warning(paste("Some columns in standardize.cont are not in continuous.vars and/or the dataset. Combining them:",
                      paste0(std.not.in.cont,collapse=", ")))
        cont.data=cbind(cont.data,data[,std.not.in.cont,drop=F])
      }
      #reorder to have standardized at the end
      new.cont.nms=c(continuous.vars[!(continuous.vars%in%standardize.cont)],standardize.cont)
      cont.data=cont.data[,c(new.cont.nms),drop=F]
      new.num.continuous=ncol(cont.data)

      if(length(standardize.cont)==1){ #if only 1 standardize
        std.idx=(colnames(cont.data)==standardize.cont)
        pre.std.data=unlist(cont.data[,std.idx])
        cont.data[,std.idx]=(pre.std.data-mean(pre.std.data,na.rm=T))/sqrt(var(pre.std.data,na.rm=T))
      }else if(length(standardize.cont)>1){ #more than 1 standardize
        std.idx=(colnames(cont.data)%in%standardize.cont)
        cont.data[,std.idx]=
          lapply(cont.data[,std.idx],
                 function(x)(x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T)))
      }else{ #no standardized
        message("standardize.cont where not variables in data. Nothing to standardize")
      }


      ##deal with cont.limits now
      if(is.null(names(continuous.limits))==FALSE){
        names.order=colnames(cont.data)[colnames(cont.data)%in%names(continuous.limits)]
        continuous.limits=continuous.limits[names.order]
      }
      ncontlim=base::length(continuous.limits)
      cont.not.std=continuous.vars[!(continuous.vars%in%c(standardize.cont))]
      num.cont.not.std=length(cont.not.std)

      #more limits supplied than continuous variables
      if(ncontlim>=new.num.continuous){
          message("countinus.limits provides a bound for each continuous column. std.limits input is ignored.")
          continuous.limits=continuous.limits[seq(1,new.num.continuous)]


      #Cases:
      # there are continuous variables that are not standardized:
      ##  only 1 limit supplied
      ##  only limits supplied for continuous that are not standardized:
      ##  only limits supplied for variables named in continuous variables
      ##  limits supplied for more than the not standardized variables

      }else if((num.cont.not.std>0)){ #there are variables that are not standardized
          if((ncontlim==1))){
            message("Only one continuous limit supplied. It will be used for all the continuous variables. That are not standardized.")
            continuous.limits=c(base::rep(continuous.limits,num.cont.not.std),rep(list(-std.limits,std.limits),new.num.continuous-num.cont.not.std))
          }else if(ncontlim==num.cont.not.std){ #add standardized limits to cont.limits
            continuous.limits=c(continuous.limits,rep(list(-std.limits,std.limits),new.num.continuous-num.cont.not.std))
          }else if(ncontlim==num.continuous){
            continuous.limits=c(continuous.limits,rep(list(-std.limits,std.limits),new.num.continuous-num.continous))
          }else if(ncontlim>num.cont.not.std){
            continuous.limits=c(continuous.limits,rep(list(-std.limits,std.limits),new.num.continuous-ncontlim))
          }else{
          stop("continuous.limits length with continuous.vars and standardize.cont inputs do not lead to a clear continuous.limits value for each continuous variable.")
        }
      }else{ #there are standardized and no variables that are not standardized
        if(ncontlim==0){
          message("No continuous limits supplied. Used std.limits for all variables.")
          continuous.limits=rep(list(-std.limits,std.limits),new.num.continuous)
        }else if(ncontlim==num.continuous){
          continuous.limits=c(continuous.limits,rep(list(-std.limits,std.limits),new.num.continuous-num.continous))
        }else{
          stop("continuous.limits length with continuous.vars and standardize.cont inputs do not lead to a clear continuous.limits value for each continuous variable.")
        }

      }
      num.continuous=new.num.continuous
      names(continuous.limits)=colnames(cont.data)
      }else{ #no standardization
        if(is.null(names(continuous.limits))==FALSE){
          names.order=colnames(cont.data)[colnames(cont.data)%in%names(continuous.limits)]
          continuous.limits=continuous.limits[names.order]
        }
        if((base::length(continuous.limits)<2)&(num.continuous>1)){
          message("Only one continuous limit supplied. It will be used for all the continuous variables.")
          continuous.limits=base::rep(continuous.limits,num.continuous)
        }else if(length(continuous.limits)>=num.continuous){
          continuous.limits=continuous.limits[seq(1,num.continuous)]
        }else{
          stop("Continuous limits not supplied for each continuous variable")
        }
        names(continuous.limits)=colnames(cont.data)
      }# end if standardized variables else no standardized variables


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
    print(num.continuous)
    print(colnames(cont.data))
    print(continuous.vars)
    print(standardize.cont)

    if(num.continuous>1){
      if(length(num.bin)==1){
        num.bin=rep(num.bin,num.continuous)
      }
    #which.cont=which(base::colnames(data)%in%base::colnames(cont.data)) #logical if continuous variable
    data[,colnames(cont.data)]=
      base::lapply(base::seq(1,ncol(cont.data)),
                                   function(i)continuous_bins(cont.data[,i],num.bins=num.bin[i],cont.limit=continuous.limits[[i]]))

    which.cont=colnames(data)%in%colnames(cont.data)
    }else if(num.continuous==1){
      data[,colnames(cont.data)]=continuous_bins(cont.data,num.bins=num.bin,cont.limit=unlist(continuous.limits))
      which.cont=colnames(data)%in%colnames(cont.data)
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

