#'Generate a frequency data.frame for multivariate histogram of mixed variable types
#'
#' @param data a data.frame to generate (epsilon,delta) differentially private version of
#' @param continuous.vars list of strings for the names, vector of indexes, or vector
#'    of logical for columns of data with continuous variables
#' @param num.bin a numeric value or vector of numeric values for the number of bins
#'    for each continuous variable. If =NA use bin.param to determine the number of bins
#' @param bin.param a numeric value in (0,1) to determine number of bins for continuous
#'    variables to be used if num.bins==NA.
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
                                 num.bin=NULL,bin.param=NA,
                                 which.cont.out=FALSE,return.time=TRUE){
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
    cont.data=data.frame(data[,continuous.vars]) #subset data to be continuous variables
    if(base::ncol(cont.data)==0){
      warning("No continuous variables detected.")
    }
    #number of continuous variables
    num.continuous=ncol(cont.data)
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
    which.cont=base::colnames(data)%in%base::colnames(cont.data) #logical if continuous variable
    data[,which.cont]=base::sapply(base::seq(1,num.continuous),
                                   function(i)continuous_bins(cont.data[,i],num.bin[i]))
  }# end if there are continuous variables
  mv.histogram=dplyr::count(data,data[seq(1,nrow(data)),],name="Freq")

  #if which.cont.out=TRUE return mv.histogram and which.cont,else only mv.histogram
  if(which.cont.out==TRUE){
    return(base::list(mv.histogram,which.cont))
  }else{
    return(mv.histogram)
  }
}

