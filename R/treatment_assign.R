#' Assign Treatments Using randomizr Package Functions
#' @importFrom randomizr simple_ra
#' @importFrom randomizr complete_ra
#' @importFrom randomizr block_ra
#' @importFrom randomizr cluster_ra
#' @importFrom randomizr block_and_cluster_ra
#' @param synth.data a data.frame with the covariates including any blocking or
#'    clustering variables to be used in assignment
#' @param assign.type a string taking values "simple","complete","block","cluster",
#'    or "block_and_cluster" to specify what kind of random assignment should be used.
#'    See Details.
#' @param treatment.colname a string to rename the treatment column. Default is "treatment"
#' @param blocks is a vector of indices, column names, or logical to specify
#'    what columns of \code{synth.data} to use as blocking variables in random
#'    treatment assignment. If \code{assign.treatment} takes value "block" or
#'    "block_and_cluster" a value for \code{blocks} must be supplied.
#' @param clusters is a vector of indices, column names, or logical to specify
#'    what columns of \code{synth.data} to use as clustering variables in random
#'    treatment assignment. If \code{assign.treatment} takes value "cluster" or
#'    "block_and_cluster" a value for \code{clusters} must be supplied.
#' @param num.arms is the number of treatment arms. If \code{num.arms=NULL}, the value
#'    may be determined by additional arguments passed to \code{treatment_assignment}.
#'    If not determined by these arguments, then there will be two treatment conditions.
#' @param ... Additional arguments will be passed to the corresponding \pkg{randomizr}
#'    function, where the function is determined by \code{assign.type}.
#' @details Additional inputs in the function are supplied to the respective \pkg{randomizr} function.
#' The randomized treatment assignment come from functions in \pkg{randomizr}.
#' If \code{assign.type}="simple", \code{\link{simple_ra}} is used.
#' If \code{assign.type}="complete", \code{\link{complete_ra}} is used.
#' If \code{assign.type}="block", \code{\link{block_ra}}
#' is used and \code{blocks} must be supplied in the input.
#' If \code{assign.type}="cluster", \code{\link{cluster_ra}}
#' is used and \code{clusters} must be supplied in the input.
#' If \code{assign.type}="block_and_cluster",
#' \code{\link{block_and_cluster_ra}}
#' is used and \code{blocks} and \code{clusters} must be supplied in the input.
#'
#'
#' @returns a data.frame that adds a treatment column to synth.data.
#'
#' @examples
#' set.seed(1)
#' data=data.frame("block.cov"=rep(c("A","A","B","B","B","C","C"),20),
#'                 "covariate"=rnorm(140,48,3))
#' treatment_assign(data,assign.type="block",num.arms=3,blocks="block.cov")
#'
#' @family syntheticData
#' @export
treatment_assign<-function(synth.data,
                           assign.type,
                           num.arms=NULL,
                           treatment.colname="treatment",
                           blocks=NULL,clusters=NULL,...){
  stopifnot(assign.type %in% c("simple","complete","block","cluster","block_and_cluster"))

  if(base::is.data.frame(synth.data)==TRUE){
    #if blocks used in treatment assignment, must have block input
    if(base::grepl("block",assign.type)==T){
      stopifnot(!base::is.null(blocks))
      if(base::is.logical(blocks)&&base::sum(blocks)==0){
        stop("No blocking variables provided.")
      }
      block.df=synth.data[,blocks,drop=F]
      if(base::ncol(block.df)==1){ #if only one blocking variable
        block.combine=unlist(c(block.df))
      }else{ #otherwise combine into one variable
        block.combine=base::apply(block.df,1,function(rw)base::paste(rw,collapse="_"))
      }
    }

    #repeat block checks above on clusters if clusters used in treatment assignment
    if(base::grepl("cluster",assign.type)==T){
      stopifnot(!is.na(clusters))
      cluster.df=synth.data[,clusters,drop=F]
      if(base::ncol(cluster.df)==1){
        cluster.combine=cluster.df
      }else{
        cluster.combine=base::apply(cluster.df,1,function(rw)paste(rw,collapse="_"))
      }
    }
    num.rows=base::nrow(synth.data)
  }else{
    num.rows=base::length(synth.data)
    synth.data=base::data.frame(synth.data,"treat"=base::rep(NA,num.rows))
  }

  # assign treatments using randomizr functions depending on assign.type
  if(assign.type=="simple"){ #simple assignment
    synth.data$treat=randomizr::simple_ra(N=num.rows,num_arms=num.arms,...)
  } else if(assign.type=="complete"){ #complete assignment
    synth.data$treat=randomizr::complete_ra(N=num.rows,num_arms=num.arms,...)
  }else if(assign.type=="block"){ #block assignment
    synth.data$treat=randomizr::block_ra(blocks=block.combine,num_arms=num.arms,...)
  }else if(assign.type=="cluster"){ #cluster assignment
    synth.data$treat=randomizr::cluster_ra(clusters=cluster.combine,num_arms=num.arms,...)
  }else if(assign.type=="block_and_cluster"){ #block and cluster assignment
    synth.data$treat=randomizr::block_and_cluster_ra(blocks=block.combine,clusters=cluster.combine,num_arms=num.arms,...)
  }else{ #OH NO!
    stop(base::paste('ERROR: Something is wrong with assign.type.',
                     'It should be "simple","complete","block","cluster", or "block_and_cluster".'))
  }


  #rename the column
  base::colnames(synth.data)=c(base::colnames(synth.data)[base::seq(1,base::ncol(synth.data)-1)],
                               treatment.colname)
  return(synth.data)
}
