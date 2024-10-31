#' Function to Produce Two Diagnostic Plots for Linear Regression
#' @param reg.model is the fitted \code{\link{glm}} or \code{\link{lm}} object
#' @param data data.frame to
### Function to Check Regression Assumptions ###
reg_assumptions=function(reg.model,data,family="gaussian",
                         response.var=NULL,pt.sz=2.5,ln.sz=1.5,bs.sz=7,
                         response.var.name=NULL,
                         include.partial.resid.plot==FALSE,
                         ...){
  mod.fit=glm(reg.model,data=data,family=family,... )
  if(is.null(response.var)==TRUE){
    response.var=names(attr(tempmod$terms,"dataClasses"))[1]
  }
  if(is.null(response.var.name)==NULL){
    response.var.name=response.var
  }
  data.temp=dplyr::bind_cols(data.frame("fitted"=mod.fit$fitted,
                                        "residuals"=mod.fit$residuals,
                                        "rstandard"=stats::rstandard(mod.fit),
                                        "zeros"=0),
                             data)
  if(family=="binomial"){
    data.temp$logit.fit=boot::logit(mod.fit$fitted)
    data.temp$residuals=stats::rstandard(mod.fit)
    plot.fit.resid=ggplot2::ggplot(data.temp,ggplot2::aes(x=logit.fit,y=residuals,col=treatment))+
      ggplot2::geom_point(size=pt.sz)+
      ggplot2::geom_smooth(se=F,col="red",linewidth=ln.sz)+
      ggplot2::geom_line(ggplot2::aes(x=fitted,y=zeros),col="black",linetype=2,size=ln.sz)+
      ggplot2::labs(x=resid.fit.x.lab,y="Standardized Residuals")+
      ggplot2::ggtitle(paste0(response.var.name,": Residual vs. Fitted"))+
      ggplot2::theme_minimal(base_size=bs.sz)
    out.plot=plot.fit.resid
  }else{
    data.temp$residuals=mod.fit$residuals
    data.temp$fitted=mod.fit$fitted
    plot.fit.resid=ggplot2::ggplot(data.temp,ggplot2::aes(x=fitted,y=residuals))+
      ggplot2::geom_point(size=pt.sz)+
      ggplot2::geom_smooth(se=F,col="red",linewidth=ln.sz)+
      ggplot2::geom_line(ggplot2::aes(x=fitted,y=zeros),col="black",linetype=2,size=ln.sz)+
      ggplot2::labs(x=resid.fit.x.lab,y="Residuals")+
      ggplot2::ggtitle(paste0(response.var.name,": Residual vs. Fitted"))+
      ggplot2::theme_minimal(base_size=bs.sz)
  }
  plot.fit.resid=ggplot2::ggplot(data.temp,ggplot2::aes(x=fitted,y=residuals,col=resid.fit.col))+
    ggplot2::geom_point(size=pt.sz)+
    ggplot2::geom_smooth(se=F,col="red",linewidth=ln.sz)+
    ggplot2::geom_line(ggplot2::aes(x=fitted,y=zeros),col="black",linetype=2,size=ln.sz)+
    ggplot2::labs(x=resid.fit.x.lab,y="StResiduals")+
    ggplot2::ggtitle(paste0(response.var.name,": Residual vs. Fitted"))+theme.add

  if(family=="gaussian"){
  data.temp$norm=stats::rnorm(nrow(data.temp))
  plot.qqnorm=ggplot2::ggplot(data=data.temp,ggplot2::aes(sample=residuals))+
    ggplot2::geom_qq(size=pt.sz,color="black")+
    ggplot2::geom_qq_line(color="red",linewidth=ln.sz)+
    ggplot2::labs(x="Theoretical Quantiles",y="Residuals Quantiles")+
    ggplot2::ggtitle(paste0(response.var.name,": Normal Q-Q Plot"))+
    ggplot2::theme_minimal(base_size=bs.sz)
    out.plot=cowplot::plot_grid(plot.fit.resid,plot.qqnorm,ncol=2)
  }
  if(include.partial.resid.plot==TRUE){
    partial.resid=car::crPlots(mod.fit)
    out.plot=list("simp.plots"=out.plot,"partial.residual"=partial.resid)
  }
  return(out.plot)
}
